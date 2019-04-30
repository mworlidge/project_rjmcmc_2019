## RANDOM WALK METROPOLIS HASTINGS ALGORITHM WITH SPACE SHUTTLE DATA

x <- c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
y <- c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)

model <- glm(y ~ x + I(x^2), family = binomial(logit))
summary(model)

a.mle <- as.numeric(model$coefficients[1])
b.mle <- as.numeric(model$coefficients[2])
c.mle <- as.numeric(model$coefficients[3])

var.a.mle <- summary(model)$cov.scaled[1, 1]
var.b.mle <- summary(model)$cov.scaled[2, 2]
var.c.mle <- summary(model)$cov.scaled[3, 3]

# LIKELIHOOD FUNCTION

likelihoodfunc <- function(theta){
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  p <- 1 - 1 /(1 + exp(a + (b*x) + c*(x^2)))
  likelihood <- exp(sum(dbinom(y,size=1,prob=p,log=TRUE)))
  return(likelihood)
}

# POSTERIOR
dpost<- function(theta){
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  p <- 1 - 1 /(1 + exp(a + (b*x) + c*(x^2)))
  likelihood <- exp(sum(dbinom(y,size=1,prob=p,log=TRUE)))
   if (likelihood == 0){
     likelihood <- .Machine$double.xmin
   }

  priora <- dnorm(a,80, 40)                  # Don't use mle, not bayesian.
  priorb <- dnorm(b,-2, 2)
  priorc <- dnorm(c,0.01, 0.01)
  
  return(likelihood*priora*priorb*priorc)
}

# PROPOSAL (with random walk)

rprop <- function(theta0){
  propa <- rnorm(1, 0, 7)      #
  propb <- rnorm(1, 0, 1)
  propc <- rnorm(1, 0, 0.001)
  proposal <- c(theta0[1]+propa,theta0[2]+propb, theta0[3]+propc)
  if (likelihoodfunc(proposal) == 0){
    return(rprop(theta0))
  }
  return(proposal)
}

# METROPOLIS-HASTING RANDOM WALK ALGORTHIM

mh <- function(x0,f,rprop,N,B){
  x <- matrix(NA, N+B, length(x0))
  fx <- rep(NA, N+B)
  x[1,] <- x0
  fx[1] <- f(x0)
  ct <- 0
  
  for(i in 2:(N+B)){
    u <- rprop(x[i-1,])
    fu <- f(u)
    r <- log(fu) - log(fx[i-1])
    R <- min(exp(r),1)
    if(runif(1) <= R){
      
      ct <- ct + 1
      x[i,] <- u
      fx[i] <- fu
      
    } else{
      x[i,] <- x[i-1,]
      fx[i] <- fx[i-1]
    }
    
  }
  return(list(x=x[-(1:B),], fx=fx[-(1:B)], yes=ct/(N+B)))}  

# RUN THE ALGORITHM

N <- 1000000
B <- 50000
mh.out<- mh(c(76,-2.04,0.01), dpost,rprop,N,B)

par(mfrow=c(1,3))
mh.out$yes*100
a.mh <- mh.out$x[,1]
hist(a.mh, freq=FALSE, col="lightblue", border="black",xlab=expression(alpha), main="") 
plot(a.mh, type="l", col="lightblue", xlab="Iteration", ylab=expression(alpha), lwd=1.5) 
lines(1:N, cumsum(a.mh) / (1:N), col="red")

b.mh <- mh.out$x[,2]
hist(b.mh, freq=FALSE, col="lightgreen", border="black" ,xlab=expression(beta), main="") 
plot(b.mh, type="l", col="lightgreen", xlab="Iteration", ylab=expression(beta)) 
lines(1:N, cumsum(b.mh) / (1:N), col="red")

c.mh <- mh.out$x[,3]
hist(c.mh, freq=FALSE, col="lightpink", border="black" , xlab=expression(gamma), main="")
plot(c.mh, type="l", col="lightpink", xlab="Iteration", ylab=expression(gamma)) 
lines(1:N, cumsum(c.mh) / (1:N), col="red")

p65 <- 1 - 1 / (1 + exp(a.mh + b.mh*65 + c.mh*(65)^2)) 
hist(p65, freq=FALSE, col="darkgreen", border="white", xlab="p(65)", main="")

p75 <- 1 - 1 / (1 + exp(a.mh + b.mh*75 + c.mh*(75)^2))
hist(p75, freq=FALSE, col="darkgreen", border="white", xlab="p(75)", main="")

p45 <- 1 - 1 / (1 + exp(a.mh + b.mh*45 + c.mh*(45)^2))
hist(p45, freq=FALSE, col="darkgreen", border="white", xlab="p(45)", main="")

t <- seq(min(x)-5,max(x)+5, length.out=50)
mean_p_t <- rep(NA, 50)
func <- function(x, alpha, beta, c){
  for (i in (1:50)){
    p_t <- 1 - 1 / (1 + exp(alpha + beta*x[i]+c*x[i]^2))
    mean_p_t[i] <- mean(p_t)
  }
  return(mean_p_t)
}

mean <- func(t,a.mh, b.mh, c.mh)
mean
plot(x,y, main="Posterior expected value of probability of defect", xlab="Temperature", ylab="Probability")
lines(t, mean, type="l", col="red")
legend(53,0.23, c("Average posterior probability of defect"), c(2,2))

#install.packages("mcmcse")
library(mcmcse)
mcse.mat(x = chain, method = "bm", g = NULL)
mcerror_bm

# Point estiamte for alpha 
a.est <- mean(a.mh) # a.est   = 77.33597
a.se = 0.7809763412

#Point estimate for beta
b.est <- mean(b.mh) # = -2.003662
b.se = 0.0226471181

#Point estimate for gamm
c.est <- mean(c.mh) # 0.01250371
c.se <- 0.0001646241

# Point estimate for Pr(t|data)
pr65.est <- mean(p65) # 0.489195
pr75.est <- mean(p75) # 0.09458911
pr45.est <- mean(p45) # 0.9989155

install.packages("BAMBI")
library(BAMBI)
DIC(chain)
