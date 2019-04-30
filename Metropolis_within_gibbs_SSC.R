## NEW METHOD FOR BLOCKWISE GIBBS SAMPLER

## MODEL SET-UP
x <- c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
y <- c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)

model <- glm(y ~ x, family = binomial(logit))
summary(model)

a.mle <- as.numeric(model$coefficients[1])
b.mle <- as.numeric(model$coefficients[2])

var.a.mle <- summary(model)$cov.scaled[1, 1]
var.b.mle <- summary(model)$cov.scaled[2, 2]

# POSTERIOR FUNCTION
LogPostPDF <- function(theta){
  a <- theta[1]
  b <- theta[2]
  p <- 1 - 1 /(1 + exp(a + (b*x)))
  likelihood <- exp(sum(dbinom(y,size=1,prob=p,log=TRUE)))
  if (likelihood == 0){
    likelihood <- .Machine$double.xmin
  }
  prior.a <- dnorm(a,0,100)                  
  prior.b <- dnorm(b,0, 100)
  return(log(likelihood*prior.a*prior.b))
}

N <- 1000000
B <- 50000
mcmc_sample <- matrix(NA, N, 2)

alpha_now <- 10.0
beta_now <- -0.3

scl_alpha <- 0.1
scl_beta <- 0.1

cta <- 0
ctb <- 0

for(i in (1:N)) {
  ## ALPHA BLOCK
  
  alpha_prop <- alpha_now + scl_alpha*rnorm(1,0,10)
  beta_prop <- beta_now
  
  pdf_now <- exp(LogPostPDF(c(alpha_now, beta_now)))
  pdf_prop <- exp(LogPostPDF(c(alpha_prop, beta_prop)))
  
  acc_prob <- min(1.0, pdf_prop / pdf_now)
  
  u <- runif(1)
  
  if(u < acc_prob) {
    cta <- cta +1
    mcmc_sample[i,1] <- alpha_prop
    mcmc_sample[i,2] <- beta_now
    alpha_now <- alpha_prop
    beta_now <- beta_now
  } else {
    mcmc_sample[i,1] <- alpha_now
    mcmc_sample[i,2] <- beta_now
  }
  
  # BETA BLOCK
  
  alpha_prop <- alpha_now
  beta_prop <- beta_now + scl_beta*rnorm(1,0,10)
  
  pdf_now <- exp(LogPostPDF(c(alpha_now, beta_now)))
  pdf_prop <- exp(LogPostPDF(c(alpha_prop, beta_prop)))
  
  acc_prob <- min(1.0, pdf_prop/pdf_now)
  
  u <- runif(1)
  
  if(u < acc_prob) {
    ctb <- ctb + 1
    mcmc_sample[i,1] <- alpha_now
    mcmc_sample[i,2] <- beta_prop
    alpha_now <- alpha_now
    beta_now <- beta_prop
  } else {
    mcmc_sample[i,1] <- alpha_now
    mcmc_sample[i,2] <- beta_now
  }
}

gibbs.out <- mcmc_sample[-(1:B),]

###########################################################
# TRACE PLOTS AND HISTOGRAMS
###########################################################

a.gibbs <- gibbs.out[,1]
hist(a.gibbs, freq=FALSE, col="lightblue", border="black", xlab=expression(paste(alpha)),main="")
plot(a.gibbs, type="l", col="lightblue", xlab="Iteration", ylab=expression(paste(alpha)), main="")
lines(1:(N-B), cumsum(a.gibbs) / (1:(N-B)), col="red")

b.gibbs <- gibbs.out[,2]
hist(b.gibbs, freq=FALSE, col="lightgreen", border="black", xlab=expression(paste(beta)),main="")
plot(b.gibbs, type="l", col="lightgreen", xlab="Iteration", ylab=expression(paste(beta)), main="")
lines(1:(N-B), cumsum(b.gibbs) / (1:(N-B)), col="red")


p65 <- 1 - 1 / (1 + exp(a.gibbs + b.gibbs*65)) 
hist(p65, freq=FALSE, col="gray", border="white", xlab="p(65)", main="Posterior distribution for failure at 65 degrees")

p75 <- 1 - 1 / (1 + exp(a.gibbs + b.gibbs*75))
hist(p75, freq=FALSE, col="gray", border="white", xlab="p(75)", main="Posterior distribution for failure at 75 degrees")

p45 <- 1 - 1 / (1 + exp(a.gibbs + b.gibbs*45))
hist(p45, freq=FALSE, col="gray", border="white", xlab="p(45)", main="Posterior distribution for failure at 45 degrees")

plot(a.gibbs, b.gibbs, xlab=expression(paste(alpha)), ylab=expression(paste(beta)))

t <- seq(min(x)-5,max(x)+5, length.out=50)
mean_p_t <- rep(NA, 50)
mean_p_t
func <- function(x, alpha, beta){
  for (i in (1:50)){
    p_t <- 1 - 1 / (1 + exp(alpha + beta*x[i]))
    mean_p_t[i] <- mean(p_t)
  }
  return(mean_p_t)
}

logistic <- func(t,a.gibbs, b.gibbs)
mean
plot(x,y, main="Posterior expected value of probability of defect", xlab="Temperature", ylab="Probability")
lines(t, logistic, type="l", col="red")
legend(53,0.23, c("Average posterior probability of defect"), c(2,2))

meanA <- rep(NA, (N-B))
for (i in (1:(N-B))){
  meanA[i] <- mean(a.gibbs[1:i])
}
a <- cbind(rep(1:(N-B)),meanA)
plot(a, type="l", main="Convergence of the average", xlab="Iterations", ylab=expression(paste("cumulative average of  ", alpha)))

meanB <- rep(NA, (N-B))
for (i in (1:(N-B))){
  meanB[i] <- mean(b.gibbs[1:i])
}
b <- cbind(rep(1:(N-B)),meanB)
plot(b, type="l", main="Convergence of the average", xlab="Iterations", ylab=expression(paste("cumulative average of  ", beta)))
###########################################################

#install.packages("mcmcse")
library(mcmcse)
mcse.mat(x = gibbs.out, method = "bm", g = NULL)
mcerror_bm

# Point estiamte for alpha 
a.g.est <- mean(a.gibbs) # a.est   = 17.0386693
a.g.se = 0.268425726

#Point estimate for beta
b.g.est <- mean(b.gibbs) # = -0.2623476
b.g.se = 0.00395404

# Point estimate for Pr(t|data)
pr65.g.est <- mean(p65) # 0.489195
pr75.g.est <- mean(p75) # 0.09458911
pr45.g.est <- mean(p45) # 0.9989155

library(coda)
par(mfrow=c(1,1))
plot(gibbs.out, trace=TRUE, density=TRUE)

