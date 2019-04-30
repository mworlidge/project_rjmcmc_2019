
######### Poisson Negative-Binomial Example #########
## adapted from: https://github.com/andeek/reversibleMCMC/tree/master/coding


### Data ###
urlfile <- 'https://raw.githubusercontent.com/andeek/reversibleMCMC/master/data/soccer.csv'
full.data <- read.csv(urlfile)
football <- full.data$TotalGoals

####### Set up functions for RJMCMC Sampler ######
### Posterior Distribution Functions ###
pk1 <- function(y, theta, alphal, betal){
  sum(log(dpois(y, theta))) + log(dgamma(theta, shape=alphal, rate=betal))
}

pk2 <- function(y, theta, alphal, betal, alphak, betak){
  r <- 1/theta[2]
  sum(log(dnbinom(y, size=r, mu=theta[1]))) + log(dgamma(theta[1], 
                                                         shape=alphal, rate=betal)) + log(dgamma(theta[2], shape=alphak,
                                                                                                 rate=betak))
}

### Draw Lambda ###
qlam_k2<-function(lambda, y, kappa, alpha, beta){
  n <- length(y)
  log(lambda)*(sum(y) + alpha - 1) - (sum(y)+n/kappa)*log(1+kappa*lambda) -
    beta*lambda 
}

mh_lam_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters on prior
  # a, b are parameters of Gamma proposal distribution for lambda
  lambda_star<-rgamma(1, a, b)
  r <- qlam_k2(lambda_star, y, kappa, alpha, beta) - 
    qlam_k2(lambda, y, kappa, alpha, beta) + log(dgamma(lambda, a, b)) - 
    log(dgamma(lambda_star, a, b))
  return(ifelse(log(runif(1)) <= r, lambda_star, lambda))
}

sim_lambda <- function(lambda, y, alpha, beta, kappa=NULL, k, a, b){
  stopifnot(k %in% c(1,2))
  n <- length(y)
  if(k==1){
    lnew <- rgamma(1, shape=alpha + sum(y), rate=beta + n) 
  }else{
    lnew <- mh_lam_k2(lambda, y, kappa, alpha, beta, a, b)
  }
  return(lnew)
}

### Draw Kappa ###
qkap_k2 <- function(lambda, y, kappa, alpha, beta){
  n <- length(y)
  -n*lgamma(1/kappa) + sum(lgamma(1/kappa + y)) + 
    (sum(y) + alpha - 1)*log(kappa) - beta*kappa - 
    (sum(y) + n/kappa)*log(1+kappa*lambda)
}

mh_kap_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters on prior
  # a, b are parameters of Gamma proposal distribution for kappa
  kappa_star<-rgamma(1, a, b)
  r <- qkap_k2(lambda, y, kappa_star, alpha, beta) - 
    qkap_k2(lambda, y, kappa, alpha, beta) + log(dgamma(kappa, a, b)) - 
    log(dgamma(kappa_star, a, b))
  return(ifelse(log(runif(1)) <= r, kappa_star, kappa))
}

sim_kappa <- function(lambda, y, alpha, beta, kappa, a, b){
  n <- length(y)
  mh_kap_k2(lambda, y, kappa, alpha, beta, a, b)
}

################## RJMCMC SAMPLER ##################
rjmcmc_sampler <- function(y, lambda0, kappa0, k0=1, p=0.5, mu, sigma, 
                           alphal, betal, alphak, betak, al, bl, ak, bk, 
                           mc.iter = 1000){
  # y is the data
  # lambda0, kappa0, and k0 are initial values
  # p is prior probability of model 1, set to 0.5 as default
  # mu and sigma are fixed parameters of Normal in the between-model step 
  # alphal, betal and alphak, betak are Gamma hyperparameters for priors on 
  #   lambda and kappa, respectively
  # al, bl and ak, bk are Gamma parameters for the proposal distributions in
  #   their respective M-H algorithms
  
  # initialize data frame to save chains
  theta_save <- as.data.frame(matrix(NA, ncol=3, nrow=mc.iter+1))
  names(theta_save)<-c("lambda", "kappa", "k")
  accept_save <- as.data.frame(matrix(NA, ncol=2, nrow=mc.iter+1))
  names(accept_save)<-c("k1", "k2")
  
  # store initial values
  theta_save[1,] <- c(lambda0, kappa0, k0)
  lambda <- lambda0
  kappa <- kappa0
  k <- k0
  c_k1 <- 0
  c_k2 <- 0
  
  for(i in 1:mc.iter){
    # Reversible JUMP MCMC to move between models
    if(k==1){
      u <- rnorm(1, 0, sigma)
      theta <- lambda
      theta_new <- c(lambda, mu*exp(u))
      accept <- log(1-p) + pk2(y, theta_new, alphal, betal, alphak, betak) -
        log(p) - pk1(y, theta, alphal, betal) - log(dnorm(u, 0, sigma)) + 
        log(theta_new[2])
      c_k1 <- c_k1 +accept*u
      if(log(runif(1)) < accept){
        k <- 2
        # within model moves
        kappa_new <- sim_kappa(lambda, y, alphak, betak, theta_new[2], ak, bk)
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, al, bl)
      }else{
        # within model moves
        kappa_new <- NA
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, al, bl)
      }
    }else{
      theta <- c(lambda, kappa)
      theta_new <- lambda
      accept <- log(p) + pk1(y, theta_new, alphal, betal) - log(1-p) - 
        pk2(y, theta, alphal, betal, alphak, betak) + 
        log(dnorm(log(theta[2]/mu), 0, sigma)) - log(theta[2])
      c_k2 <- c_k2 + accept
      if(log(runif(1)) < accept){
        k <- 1
        # within model moves
        kappa_new <- kappa
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, al, bl)
      }else{
        # within model moves
        kappa_new <- sim_kappa(lambda, y, alphak, betak, kappa, ak, bk) 
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, al, bl)
      }
    }
    theta_save[i+1,] <- c(lambda_new, kappa_new, k)
    lambda <- lambda_new
    kappa <- kappa_new
    cat("\r", i)
  }
  print(c_k1)
  print(c_k2)
  return(theta_save)
}

##### Try it out #####
y <- football

# Parameters on priors
alphal <- 25
betal <- 10
alphak <- 1
betak <- 10

# Parameters for between-model step
mu = 0.015
sigma <- 1.5

# Proposal parameters for MH algorithms for lambda and kappa, k=2
al <- 30
bl <- 15
ak <- 2
bk <- 10

# Run MCMC
test <- rjmcmc_sampler(football, lambda0=2, kappa0=1, k0=1, mu=mu, 
                       sigma=sigma, alphal=alphal, betal=betal, alphak=alphak, 
                       betak=betak, al=al, bl=bl, ak=ak, bk=bk, mc.iter=50000)

test.s <- test[-c(1:5000),]

testk1 <- subset(test.s, k==1)

testk2 <- subset(test.s, k==2)

p1 <- round(nrow(testk1)/nrow(test.s),3)
p2 <- round(nrow(testk2)/nrow(test.s),3)

mk <- round(mean(testk2$kappa),3)
ml1 <- round(mean(testk1$lambda),3)
ml2 <- round(mean(testk2$lambda),3)

lambdak1 <- testk1$lambda
hist(lambdak1, freq=FALSE, col="lightblue", border="black", xlab=expression(lambda), main="")

lambdak2 <- testk2$lambda
hist(lambdak2, freq=FALSE, col="lightblue", border="black", xlab=expression(lambda),main="")

kappak2 <- testk2$kappa
hist(kappak2, freq=FALSE, col="lightblue", border="black", xlab=expression(kappa), main="")

p.hatk1 <- length(testk1[,1])/length(test[,1])
p.hatk2 <- length(testk2[,1])/length(test[,1])
bayes.factor.1 <- p.hatk1/p.hatk2
bayes.factor.2 <- p.hatk2/p.hatk1

bayes.factor.1
bayes.factor.2

barplot(c(p.hatk1, p.hatk2),names.arg=c("p(k=1 | data)", "p(k=2 | data)"), col= c("lightblue", "darkblue"), space=0.5, ylim=c(0,1), xlab="posterior probability of the models", ylab="probability") #cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

plot(lambdak1, type="l", col="lightblue", xlab="Iteration", ylab=expression(lambda))
lines(1:length(lambdak1), cumsum(lambdak1)/(1:length(lambdak1)), col="red")

plot(lambdak2, type="l", col="lightblue", xlab="Iteration", ylab=expression(lambda)) 
lines(1:length(lambdak2), cumsum(lambdak2)/(1:length(lambdak2)), col="red")

plot(kappak2, type="l", col="lightblue", xlab="Iteration", ylab=expression(kappa))
lines(1:length(kappak2), cumsum(kappak2)/(1:length(kappak2)), col="red")

p1 <- dpois(3, lambdak1)
p2 <- dnbinom(3, size=1/kappak2, mu=lambdak2)

hist(p1, freq=FALSE, col="gray", border="white", xlab="p, k=1, y=1", main="")
hist(p2, freq=FALSE, col="gray", border="white", xlab="p, k=2, y=1", main="")

