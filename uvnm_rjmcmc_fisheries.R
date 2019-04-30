##############################
# uvnm.rjmcmc with fisheries dataset
##############################

rdirichlet <- function (n, alpha) {
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rtriang")
  
  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

rMultinom <- function(p){
  if(any(p < 0))
    stop("non-positive probability")
  if(!isTRUE(all.equal(rowSums(p),rep(1, nrow(p)))))
    stop("The sums of probabilities have to be 1.")
  k <- ncol(p)
  n <- nrow(p)
  csm <- ifelse(outer(1:k,1:(k - 1),FUN=`<=`), 1, 0)
  cp <- p %*% csm
  U <- runif(n)
  rowSums(sweep(cp, 1, U, `<=`)) + 1
}



######################################
## decide whether to split-combine and birth-death
######################################

decideScBd <- function(k, kmax){
  if(k==1)
    action <- 1
  else{
    if(k==kmax)
      action <- 0
    else
      action <- ifelse(runif(1) <0.5, 1, 0)
  }
}

################################################
## get the probabilities for z in the split and combine cases
################################################

getP <- function(yjstar, muj1, muj2, sigma2j1, sigma2j2, wj1, wj2){
  p <- cbind(dnorm(yjstar, muj1, sqrt(sigma2j1)),
             dnorm(yjstar, muj2, sqrt(sigma2j2)))
  p[,1] <- p[,1] * wj1
  p[,2] <- p[,2] * wj2                 
  s <- rowSums(p)
  p[which(s == 0 | !is.finite(s)),] <- 1
  p / rowSums(p)
}

################################################
## get acceptance probabilities for split-combine
################################################

getAsc <- function(yjstar, j1.pos, j2.pos, k,
                   wj1, wj2, wjstar,
                   muj1, muj2, mujstar,
                   sigma2j1, sigma2j2, sigma2jstar,
                   delta, xi, kappa, alpha, beta,
                   p, b, d, u1, u2, u3){
  
  ## calculate acceptance probability
  ## take log first then transform to exponential. Cause there are
  ## multiplication of larger number of densities and so on.
  
  #likelihood ratio
  log.lr <- sum(dnorm(yjstar[j1.pos], mean=muj1, sd=sqrt(sigma2j1), log=TRUE)) +
    sum(dnorm(yjstar[j2.pos], mean=muj2, sd=sqrt(sigma2j2), log=TRUE)) -
    sum(dnorm(yjstar, mean=mujstar, sd=sqrt(sigma2jstar), log=TRUE))
  l1 <- length(j1.pos)
  l2 <- length(j2.pos)
  
  #term 1-3: ratio between two states
  #term1 includes various ratios
  log.term1 <- log.lr +    #likelihood ratio
    log(1) +               #ratio of p(k+1)/p(k)
    (delta-1+l1)*log(wj1) + (delta-1+l2)*log(wj2) - (delta-1+l1+l2)*log(wjstar) - #ratio of z
    lbeta(delta, k*delta) #ratio of w
  
  #term2 is the ratio of mu
  log.term2 <- log(k+1) +
    0.5*log(kappa/(2*pi)) -
    0.5*kappa*((muj1-xi)*2 + (muj2-xi)*2 - (mujstar-xi)*2) 
  
  #term3 is the ratio of sigma
  log.term3 <- alpha*log(beta) - lgamma(alpha) +
    (-alpha-1)*(log(sigma2j1) + log(sigma2j2) - log(sigma2jstar)) -
    beta*(1/sigma2j1 + 1/sigma2j2 - 1/sigma2jstar)
  
  #term4 is the transform ratio between two states
  log.Palloc <- sum(log(p[,1][j1.pos])) + sum(log(p[,2][j2.pos]))
  log.term4 <- log(d[k+1]) - log(b[k]) - log.Palloc -
    dbeta(u1, 2, 2, log=TRUE) - dbeta(u2, 2, 2, log=TRUE) - dbeta(u3, 1, 1, log=TRUE)
  
  #term5 is the Jacobian
  log.term5 <- log(wjstar) +log(abs(muj1-muj2)) + log(sigma2j1) + log(sigma2j2) -
    log(u2) -log(1-u2^2) - log(u3) - log(1-u3) - log(sigma2jstar)
  
  exp(log.term1 + log.term2 + log.term3 + log.term4 + log.term5)
}

######################################
## SPLIT MOVE
######################################

split <- function(y, k, w, mu, sigma2, Z, b, d, delta, xi, kappa, alpha, beta){
  #choose a component to split
  jstar <- sample(1:k, 1)
  
  s.acc <-0
  
  #generate intermidiate parameters
  u1 <- rbeta(1, 2, 2)
  u2 <- rbeta(1, 2, 2)
  u3 <- rbeta(1, 1, 1)
  
  #generate two new ws
  wjstar <- w[jstar]
  wj1 <- wjstar * u1
  wj2 <- wjstar * (1-u1)
  
  
  #generate two new mus
  mujstar <- mu[jstar]
  sigma2jstar <- sigma2[jstar]
  muj1 <- mujstar - u2*sqrt(sigma2jstar)*sqrt(wj2/wj1)
  muj2 <- mujstar + u2*sqrt(sigma2jstar)*sqrt(wj1/wj2)
  
  #check order of mu
  #newmu.part <- c(mu[jstar-1], muj1, muj2, mu[jstar+1])
  #if(all(order(newmu.part) == 1:length(newmu.part))){
  #generate two new sigma2s
  sigma2j1 <- u3 * (1-u2^2) * sigma2jstar * wjstar / wj1
  sigma2j2 <- (1-u3) * (1-u2^2) * sigma2jstar * wjstar / wj2
  
  #allocate z_i=jstar
  jstar.pos <- which(Z==jstar)
  if(length(jstar.pos) > 0){
    yjstar <- y[jstar.pos]
    p <- getP(yjstar, muj1, muj2, sigma2j1, sigma2j2, wj1, wj2)
    zj12 <- rMultinom(p)
    j1.pos <- which(zj12==1)
    j2.pos <- which(zj12==2)
    
    A <- getAsc(yjstar, j1.pos, j2.pos, k,
                wj1, wj2, wjstar,
                muj1, muj2, mujstar,
                sigma2j1, sigma2j2, sigma2jstar,
                delta, xi, kappa, alpha, beta,
                p, b, d, u1, u2, u3)
    
    if(runif(1) < min(1, A)){
      s.acc <- 1
      indicator <- rep(0,k)
      indicator[jstar] <- 1
      indicator[which((1:k) > jstar)] <- 2
      ind0 <- which(indicator==0)
      ind2 <- which(indicator==2)
      
      #generate new w vector
      w <- c(w[ind0], wj1, wj2, w[ind2])
      #generate new mu vector
      
      #TRACK CHANGE 'sample'
      mu <- sample(c(mu[ind0], muj1, muj2, mu[ind2]))
      #generate new sigma2 vector
      sigma2 <- c(sigma2[ind0], sigma2j1, sigma2j2, sigma2[ind2])
      #gernerate new Z matrix
      larger <- which(Z > jstar)
      Z[larger] <- Z[larger] + 1
      Z[jstar.pos[j1.pos]] <- jstar
      Z[jstar.pos[j2.pos]] <- jstar + 1
      
      
    }
  }
  
  list(w=w, mu=mu, sigma2=sigma2, Z=Z, s.acc=s.acc)
}

######################################
## COMBINE MOVE
######################################
combine <- function(y, k, w, mu, sigma2, Z, b, d, delta, xi, kappa, alpha, beta){
  #choose a pair of components to combine
  j1 <- sample(1:(k-1), 1)
  j2 <- j1 + 1
  
  # combine accepted
  c.acc <- 0
  
  #generate new parameters
  wj1 <- w[j1]
  wj2 <- w[j2]
  muj1 <- mu[j1]
  muj2 <- mu[j2]
  sigma2j1 <- sigma2[j1]
  sigma2j2 <- sigma2[j2]
  wjstar <- w[j1] + w[j2]
  mujstar <- (wj1*muj1 + wj2*muj2) / wjstar
  sigma2jstar = wj1*wj2*((muj1-muj2)/wjstar)**2 +
    (wj1*sigma2j1+wj2*sigma2j2)/wjstar
  
  #calculate acceptance probability
  #likelihood ratio
  jstar.pos <- which(Z==j1 | Z==j2)
  j1.pos <- match(which(Z==j1), jstar.pos)
  j2.pos <- match(which(Z==j2), jstar.pos)
  yjstar <- y[jstar.pos]
  p <- getP(yjstar, muj1, muj2, sigma2j1, sigma2j2, wj1, wj2)
  
  #generate intermidiate parameters
  #u1 u2 u3 are derived from split move (equations below eqn (10) on page 739
  u1 <- rbeta(1, 2, 2)
  u2 <- rbeta(1, 2, 2)
  u3 <- rbeta(1, 1, 1)
  #u1 <- wj1 / wjstar
  #u2 = (muj2-muj1) * sqrt(wj1*wj2/sigma2jstar) / wjstar
  #u2 = max(u2,1e-12)
  #u2 = min(u2,1.0-1e-4)
  #u3 = wj1*sigma2j1 / (wj1*sigma2j1+wj2*sigma2j2)
  
  A <- getAsc(yjstar, j1.pos, j2.pos, k-1,
              wj1, wj2, wjstar,
              muj1, muj2, mujstar,
              sigma2j1, sigma2j2, sigma2jstar,
              delta, xi, kappa, alpha, beta,
              p, b, d, u1, u2, u3)
  
  # ORDERING
  if(runif(1) < min(1, 1/A)){
    c.acc <- 1
    indicator <- rep(0,k)
    indicator[c(j1,j2)] <- 1
    indicator[which((1:k) > j2)] <- 2
    ind0 <- which(indicator==0)
    ind2 <- which(indicator==2)
    #generate new w vector
    w <- c(w[ind0], wjstar, w[ind2])
    #generate new mu vector
    
    # TRACK CHANGE 'sample'
    mu <- sample(c(mu[ind0], mujstar, mu[ind2]))
    #generate new sigma2 vector
    sigma2 <- c(sigma2[ind0], sigma2jstar, sigma2[ind2])
    #gernerate new Z matrix
    Z[which(Z==j2)] <- j1
    large <- which(Z > j2)
    Z[large] <- Z[large] - 1
  }
  list(w=w, mu=mu, sigma2=sigma2, Z=Z, c.acc=c.acc)
}


######################################
## Acceptance probability for birth-death move
######################################

getAbd <- function(n, k, k0, delta, wjstar, b, d){
  log.term1 <- log(1) - lbeta(k*delta, delta) + (delta-1)*log(wjstar) +
    (n+k*delta-k)*log(1-wjstar) + log(k+1)
  #note that there is an error in the original paper:
  # (1-wjstar)^(k-1) instead of (1-wjstar)^k
  log.term2 <- log(d[k+1]) - log(k0+1) - log(b[k]) - dbeta(wjstar,1,k, log=TRUE) + (k-1)*log(1-wjstar)
  exp(log.term1 + log.term2)
}


######################################
## BIRTH MOVE
######################################
birth <- function(n, k, w, mu, sigma2, Z, delta, xi, kappa, alpha, beta, b, d){
  
  b.acc <-0
  wjstar <- rbeta(1, 1, k)
  k0 <- sum(unlist(lapply(1:k, function(i) sum(Z==i)))==0)
  A <- getAbd(n, k, k0, delta, wjstar, b, d)
  if(runif(1) < min(1, A)){
    
    b.acc <- 1
    
    mujstar <- rnorm(1, mean=xi, sd=sqrt(1/kappa))
    sigma2jstar <- rgamma(1, shape=alpha, rate=beta)
    sigma2jstar <- 1 / sigma2jstar
    w <- w*(1-wjstar)
    jstar.pos <- which(mu > mujstar)[1]
    if(is.na(jstar.pos)){
      w <- c(w, wjstar)
      mu <- c(mu, mujstar)
      sigma2 <- c(sigma2, sigma2jstar)
    }
    else{
      indicator <- rep(0,k)
      indicator[which((1:k) >= jstar.pos)] <- 1
      ind0 <- which(indicator==0)
      ind1 <- which(indicator==1)
      w <- c(w[ind0], wjstar, w[ind1])
      
      #TRACK CHANGE 'sample'
      mu <- sample(c(mu[ind0], mujstar, mu[ind1]))
      sigma2 <- c(sigma2[ind0], sigma2jstar, sigma2[ind1])
      larger <- which(Z >= jstar.pos)
      Z[larger] <- Z[larger] + 1
    }
  }
  list(w=w, mu=mu, sigma2=sigma2, Z=Z, b.acc=b.acc)
}


######################################
## DEATH MOVE
######################################
death <- function(n, k, w, mu, sigma2, Z, delta, xi, kappa, alpha, beta, b, d){
  
  d.acc <-0
  d.candidate <- which(unlist(lapply(1:k, function(i) sum(Z==i)==0)))
  if(length(d.candidate) > 0){
    d.pos <- sample(1:length(d.candidate),1)
    d.pos <- d.candidate[d.pos]
    wjstar <- w[d.pos]
    k0 <- length(d.candidate) - 1
    A <- getAbd(n, k-1, k0, delta, wjstar, b, d)
    if(runif(1) < min(1, 1/A)){
      
      d.acc <- 1
      w <- w[-d.pos]
      w <- w / sum(w)
      mu <- mu[-d.pos]
      sigma2 <- sigma2[-d.pos]
      larger <- which(Z > d.pos)
      Z[larger] <- Z[larger] - 1
    }
  }
  list(w=w, mu=mu, sigma2=sigma2, Z=Z, d.acc=d.acc)
}

############################
# RJ MCMC sampler
###########################
uvnm.rjmcmc <- function(y, nsweep, kmax, k, w, mu, sigma2, Z,
                        delta=1, xi=NULL, kappa=NULL, alpha=2,
                        beta=NULL, g=0.2, h=NULL, verbose=TRUE){
  #Error checking
  if(nsweep <= 0)
    stop("The number of sweeps has to be positive.")
  if(kmax < 0 || k < 0)
    stop("The number of components have to be positive.")
  if(kmax < k)
    stop("The maximum number of components allowed is larger than
         the intitial value of components.")
  if(length(w) != k){
    w <- rep(w, len=k)
    warning("The length of 'w' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(any(w <=0))
    stop()
  if(sum(w) != 1){
    
  }
  if(length(mu) != k){
    mu <- rep(mu, len=k)
    warning("The length of 'mu' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(length(sigma2) != k){
    sigma2 <- rep(sigma2, len=k)
    warning("The length of 'sigma2' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(any(sigma2 <= 0))
    stop
  if(length(Z) != length(y)){
    Z <- rep(Z, len=length(y))
    warning("The length of 'Z' was not equal to the length of 'y' and
            was forced to be equal by being cut off or recycled.")
  }
  
  R <- diff(range(y))
  if(is.null(xi))
    xi <- median(y)
  if(is.null(kappa))
    kappa <- 1/R^2
  if(is.null(alpha))
    alpha <- 2
  if(is.null(g))
    g <- 0.2
  if(is.null(h))
    h <- 20/R^2
  if(is.null(beta))
    beta <- rgamma(1, shape=g, rate=h)
  
  n=length(y)
  
  cat(paste("R=",R,"\n", sep=""))
  cat(paste("xi=", xi,"\n", sep=""))
  cat(paste("kappa=",kappa,"\n", sep=""))
  cat(paste("alpha=",alpha, "\n", sep=""))
  cat(paste("g=",g, "\n", sep=""))
  cat(paste("h=",h, "\n",sep=""))
  
  split.move <- 0
  combine.move <- 0
  birth.move <- 0
  death.move <- 0
  
  s.acc <- 0
  c.acc <- 0
  b.acc <- 0
  d.acc <- 0
  
  
  #split probabilities
  b <- rep(0.5, kmax)
  b[kmax] <- 0
  #combine probabilities
  d <- rep(0.5, kmax)
  d[1] <- 0
  
  k.save <- rep(0, nsweep)
  w.save <- mu.save <- sigma2.save <- vector("list", nsweep)
  Z.save <- matrix(0, nrow=n, ncol=nsweep)
  for(i in 1:nsweep){
    k <- length(mu)
    
    #update w|...
    Z.expand <- do.call(cbind, lapply(1:k, function(i) ifelse(Z==i, 1, 0)))
    Nj <- colSums(Z.expand)
    w <- rdirichlet(1, delta + Nj)
    
    #update mu|...
    sumYj <- colSums(y * Z.expand)
    precision <- Nj/sigma2 + kappa
    mean <- (sumYj/sigma2 + kappa*xi) / precision
    mu.new <- rnorm(k, mean=mean, sd=sqrt(1/precision))
    if(all(order(mu.new) == 1:k)){
      mu <- mu.new
    }
    
    #update simga2|...
    Diff2j <- outer(y, mu, `-`) ^ 2
    sumDiff2j <- colSums(Diff2j * Z.expand)
    sigma2 <- rgamma(k, shape=alpha + Nj/2, rate=beta + sumDiff2j/2)
    sigma2 <- 1 / sigma2
    
    #update Z
    #p <- exp(- Diff2j %*% diag(1/(2*sigma2), nrow=k, ncol=k))
    #p <- p %*% diag(w/sqrt(sigma2), nrow=k, ncol=k)
    p <- do.call(cbind, lapply(1:k, function(i)
      w[i] * dnorm(y, mu[i], sqrt(sigma2[i]))))
    s <- rowSums(p)
    p[which(s == 0 | !is.finite(s)),] <- 1
    p <- p / rowSums(p)
    if(ncol(p) > 1){
      Z <- rMultinom(p)
    }
    else{
      Z <- rep(1, n)
    }
    
    #update beta
    beta <- rgamma(1, shape=g + k*alpha, rate=h + sum(1/sigma2))
    
    #combine or split
    action <- decideScBd(k, kmax)
    if(action==1){
      split.results <- split(y, k, w, mu, sigma2, Z, b, d, delta, xi, kappa, alpha, beta)
      split.move <- split.move +1
      s.acc <- s.acc + split.results$s.acc
      w <- split.results$w
      mu <- split.results$mu
      sigma2 <- split.results$sigma2
      Z <- split.results$Z
      k <- length(w)
    }
    else{
      combine.results <- combine(y, k, w, mu, sigma2, Z, b, d, delta, xi, kappa, alpha, beta)
      combine.move <- combine.move + 1
      c.acc <- c.acc + combine.results$c.acc
      w <- combine.results$w
      mu <- combine.results$mu
      sigma2 <- combine.results$sigma2
      Z <- combine.results$Z
      k <- length(w)
    }
    
    #birth-death
    action <- decideScBd(k, kmax)
    if(action==1){
      birth.results <- birth(n, k, w, mu, sigma2, Z, delta,
                             xi, kappa, alpha, beta, b, d)
      birth.move <- birth.move +1
      b.acc <- b.acc + birth.results$b.acc
      w <- birth.results$w
      mu <- birth.results$mu
      sigma2 <- birth.results$sigma2
      Z <- birth.results$Z
      k <- length(w)
    }
    else{
      
      death.results <- death(n, k, w, mu, sigma2, Z, delta,
                             xi, kappa, alpha, beta, b, d)
      death.move <- death.move + 1
      d.acc <- d.acc + death.results$d.acc
      w <- death.results$w
      mu <- death.results$mu
      sigma2 <- death.results$sigma2
      Z <- death.results$Z
      k <- length(w)
    }
    
    k.save[i] <- k
    w.save[[i]] <- w
    mu.save[[i]] <- mu
    sigma2.save[[i]]<- sigma2
    Z.save[,i] <- Z
    
    if (verbose && i %% 10000 == 0){
      cat(paste(i, " sweeps", " have finished.\n", sep=""))
    }
    
    
  }
  
  list(k.save=k.save, w.save=w.save, mu.save=mu.save,
       sigma2.save=sigma2.save, Z.save=Z.save, s.move=split.move, 
       b.move=birth.move, c.move=combine.move, d.move=death.move,
       s.acc=s.acc, c.acc=c.acc,b.acc=b.acc,d.acc=d.acc)
}


#################################
# FISHERIES DATA
#################################

library(mixAK)
library(coda)

N <- 200
components <- sample(1:3, prob=c(0.3,0.4,0.3), size=N, replace=TRUE)
mu <- c(2,5,10)
sd <- sqrt(c(0.5,0.1,2))
samples <- rnorm(n=N, mean=mu[components], sd=sd[components])

hist(samples, breaks=50, main="", xlab="Simulated data", col="lightblue", freq = FALSE)
summary(samples)
fish
data("fish", package="bayesmix")
f <- fish[,1]
f
samples
################
# save
################

hist(f, breaks=30, col="orange")

################################
# INITALISE RJMCMC SAMPLER
################################

w <- c(1/3,1/3,1/3)
mu <- c(3,5,9)
sigma2 <- c(1,1,1)

Z <- do.call(cbind, lapply(1:3, function(i)
  w[i]*dnorm(f, mu[i], sqrt(sigma2[i]))))
Z <- apply(Z, 1, function(x) which(x==max(x))[1])

Z




result <- uvnm.rjmcmc(y=f, nsweep=100000, kmax=30, k=3, w=w,
                      mu=mu, sigma2=sigma2, Z=Z,
                      delta=1, xi=NULL, kappa=NULL, alpha=2,
                      beta=NULL, g=0.2, h=NULL, verbose=TRUE)


####################################
# CALCULATE MOVES
####################################

result$b.move
result$d.move
result$s.move
result$c.move

result$b.acc
result$d.acc
result$s.acc
result$c.acc


acc.p <- c(result$s.move/result$s.acc, result$c.move/result$c.acc, result$b.move/result$b.acc, result$d.move/result$d.acc)
acc.p

###############################
# Raw plot and bar plot for K
###############################

k <- c()
for (i in 1:12){
  k <- cbind(k,length(which(result$k.save==i)))
}
k/100000

# 2 BARPLOT OF K
barplot(k, width = 1, names.arg=c("1","2","3","4","5","6","7","8","9","10","11","12"), xlab=("Number of components"), ylab=("Frequency"), col="orange")

# TRACE PLOT OF K
plot(result$k.save, xlab= "Iterations", ylab= "k", type="l", col="orange") # Does this model favour a larger number of compoents

######################################
# PLOTS FOR EACH COMPONENT 
######################################
w1 <- c()
w2 <- c()
w3 <- c()
w4 <- c()

mu1 <- c()
mu2 <- c()
mu3 <- c()
mu4 <- c()

sdev1<-c()
sdev2<-c()
sdev3<-c()
sdev4<-c()


for (i in (1:100000)) {
  if (length(result$mu.save[[i]])==4) {
    w1 <- cbind(w1,result$w.save[[i]][1])
    w2 <- cbind(w2,result$w.save[[i]][2])
    w3 <- cbind(w3,result$w.save[[i]][3])
    w4 <- cbind(w4,result$w.save[[i]][4])
    
    mu1 <- cbind(mu1, result$mu.save[[i]][1])
    mu2 <- cbind(mu2, result$mu.save[[i]][2])
    mu3 <- cbind(mu3, result$mu.save[[i]][3])
    mu4 <- cbind(mu4, result$mu.save[[i]][4])
    
    sdev1 <- cbind(sdev1,sqrt(result$sigma2.save[[i]][1]))
    sdev2 <- cbind(sdev2,sqrt(result$sigma2.save[[i]][2]))
    sdev3 <- cbind(sdev3,sqrt(result$sigma2.save[[i]][3]))
    sdev4 <- cbind(sdev4,sqrt(result$sigma2.save[[i]][4]))
  }
  else {
    
  }
}


# Reversible Jump MCMC output on the parameters of the model for simulating data, condition on k=4. 

# TRACE PLOTS FOR EACH PARAMETER

plot(w1[1:length(w1)], type='l',ylim=c(0,1))
lines(w2[1:length(w2)], type='l')
lines(w3[1:length(w3)], type='l')
lines(w4[1:length(w4)], type='l')

plot(mu1[1:length(mu1)], type='l',ylim=c(1,13), ylab="Means", col="orange")
lines(mu2[1:length(mu2)], type='l',col="orange")
lines(mu3[1:length(mu3)], type='l',col="orange")
lines(mu4[1:length(mu4)], type='l',col="orange")

#plot(sdev1[1:length(sdev1)], type='l',ylim=c(-1,12))
#lines(sdev2[1:length(sdev2)], type='l')
#lines(sdev3[1:length(sdev3)], type='l')
#lines(sdev4[1:length(sdev4)], type='l')


# CUMULAIVE AVERAGE OF WEIGTHS CONDITIONED ON K=4 
plot(1, type="n", xlim=c(1,length(w1)), ylim=c(0,1), ylab="Weights")
lines(1:length(w1), cumsum(w1) / (1:length(w1)), col="red")
lines(1:length(w2), cumsum(w2) / (1:length(w2)), col="red")
lines(1:length(w3), cumsum(w3) / (1:length(w3)), col="red")
lines(1:length(w4), cumsum(w4) / (1:length(w4)), col="red")

# CUMULAIVE AVERAGE OF MEANS CONDITIONED ON K=4
plot(1, type="n", xlim=c(1,length(mu1)), ylim=c(0,10), ylab="Means")
lines(1:length(mu1), cumsum(mu1) / (1:length(mu1)), col="red")
lines(1:length(mu2), cumsum(mu2) / (1:length(mu2)), col="red")
lines(1:length(mu3), cumsum(mu3) / (1:length(mu3)), col="red")
lines(1:length(mu4), cumsum(mu4) / (1:length(mu4)), col="red")

# CUMULAIVE AVERAGE OF WEIGTHS CONDITIONED ON K=4 
plot(1, type="n", xlim=c(1,length(sdev1)), ylim=c(0,2), ylab="Standard deviation")
lines(1:length(sdev1), cumsum(sdev1) / (1:length(sdev1)), col="red")
lines(1:length(sdev2), cumsum(sdev2) / (1:length(sdev2)), col="red")

lines(1:length(sdev3), cumsum(sdev3) / (1:length(sdev3)), col="red")
lines(1:length(sdev4), cumsum(sdev4) / (1:length(sdev4)), col="red")


# HISTOGRAMS FOR EACH PARAMETER
# Density estimates conditioned on k=4 for mean, weights and standard dev

ksave <- result$k.save
focus.k <- 12
pick.k <- which(ksave==4)


mu_hist <- unlist(result$mu.save[pick.k])
sigma2_hist <- unlist(result$sigma2.save[pick.k])
w_hist <- unlist(result$w.save[pick.k])

hist(w_hist, breaks=50, main="", xlab="Weights", col="orange")
hist(mu_hist, breaks=100,xlim=c(-1,12),main="", xlab="Means", col="orange")
hist(sigma2_hist, breaks=50, main="", xlab="Standard deviation",col="orange")

#DESITY ESTIMATE

post <- dnorm(n=length(mu_hist), mean=mu_hist, sd=sqrt(sigma2_hist))
hist(f, breaks=30, main="", xlab="Fisheries data", col="orange", freq = FALSE)
#lines(density(post),na.rm=TRUE, lwd=2, lty=1)
lines(density(post),na.rm=TRUE, lwd=1, lty=7)

w.mat <- matrix(c(w1,w2,w3,w4), ncol=4, nrow=length(w1), byrow=FALSE)
mu.mat <- matrix(c(mu1,mu2,mu3,mu4), ncol=4, nrow=length(mu1), byrow=FALSE)
sdev.mat <- matrix(c(sdev1,sdev2,sdev3,sdev4), ncol=4, nrow=length(sdev1), byrow=FALSE)  


colMeans(w.mat)
colMeans(mu.mat)
colMeans(sdev.mat)^2

#conditional density estimation
focus.k <- 8
pick.k <- which(ksave==focus.k)
w <- unlist(result$w.save[pick.k])
mu <- unlist(result$mu.save[pick.k])
sigma2 <- unlist(result$sigma2.save[pick.k])
den.estimate <- rep(w, each=length(f)) *
  dnorm(rep(f, length(w)), mean=rep(mu, each=length(f)),
        sd=rep(sqrt(sigma2), each=length(f)))
den.estimate <- rowMeans(matrix(den.estimate, nrow=length(f)))*focus.k

#within-sample classification
class <- apply(result$Z.save[,pick.k], 1,
               function(x) c(sum(x==1), sum(x==2), sum(x==3)))
class <- max.col(t(class))
#par(mfrow=c(1,1))

#visualize the results
#hist(f, freq=FALSE, breaks=20, main="", ylim=c(0,0.35), col="orange")
lines(sort(f), den.estimate[order(f)], na.rm=TRUE, lwd=1, lty=3)



