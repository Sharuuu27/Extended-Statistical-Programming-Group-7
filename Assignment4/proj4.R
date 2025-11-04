# Team Members: 
# 1. Shivasharini Govindasamy (Sharuuu27) - S2766392
# 2. Yu-Hsuan Hung (YuHsuan07) - S2793274
# 3. Yasuhiro Hara (hiroh-git) - S2826059

# Team Contributions:
# Shivasharini -  
# Yu-Hsuan -  
# Yasuhiro- 

# Github repo link:
#

## OVERVIEW: 


## CODE STRUCTURE:
# 1. 
# 2. 
# 3. 
# 4. 
# 5. 
# 6. 


#Question1
#setwd("file")
matrixes <- function(data, k = 80){
  n <- nrow(data)
  t1 <- data$julian[1] - 30
  tn <- data$julian[n]
  t_coverage <- t1:tn
  mid_knots <- seq(from = t1, to = tn, length.out = k - 2)
  knot <- c(rep(min(t_coverage),3), mid_knots, rep(max(t_coverage),3))
  x_tilde <- splineDesign(knots=knot,x=t_coverage,ord=4)
  s <- crossprod(diff(diag(k), diff=2))
  x <- matrix(0, nrow = n, ncol = k)
  d <- 1:80
  edur <- 3.151
  sdur <- 0.469
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  for (i in 1:n) {
    j_max <- min(29 + i, 80)
    for (j in 1:j_max) {
      y <- 30+i-j 
      x[i,] <- x[i,]+x_tilde[y,]*pd[j]
    }
  }
  list(X_tilde = x_tilde, X = x, S = s)
}

m <- matrixes(data)
x_tilde <- m$X_tilde
x <- m$X
s <- m$S
y <- data$deaths
#Question2

pnll <- function(gamma,lambda,x,y,s){
  beta <- exp(gamma)
  mu <- x %*% beta
  ll <- y*log(mu)-mu-lgamma(y+1)
  p <- lambda*(t(beta)%*%s%*%beta)/2
  -sum(ll) + p
}
gpnll <- function(gamma,lambda,x,y,s){
  beta <- exp(gamma)
  mu <- as.vector(x %*% beta)
  z <- y/mu-1
  dll <- as.vector(t(x) %*% z ) * beta
  dp <- beta * as.vector(s %*% beta)
  -dll+dp
}

#check 
fd <- gamma0 <- rep(log(1e-3),k)
lambda0 <- 5e-5
pnll0 <- pnll(gamma0,lambda0,x=x,y=y,s=s)
eps <- 1e-7
for (i in seq_along(gamma0)) {
  gamma1 <- gamma0
  gamma1[i] <- gamma1[i] + eps          
  pnll1 <- pnll(gamma1, lambda0, x = x, y = y, s = s)
  fd[i] <- (pnll1 - pnll0) / eps         
}
head(fd);head(gpnll(gamma0,lambda0,x=x,y=y,s=s))


#Question3
fit <- optim(gamma0,pnll,gr=gpnll,lambda=5e-5,x=x,y=y,s=s,method = "BFGS",hessian=TRUE)
gamma_h <- fit$par
beta_h <- exp(gamma_h)
mu_h <- x %*% beta_h
f_h<- x_tilde %*% beta_h


plot(t_coverage, f_h, type = "l", col = 2,
     xlab = "time", ylab = "deaths",
     ylim = range(c(f_h, data$deaths)))
lines(data$julian, data$deaths, col = "black")