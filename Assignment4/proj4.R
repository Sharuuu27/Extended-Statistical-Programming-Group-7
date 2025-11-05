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

library(splines)
#setwd("file")
data <- read.table("engcov.txt",header = TRUE)

#Question1
# Define a function to compute X_tilde, X and S
matrixes <- function(data, k = 80){
  n <- nrow(data) # number of observation days, n
  t1 <- data$julian[1] - 30  # earliest infection day, setting as t1-30
  tn <- data$julian[n]       # latest infection day
  t_coverage <- t1:tn        # the scope of infection days
  mid_knots <- seq(from = t1, to = tn, length.out = k - 2) 
  # k-2 knots cover the scope of infection days
  unit <- (tn-t1)/(k-2-1) 
  knot <- seq(t1-unit*3, tn+unit*3,length.out=(k+4))
  # Adding 6 points beyond the original k+2 knots, ensuring a total of k+4 knots
  x_tilde <- splineDesign(knots=knot,x=t_coverage,ord=4,outer.ok=TRUE)
  # design matrix, X_tilde, for the B-splines defined by knots at the values 
  # in the scope of infection days
  s <- crossprod(diff(diag(k), diff=2)) 
  # smoothing penalty matrix S
  x <- matrix(0, nrow = n, ncol = k)
  # Initialise the death model matrix X
  d <- 1:80 # the range of days from infection to death.
  edur <- 3.151 # mean parameter for the log-normal distribution.
  sdur <- 0.469 # Standard deviation parameter for the log-normal distribution
  pd <- dlnorm(d, edur, sdur) # probability density
  pd <- pd / sum(pd) # Normalize probabilities pd
  for (i in 1:n) {
    j_max <- min(29 + i, 80) 
    #Determine the maximum interval j_max for the convolution limit.
    for (j in 1:j_max) {
      y <- 30+i-j 
      x[i,] <- x[i,]+x_tilde[y,]*pd[j]
      # X_{i,k}=\sum_{j=1}^{min(29+i,80)}\tilde{X}_{(30+i-j),k}\times \pi_{j}
    }
  }
  list(X_tilde = x_tilde, X = x, S = s)
}

data <- read.table("engcov.txt", header=TRUE)
data$date <- as.Date(data$date, format = "%d/%m/%y")

m <- matrixes(data)
x_tilde <- m$X_tilde
x <- m$X
s <- m$S
y <- data$deaths


#Question2

pnll <- function(gamma,x,y,s,lambda){
  beta <- exp(gamma)
  mu <- x %*% beta
  ll <- y*log(mu)-mu-lgamma(y+1)
  p <- lambda*(t(beta)%*%s%*%beta)/2
  -sum(ll) + p
}

gpnll <- function(gamma,x,y,s,lambda){
  beta <- exp(gamma)
  mu <- as.vector(x %*% beta)
  z <- y/mu-1
  dll <- as.vector(t(x) %*% z) * beta
  dp <- lambda * beta * as.vector(s %*% beta)
  -dll+dp
}

#check 
k <- 80

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
fit <- optim(par=gamma0,
             fn=pnll,
             gr=gpnll,
             lambda=5e-5,
             x=x,y=y,s=s,
             method="BFGS",
             hessian=TRUE)
gamma_h <- fit$par
beta_h <- exp(gamma_h)
mu_h <- x %*% beta_h
f_h<- x_tilde %*% beta_h


n <- nrow(data)
t1 <- data$julian[1] - 30
tn <- data$julian[n]
t_coverage <- t1:tn

plot(data$julian, y, type="p", pch=19, col="darkgray", cex=.8,
     xlab="Day", ylab="Deaths/Infection", 
     xlim=c(t1,tn), ylim=c(0,1750),
     main="")
lines(data$julian, mu_h, col="black", lwd=2)

lines(t_coverage, f_h, type="l", col="blue", lwd=2)
legend("topright", 
       legend=c("Observed Deaths", "Fitted Deaths", "Estimated Infection"), 
       col=c("darkgray", "black", "blue"), 
       pch=c(19, NA, NA), lty=c(NA, 1, 1), cex=0.8)


#Question4

