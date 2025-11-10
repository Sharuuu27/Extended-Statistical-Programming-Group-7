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
data <- read.table("engcov.txt", header=TRUE)
data$date <- as.Date(data$date, format="%d/%m/%y")


## Evaluate X_tilde, X and S

matrixes <- function(data, k=80){
# define a function to compute X_tilde, X and S 
  n <- nrow(data) # number of observation days, n
  t1 <- data$julian[1] - 30  # earliest infection day, setting as t1-30
  tn <- data$julian[n]       # latest infection day
  t_coverage <- t1:tn        # the scope of infection days
  # k-2 knots cover the scope of infection days
  unit <- (tn-t1)/(k-2-1) 
  knot <- seq(t1-unit*3, tn+unit*3, length.out=(k+4))
  # adding 6 points beyond the original k+2 knots, ensuring a total of k+4 knots
  x_tilde <- splineDesign(knots=knot, x=t_coverage, ord=4, outer.ok=TRUE)
  # design matrix, X_tilde, for the B-splines defined by knots at the values 
  # in the scope of infection days
  s <- crossprod(diff(diag(k), diff=2)) 
  # smoothing penalty matrix S
  x <- matrix(0, nrow=n, ncol=k)
  # initialise the death model matrix X
  d <- 1:80 # the range of days from infection to death.
  edur <- 3.151 # mean parameter for the log-normal distribution.
  sdur <- 0.469 # standard deviation parameter for the log-normal distribution
  pd <- dlnorm(d, edur, sdur) # probability density
  pd <- pd / sum(pd) # normalize probabilities pd
  for (i in 1:n) { #loop over each observation days
    j_max <- min(29+i, 80) 
    # determine the maximum interval j_max for the convolution limit.
    for (j in 1:j_max) { # loop over all possible interval j
      y <- 30+i-j 
      x[i,] <- x[i,] + x_tilde[y,]*pd[j]
      # X_{i,k}=\sum_{j=1}^{min(29+i,80)}\tilde{X}_{(30+i-j),k}\times \pi_{j}
    }
  }
  list(X_tilde=x_tilde, X=x, S=s) # return the three computed matrices
}

# extract x_tilde, x, and s 
k <- 80
m <- matrixes(data, k=k)
x_tilde <- m$X_tilde
x <- m$X
s <- m$S
y <- data$deaths


## Preparation for optim()

pnll <- function(gamma,x,y,s,lambda){
# define the penalised negative log likelihood (pnll) function
  beta <- exp(gamma) # beta = exp(gamma), ensuring positive f
  mu <- x %*% beta   # mu = X %*% beta
  ll <- y*log(mu) - mu - lgamma(y+1) # log-likelihood
  p <- lambda*(t(beta) %*% s %*% beta)/2 # penalty
  -sum(ll) + p # return pnll= -sum(ll) + p
}

gpnll <- function(gamma,x,y,s,lambda){
# Define the gradient of pnll
  beta <- exp(gamma) # beta = exp(gamma), ensuring positive f
  mu <- as.vector(x %*% beta) # mu = X %*% beta
  z <- y/mu-1 # z = y/mu - 1
  dll <- as.vector(t(x) %*% z) * beta # d(sum l_i)/d(gamma)
  dp <- lambda * beta * as.vector(s %*% beta)# d(p)/d(gamma)
  -dll+dp # return total gradient:-[d(sum l_i)/d(gamma)]+[d(p)/d(gamma)]
}

#check 
fd <- gamma0 <- rep(log(10),k) # set gamma0 for testing
lambda0 <- 0.1 # set test lambda value.
pnll0 <- pnll(gamma0, lambda0, x=x, y=y, s=s) # pnll at gamma0 and lambda0
eps <- 1e-7 # finite difference interval 
for (i in seq_along(gamma0)) { # loop over gamma0
  gamma1 <- gamma0
  gamma1[i] <- gamma1[i] + eps  # increase gamma1[i] by eps         
  pnll1 <- pnll(gamma1, lambda0, x=x, y=y, s=s) # compute resulting pnll
  fd[i] <- (pnll1 - pnll0) / eps   # approximate -dl/dgamma[i]      
}
head(fd);head(gpnll(gamma0, lambda0, x=x, y=y, s=s)) 
# the result indicated that gpnll is coded correctly


## Plot the actual and fitted deaths and fitted daily infection curve

fit <- optim(par=gamma0, fn=pnll, gr=gpnll,
             lambda=5e-5, x=x, y=y, s=s,
             method="BFGS") #, hessian=TRUE
# Use optim function with BFGS method to optimise
gamma_h <- fit$par       # Get the fitted optimal parameters gamma_h
beta_h <- exp(gamma_h)   # Calculate the fitted B-spline coefficients beta_h
mu_h <- x %*% beta_h     # Calculate the fitted expected deaths mu_h
f_h<- x_tilde %*% beta_h # Calculate the fitted new infection curve f_h

# Set parameters
n <- nrow(data)
t1 <- data$julian[1] - 30
tn <- data$julian[n]
t_coverage <- t1:tn

# Plot actual deaths data (gray dots)
plot(data$julian, y, 
     type="p", 
     pch=19, 
     col="darkgray", 
     cex=.8,
     xlab="Day", 
     ylab="Deaths/Infection", 
     xlim=c(t1,tn), 
     ylim=range(c(f_h, data$deaths)))
# Plot fitted deaths data (black line)
lines(data$julian, mu_h, col="black", lwd=1.5)
# Plot fitted infection data (blue line)
lines(t_coverage, f_h, type="l", col="blue", lwd=1.5)
# Plot legend
legend("topright", 
       legend=c("Actual Deaths", "Fitted Deaths", "Fitted Infection"), 
       col=c("darkgray", "black", "blue"),
       pch=c(19, NA, NA), lty=c(NA, 1, 1), lwd=2, cex=.8)


## Choosing an approriate lambda using BIC

lsp <- seq(-13, -7, length=50) # Defining the grid
lambdas <- exp(lsp)
#n <- nrow(data) # number of observation
#k <- ncol(x) # number of parameters: 80

BIC <- EDF <- rep(0, length(lsp))
gamma_start <- gamma_h # Initial parameter（from Step 3）

for (i in 1:length(lambdas)) {
  # Grid search
  lambda_i <- lambdas[i]
 
  ## 1. Fit model
  fit <- optim(par=gamma_start, fn=pnll, gr=gpnll,
               lambda=lambda_i, x=x, y=y, s=s,
               method="BFGS") #, hessian=TRUE
  
  gamma_h <- fit$par
  beta_h <- exp(gamma_h)
  mu_h <- as.vector(x %*% beta_h) # fitted expected deaths mu_h (μ̂_i)
  
  ## 2. Calculate EDF
  W <- diag(y/(mu_h^2))
  
  XWX <- t(x) %*% W %*% x # H_0 = X^T W X + 0*S = X^T W X
  
  H_lambda <- XWX + lambda_i*s # H_lambda = X^T W X + lambda S
  
  H_inv_H0 <- solve(H_lambda, XWX) # (H_lambda)^(-1) %*% H_0
  
  EDF[i] <- sum(diag(H_inv_H0)) # trace(A) = sum(diag(A))
  
  ## 3. Calculate l(beta_h)
  P <- lambda_i*(t(beta_h) %*% s %*% beta_h)/2 
  # P = lambda * beta_h^T S * beta_h / 2
  
  ll <- -fit$value + P # ll = -nll = -pnll + P
  
  ## 4. BIC = -2l(beta_h) + log(n)*EDF
  BIC[i] <- -2*ll + log(n)*EDF[i]
}

## 5. Find optimal λ
i_opt <- which.min(BIC)      # index where BIC is minimum
lambda_opt <- lambdas[i_opt] # optimal lambda!
lsp_opt <- lsp[i_opt]        # optimal log(lambda)
EDF_opt <- EDF[i_opt]        # EDF at minimum BIC
BIC_min <- BIC[i_opt]        # minimum BIC


plot(lsp, BIC, type="l",  # Visualise the results
     xlab=expression(log(lambda)), 
     ylab="BIC", 
     main="BIC Optimization Results (Grid Search)")
points(lsp_opt, BIC_min, col="red", pch=19)

# Print the results
cat("--- Optimal Lambda Selection ---\n")
cat("Optimal log(lambda):", lsp_opt, "\n")
cat("Optimal lambda:", lambda_opt, "\n")
cat("Minimum BIC:", BIC_min, "\n")
cat("EDF at minimum BIC:", EDF_opt)

# Fit model with optimal lambda
fit_opt <- optim(par=gamma_start, fn=pnll, gr=gpnll,
                 lambda=lambda_opt, x=x, y=y, s=s,
                 method="BFGS") #, hessian=TRUE
# Use optim function with BFGS method to optimise
gamma_h_opt <- fit_opt$par       # Get the fitted optimal parameters gamma_h
beta_h_opt <- exp(gamma_h_opt)   # Calculate the fitted B-spline coefficients beta_h
mu_h_opt <- x %*% beta_h_opt     # Calculate the fitted expected deaths mu_h
f_h_opt <- x_tilde %*% beta_h_opt # Calculate the fitted new infection curve f_h

# Plot with optimal parameters
plot(data$julian, y, # Plot actual deaths data (gray dots)
     type="p", 
     pch=19, 
     col="darkgray", 
     cex=.8,
     xlab="Day", 
     ylab="Deaths/Infection", 
     xlim=c(t1,tn), 
     ylim=range(c(f_h_opt, data$deaths)))

lines(data$julian, mu_h_opt,  # Plot fitted deaths data (black line)
      col="black", 
      lwd=1.5)

lines(t_coverage, f_h_opt, # Plot fitted infection data (blue line)
      type="l", 
      col="blue", 
      lwd=1.5)

legend("topright", # Plot legend
       legend=c("Actual Deaths", "Fitted Deaths", "Fitted Infection"), 
       col=c("darkgray", "black", "blue"), 
       pch=c(19, NA, NA), lty=c(NA, 1, 1), lwd=2, cex=.8)


# Non-parametric bootstrapping for uncertainty

pnll_weight <- function(gamma, x, y, s, lambda, w) {
  # define the penalised negative log likelihood (pnll) function with weights
  beta <- exp(gamma) # ensuring f(t) is positive
  mu <- x %*% beta   
  ll <- y*log(mu) - mu - lgamma(y+1) 
  p <- lambda*(t(beta) %*% s %*% beta)/2 # penalty 
  -sum(w * ll) + p # return weighted pnll
}

gpnll_weight <- function(gamma, x, y, s, lambda, w) {
  # Define the gradient of weighted pnll
  beta <- exp(gamma) 
  mu <- as.vector(x %*% beta)
  z <- y/mu - 1 
  dll <- as.vector(t(x) %*% (w * z)) * beta #(z is multiplied with weights)
  dp <- lambda * beta * as.vector(s %*% beta) 
  -dll + dp # return total weighted gradient
}

n_bootstrap <- 200 # to generate 200 bootstrap replicates

fhat_bootstrap <- matrix(0, nrow = nrow(x_tilde), ncol = n_bootstrap) 
# creating matrix for the bootstrap

for (i in 1:n_bootstrap) {
  # Generate bootstrap weights for replicate = i
  wb <- tabulate(sample(n, replace=TRUE), n)
  
  # Fit the model using the weighted functions
  fit_b_hat <- optim(par=gamma_h_opt,  # Optimal parameter
                 fn=pnll_weight,       # New weighted pnll
                 gr=gpnll_weight,      # New weighted gradient
                 lambda=lambda_opt,    # Fixed optimal lambda from Q4
                 x=x, y=y, s=s, w=wb,
                 method="BFGS")
  
  gamma_b <- fit_b_hat$par # parameters from this bootstrap fit
  beta_b <- exp(gamma_b)
  
  f_hat_b <- x_tilde %*% beta_b # fitted infection curve for this replicate
  
  fhat_bootstrap[, i] <- f_hat_b # store the result of infection curve
}


# Final results plot with confidence intervals

fitted_lower <- apply(fhat_bootstrap, 1, quantile, probs = 0.025) 
# lower confidence bound for infection in each day

fitted_upper <- apply(fhat_bootstrap, 1, quantile, probs = 0.975)
# upper confidence bound for infection in each day

ylim_max <- max(c(y, fitted_upper), na.rm = TRUE) 
# ensure the canvas is big enough to fit the entire plot

plot(data$julian, y, 
     type = "p", 
     pch = 19, 
     col = "darkgray", 
     cex = 0.8,
     xlab = "Day", 
     ylab = "Deaths / Infections",
     main = "Fitted COVID-19 Deaths and Infections with 95% CI",
     xlim = c(t1, tn), 
     ylim = c(0, ylim_max))


polygon(c(t_coverage, rev(t_coverage)), 
        # outlines the confidence band width along the t-coverage
        c(fitted_upper, rev(fitted_lower)),
        col = rgb(0.2, 0.2, 1, 0.2),
        border = NA)

lines(data$julian, mu_h_opt, # Plot fitted deaths data (black line)
      col = "black", 
      lwd = 1.5)

lines(t_coverage, f_h_opt, # Plot fitted infection data (blue line)
      col = "blue", 
      lwd = 1.5)

legend("topright", # plot legend
       legend = c("Actual Deaths", "Fitted Deaths", "Fitted Infections",
                  "Infections (95% CI))"),
       col = c("darkgray", "black", "blue", rgb(0.2, 0.2, 1, 0.2)),
       pch = c(19, NA, NA, NA), 
       lwd = c(NA, 2, 2, 6),
       fill = c(NA, NA, NA, NA), 
       border = NA, 
       cex = 0.8) 