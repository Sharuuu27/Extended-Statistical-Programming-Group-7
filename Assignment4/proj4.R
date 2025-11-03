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
matrixes <- function(data, k = 80,...){
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
}