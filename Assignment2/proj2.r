# Team Members: 
# 1. Shivasharini Govindasamy (Sharuuu27) - S2766392
# 2. Yu-Hsuan Hung (YuHsuan07) - S2793274
# 3. Yasuhiro Hara (hiroh-git) - S2826059

# Team Contributions:
# Yasuhiro - 
# Yu-Hsuan - 
# Shivasharini - 

## Code to [...].
## The challenge is that [...].
## 3 strategies are suggested:
## 1. [...].
## 2. [...].
## 3. [...].
## [...]. 

set.seed(0)
## Question 1
# To assign the n people to thier unique household ID

household_id <- function(n = 1000, hmax = 5) {
  # n = total population
  # hmax = maximum members in a household
  h <- rep(1:n, times = sample(1:hmax, n, replace = TRUE))
  # the size of households are uniformly distributed
  return(h[1:n])
  # to make sure the population = n
}

h<-household_id()

## Question 2
# To assign Sociability Parameter (Beta-i)
generate_beta <- function(n = 1000, bmu = 5e-5, bsc = 1e-5) {
  beta <- rgamma(n, shape=bmu/bsc, scale=bsc)
  return(beta)
}  

beta <- generate_beta()

# To assign regular network relations

get_net <- function(beta, h, nc = 15) {
  
  n <- length(beta)
  mean_beta <- mean(beta)
  prob_denominator <- mean_beta^2 * (n - 1)
  
  network <- replicate(n, integer(0), simplify = FALSE) 
  #creating network list of all n people
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) { # to make sure the pairs are not itself
      if (h[i] == h[j]) { 
        next # if i and j of the same household, go to the next pair
      }
      
      link <- (nc * beta[i] * beta[j]) / prob_denominator 
      # probability of contact network
      
      link <- min(link, 1) # to ensure the probability does not exceed 1
      
      if (runif(1) < link) { # to check weather i and j are regular contacts
        network[[i]] <- c(network[[i]], j) # assign i to j's contact list
        network[[j]] <- c(network[[j]], i) # assign j to i's contact list
      }
    }
  }
  
  return(network)
}

get_net(beta, h)


