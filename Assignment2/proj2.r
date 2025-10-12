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

get.net <- function(beta, h, nc = 15) {
  
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

alink <- get.net(beta, h)


## Question3
## SEIR simulation model with social structure
nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, 
                  gamma = .4, nc = 15,nt = 100, pinf = .005) {
  ## alpha = the daily probs I[j] -> S[i]
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## nc = Average contacts; nt = number of days
  ## pinf = initial proportion infected.
  
  ##Initialization 
  x <- rep(0, n)
  ni <- round(n * pinf)
  
  ## Randomly set the initial number of infected individuals (I=2)
  x[sample(n, ni)] <- 2
  
  ## Pre-calculate constant denominator for 3rd way
  mean_beta <- mean(beta)
  prob_deno <- mean_beta^2 * (n - 1)
  
  ## Set up storage for each state
  S <- E <- I <- R <- rep(0, nt)
  
  ## Record initial population counts (Day 1)
  S[1] <- sum(x == 0);E[1] <- sum(x == 1)
  I[1] <- sum(x == 2);R[1] <- sum(x == 3)
  
  ## Loop over the remaining days
  for (i in 2:nt) {
    
    u <- runif(n) ## generate uniform random deviates
    
    x[x == 2 & u < delta] <- 3  ## I -> R with prob delta
    x[x == 1 & u < gamma] <- 2 ## E -> I with prob gamma
    
    ## S -> E (Three ways)
    I_idx <- which(x == 2) ## get infection indices
    
    ## Loop through each infector j.
    for (j in I_idx) {
      h_id <- h[j] ## get the household ID of j
      net <- alink[[j]] ## get the regular network contacts (indices) of j
      
      ## 1st way: household transmission 
      ## find people who are S(0) and in j's household
      S_in_hh_idx <- which(x == 0 & h == h_id)
      ## check if there are any S targets in the household
      if (length(S_in_hh_idx) > 0) {
        ## generate uniform random deviates
        u_hh <- runif(length(S_in_hh_idx)) 
        ## determine who is infected (prob alpha[1]).
        x[S_in_hh_idx[u_hh < alpha[1]]] <- 1 ## S-> E with prob alpha[1]
      }
      
      ## 2nd way: Regular Network Transmission 
      ## Filter out contacts who are still S(0)
      S_net_idx <- net[x[net] == 0] 
      ## check if there are any S targets in the regular network
      if (length(S_net_idx) > 0) {
        ## generate uniform random deviates
        u_net <- runif(length(S_net_idx)) 
        ## determine who is infected (prob alpha[2])
        x[S_net_idx[u_net < alpha[2]]]<- 1 ## S-> E with prob alpha[2]
      }
      
      ## 3rd way: Irrespective Transmission 
      S_all_idx <- which(x == 0)
      
      ## Build exclusion list: j itself, all household members, and all regular contacts
      excluded_idx <- unique(c(j, which(h == h_id), net)) 
      
      ## S_random_indices: S individuals NOT in j's household or regular net
      S_random_idx <- S_all_idx[!(S_all_idx %in% excluded_idx)]
      ## check if there are any S targets
      if (length(S_random_idx) > 0) {
        
        ## Setting P_random
        P_random<- pmin(alpha[3] * nc * beta[S_random_idx]* beta[j] / prob_deno, 1) 
        ## Cap probability at 1
        ## generate unifrom random deviates
        u_random <- runif(length(S_random_idx)) 
        ## determine who is infected (prob P_random)
        x[S_random_idx[u_random < P_random]] <- 1 ## S-> E with prob P_random
      }
    }
    
    ## Store daily results
    S[i] <- sum(x == 0); E[i] <- sum(x == 1)
    I[i] <- sum(x == 2); R[i] <- sum(x == 3)
  }
  
  return(list(S = S, E = E, I = I, R = R))
}





