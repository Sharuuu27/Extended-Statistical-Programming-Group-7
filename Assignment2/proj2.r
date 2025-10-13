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
### To assign the n people to their unique household ID

n <- 1000
hmax <- 5

h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])


## Question 2
### To assign Sociability Parameter (Beta-i)
generate_beta <- function(n = 1000, bmu = 5e-5, bsc = 1e-5) {
  beta <- rgamma(n, shape=bmu/bsc, scale=bsc)
  return(beta)
}  

beta <- generate_beta(n=n)

### To assign regular network relations

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
  ## beta = the transmission rate parameter;
  ## h = household each person belongs to;
  ## alpha = the daily probs I[i] -> S[j]
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## nc = Average contacts; nt = number of days
  ## pinf = initial proportion infected.
  
  ##Initialization 
  set.seed(0)
  n <- length(beta)
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
  S[1] <- n-ni;I[1] <- ni
  
  ## Loop over the remaining days
  for (t in 2:nt) {
    
    u <- runif(n)          ## generate uniform random deviates
    I_idx <- which(x == 2) ## get infection indices
    
    x[x == 2 & u < delta] <- 3 ## I -> R with prob delta
    x[x == 1 & u < gamma] <- 2 ## E -> I with prob gamma
    
    # S -> E (Three ways)
    for (i in I_idx) {  ## Loop through each infector i
      h_id <- h[i]      ## get the household ID of i
      net <- alink[[i]] ## get the regular network contacts (indices) of i
      
      ## 1st way: household transmission
      ## find people who are S(0) and in i's household
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
      ### check if there are any S targets in the regular network
      if (length(S_net_idx) > 0) {
        ## generate uniform random deviates
        u_net <- runif(length(S_net_idx))
        ## determine who is infected (prob alpha[2])
        x[S_net_idx[u_net < alpha[2]]] <- 1 ## S-> E with prob alpha[2]
      }
      
      ## 3rd way: Irrespective Transmission
      S_all_idx <- which(x == 0)
      ## check if there are any S targets
      if (length(S_all_idx) > 0) {
        ## Setting P_random
        ## Use pmin to cap probability at 1
        P_random <- pmin(alpha[3] * nc * beta[i] * beta[S_all_idx] / prob_deno, 1)
        ## generate uniform random deviates
        u_random <- runif(length(S_all_idx))
        ## determine who is infected (prob P_random)
        x[S_all_idx[u_random < P_random]] <- 1 ## S-> E with prob P_random
      }
    }
    
    ## Store daily results
    S[t] <- sum(x == 0); E[t] <- sum(x == 1)
    I[t] <- sum(x == 2); R[t] <- sum(x == 3)
  }
  
  return(list(S = S, E = E, I = I, R = R))
}


## Question4
## A function to nicely plot the dynamics of the simulated population states
plot_graphs <- function(alphas, betas, labels) {
  ## alphas = the daily probs I[i] -> S[j] in matrix format;
  ## betas = the transmission rate parameter in matrix format;
  ## labels = title of each parameter settings
  
  ## setting for plotting multiple graphs
  ncol <- length(alphas[1,])    ## set number of graphs to plot
  par(mfcol=c(2,ncol/2),        ## layout setting of the plot
      mar=c(4,4,2,1),           ## margin settings
      oma=c(.2,.2,.2,.2),       ## outer margin settings
      family="sans",            ## font style settings (Arial)
      cex.axis=.6, cex.lab=.9)  ## font size settings
  
  for (i in 1:ncol){
    ## draw y axis labels only on the first graph
    if (i<=2) scY <- "N" else scY <- ""      ## draw Y-axis label on the left graphs
    if (i%%2==0) scX <- "day" else scX <- "" ## draw X-asix label on the bottom graphs
    
    ## Simulate SEIR model with nseir function
    epi <- nseir(beta=betas[,i],h,alink,alpha=alphas[,i])
    
    ## plot simulated data
    plot(epi$S,ylim=c(0,max(epi$S)),xlab=scX,ylab=scY,
         main=labels[i],
         type="l",cex=.5,cex.main=.8,lwd=2)     ## S (black)
    lines(epi$E,col=4,type="l",cex=.5,lwd=2)    ## E (blue)
    lines(epi$I,col=2,type="l",cex=.5,lwd=2)    ## I (red)
    lines(epi$R,col=3,type="l",cex=.5,lwd=2)    ## I (red)
    legend(x="right",y="center",                ## legend
           legend=c("S", "E", "I", "R"), col=c(1, 4, 2, 3),
           lty=1, bty="n",cex=.7,lwd=1.5)
  }
}


## Question5
## Compare 4 scenarios and plot
### parameter inputs
alphas <- array(c(c(.1,.01,.01),
                  c(0,0,.04),
                  c(.1,.01,.01),
                  c(0,0,.04)),
                dim=c(3,4))

set.seed(0)
beta_u <- runif(n)
beta_mean <- rep(mean(beta_u),n)
betas <- array(c(beta_u,beta_u,beta_mean,beta_mean),
               dim=c(n,4))

labels <- c("Default parameters\n(αh=.1,αc=αr=.01)",
            "Without household and regular network\n(αh=αc=0,αr=.04)",
            "Constant β\n(αh=.1,αc=αr=.01)",
            "Constant β & random mixing\n(αh=αc=0,αr=.04)")

## Plot graphs
plot_graphs(alphas, betas, labels)
