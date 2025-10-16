# Team Members: 
# 1. Shivasharini Govindasamy (Sharuuu27) - S2766392
# 2. Yu-Hsuan Hung (YuHsuan07) - S2793274
# 3. Yasuhiro Hara (hiroh-git) - S2826059

# Team Contributions:
# Shivasharini - Developed the household vector h generation (Q1) 
                #and the get.net function for the regular contact network (Q2)
# Yu-Hsuan - Implemented the core SEIR simulation function, nseir, 
            #which includes the social structure elements (Q3)
# Yasuhiro- Developed the function for plotting the simulation dynamics (Q4) 
            #and performed the comparative analysis of the four model scenarios (Q5)

# Github repo link:
# https://github.com/Sharuuu27/Extended-Statistical-Programming-Group-7-Assignment-Submission

## OVERVIEW: STOCHASTIC SEIR SIMULATION WITH SOCIAL STRUCTURE
# This script implements a stochastic SEIR epidemiological model.
# The model incorporates explicit social structure by simulating contact based on 
# household membership, individual social networks, and random mixing.

## CODE STRUCTURE:
# 1. Setup (N, hmax) and Household Assignment
# 2. Sociability Parameter generation and Network Generation function
# 3. SEIR simulation function (nseir)
# 4. Plotting function (plot_graphs)
# 5. Scenario Comparison and Plotting
# 6. Interpretation of the plots

set.seed(0)

## Step 1: Household Assignment
### To assign the n people to their unique household ID
n <- 1000 ## Define population size N
hmax <- 5  ## Max household size

## Create household vector 'h' such that household sizes are 
## uniformly distributed between 1 and hmax.
## This one-liner uses: 
## 1. sampling sizes (sample(1:hmax, n, replace = TRUE))
## 2. repeating household IDs (rep(1:n, times = ...))
## 3. sampling the first N elements and shuffling (sample(...)[1:n])
h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])


## Step 2: Assignment of Sociability Parameter (beta_i) and Network Generation.
### To assign Sociability Parameter (beta_i)
# Generates the vector of sociability parameters (beta_i) for each person.
beta <- runif(n)

get.net <- function(beta, h, nc = 15) {
  ## Function to generate the regular contact network.
  ## The network excludes household members
  n <- length(beta) 
  mean_beta <- mean(beta)
  ## Calculate the denominator used in the link probability formula
  prob_denominator <- mean_beta^2 * (n - 1)
  # creating network list of all n people
  network <- replicate(n, integer(0), simplify = FALSE) 
  # creating all possible pairs
  row1 <- c(); row2 <- c()
  row1 <- unlist(sapply(1:n, function(i) rep(i, n-i+1)))
  row2 <- unlist(sapply(1:(n-1), function(i) (i+1):n))
  unique_pairs <- array(c(row1, row2),dim=c(length(row1), 2))
  unique_pairs <- t(unique_pairs)
  # to identify which pairs are household members
  household_member <- h[unique_pairs[1, ]] == h[unique_pairs[2, ]]
  
  # to store the pairs which are NOT household members
  non_household <- unique_pairs[, !household_member]
  # contact network probability and cap at 1
  link <- ((nc * beta[non_household[1, ]] * beta[non_household[2, ]]) 
           / prob_denominator)
  
  # generate uniform random numbers to compare with contact network probability
  random_u <- runif(ncol(non_household))
  regular_contact <- random_u < link
  # to store the regular contacts only
  network_list <- non_household[, regular_contact]
  # To add regular contacts to each person's list
  for (j in 1:ncol(network_list)) {
    person_1 <- network_list[1, j] 
    person_2 <- network_list[2, j] 
    
    
    network[[person_1]] <- c(network[[person_1]], person_2)
    network[[person_2]] <- c(network[[person_2]], person_1)
  }
  
  return(network)
  
}

alink <- get.net(beta, h)


## Step 3: Implementation of the SEIR Simulation Model
## SEIR simulation model with social structure
nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, 
                  gamma = .4, nc = 15,nt = 100, pinf = .005) {
  ## beta = socialibility parameter;
  ## h = household each person belongs to;
  ## alpha = the daily probs I[i] -> S[j]
  ## Note: Each alpha stands for 
  ## - household transmission prob
  ## - contact network transmission prob
  ## - random mixing base prob
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## nc = Average contacts; nt = number of days
  ## pinf = initial proportion infected.
  
  ##Initialization 
  set.seed(0)
  n <- length(beta)
  x <- rep(0, n)
  ni <- round(n * pinf)
  
  ## Randomly set the initial number of infected individuals (I=2)
  x[1:ni] <- 2
  
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
    
    ## Fixed Probability Transitions (I -> R & E -> I)
    x[x == 2 & u < delta] <- 3 ## I -> R with prob delta
    x[x == 1 & u < gamma] <- 2 ## E -> I with prob gamma
    
    # S -> E Transitions (Three ways)
    for (i in I_idx) {  ## Loop through each infector i
      h_id <- h[i]      ## get the household ID of i
      net_idx <- alink[[i]] ## get the regular network contacts (indices) of i
      
      ## 1st way: household transmission
      ## Find people who are S(0) and in i's household
      S_same_hh_idx <- which(x == 0 & h == h_id)
      S_same_hh_idx_len <- length(S_same_hh_idx)
      ## check if there are any S targets in the household
      if (S_same_hh_idx_len > 0) {
        ## generate uniform random deviates
        u_hh <- runif(S_same_hh_idx_len)
        
        ## determine who is infected (prob alpha[1]).
        x[S_same_hh_idx[u_hh < alpha[1]]] <- 1 ## S-> E with prob alpha[1]
      }
      
      ## 2nd way: Regular Network Transmission
      ## Filter out contacts who are still S(0)
      S_net_idx <- net_idx[x[net_idx] == 0]
      S_net_idx_len <- length(S_net_idx)
      ### check if there are any S targets in the regular network
      if (S_net_idx_len > 0) {
        ## generate uniform random deviates
        u_net <- runif(S_net_idx_len)
        ## determine who is infected (prob alpha[2])
        x[S_net_idx[u_net < alpha[2]]] <- 1 ## S-> E with prob alpha[2]
      }
      
      ## 3rd way: Irrespective Transmission 
      # Find all susceptible individuals
      S_noinf_idx <- which(x == 0)
      S_noinf_idx_len <- length(S_noinf_idx)
      ## check if there are any S targets
      if (S_noinf_idx_len > 0) {
        ## Setting P_random
        ## Note: Use pmin to cap probability at 1
        P_random <- alpha[3] * nc * beta[i] * beta[S_noinf_idx] / prob_deno
        ## generate uniform random deviates
        u_random <- runif(S_noinf_idx_len)
        ## determine who is infected (prob P_random)
        x[S_noinf_idx[u_random < P_random]] <- 1 ## S-> E with prob P_random
      }
    }
    
    ## Store daily results by counting individuals in each state
    S[t] <- sum(x == 0); E[t] <- sum(x == 1)
    I[t] <- sum(x == 2); R[t] <- sum(x == 3)
  }
  
  return(list(S = S, E = E, I = I, R = R, t=1:nt))
}


## Step 4: Function for plotting simulation results
## A function to nicely plot the dynamics of the simulated population states
## Uses base R plotting functions (plot, lines) and graphical parameters (par)
plot_graphs <- function(alphas, betas, labels, h, alink) {
  ## alphas = the daily probs I[i] -> S[j] in matrix format;
  ## betas = the transmission rate parameter in matrix format;
  ## labels = title of each parameter settings
  
  ## setting for plotting multiple graphs
  ncol <- length(alphas[1,])    ## set number of graphs to plot
  
  # Set up graphical parameters for multi-panel plot (2 rows, 2 columns if ncol=4).
  # mfcol fills plots column by column.
  par(mfcol=c(2,ncol/2),        ## layout setting of the plot
      mar=c(4,4,2,1),           ## margin settings
      oma=c(.2,.2,.2,.2),       ## outer margin settings
      family="sans",            ## font style settings (Arial)
      cex.axis=.6, cex.lab=.9)  ## font size settings
  
  for (i in 1:ncol){
    ## Determine axis labels based on panel position
    ## Y label 'N' is drawn only on the first column (i=1, 2).
    if (i<=2) scY <- "N" else scY <- ""      ## draw Y-axis label on the left graphs
    if (i%%2==0) scX <- "day" else scX <- "" ## draw X-asix label on the bottom graphs
    
    ## Simulate SEIR model with nseir function
    epi <- nseir(beta=betas[,i],h,alink,alpha=alphas[,i])
    
    ## plot simulated data
    plot(x=epi$t, y=epi$S,
         ylim=c(0,max(epi$S)),xlab=scX,ylab=scY,
         main=labels[i],
         type="l",cex=.5,cex.main=.8,lwd=2) ## S (black)
    lines(x=epi$t, y=epi$E,
          col=4,type="l",cex=.5,lwd=2)      ## E (blue)
    lines(x=epi$t, y=epi$I,
          col=2,type="l",cex=.5,lwd=2)      ## I (red)
    lines(x=epi$t, y=epi$R,
          col=3,type="l",cex=.5,lwd=2)      ## I (red)
    legend(x="right",y="center",            ## legend
           legend=c("S", "E", "I", "R"), col=c(1, 4, 2, 3),
           lty=1, bty="n",cex=.7,lwd=1.5)
  }
}


## Step 5: Compare 4 scenarios and plot
## Prepare parameters for the 4 required scenarios: 
## 1. Full model, varied beta (Default)
## 2. Random mixing, varied beta
## 3. Full model, constant beta
## 4. Random mixing, constant beta
alphas <- array(c(c(.1,.01,.01),
                  c(0,0,.04),
                  c(.1,.01,.01),
                  c(0,0,.04)),
                dim=c(3,4))

set.seed(0)

## Generate constant beta set to the average of beta_u for every element
beta_mean <- rep(mean(beta),n)
betas <- array(c(beta,beta,beta_mean,beta_mean),
               dim=c(n,4))

## Labels describing each scenario for plotting
labels <- c("Default parameters\n(αh=.1,αc=αr=.01)",
            "Without household and regular network\n(αh=αc=0,αr=.04)",
            "Constant β\n(αh=.1,αc=αr=.01)",
            "Constant β & random mixing\n(αh=αc=0,αr=.04)")

## Plot graphs
plot_graphs(alphas, betas, labels, h=h, alink=alink)


## The apparent effect of the household and network structure, relative to random mixing

### Overview:
# The epidemic spreads faster when the model is computed without the household 
# and regular network structure (bottom plots). 
# In addition, it has a bigger outbreak when sociability parameter beta is constant.

# When assuming random mixing only, not only the social structure, but also 
# everyone has a higher probability to infect others in the population (due to
# bigger numerator of the daily probability of i infecting j), which may
# accelerating the expand of the epidemic.
# This can be supported with the fact that the less number of people remaining
# in state S(susceptible) at the end of the simulation compared to other 2 scenarios. 

# Furthermore, differences in the settings of the sociability parameter may
# suggest that individuals with a relatively low level of sociability may suppress
# the infection with a variable sociability parameter, resulting in slower 
# overall expansion compared to a scenario in which everyone has an average level
# of sociability.


### Details of each graph
## 1. Default parameters
# (alpha_h =.1, alpha_c = alpha_r =.01; beta varied) (left top)
#
# The number of infected (I, red line) and exposed (E, blue line) peaks 
# around days 38 and 30, respectively. The susceptible (S, black line) 
# population declines, reaching a minimum of approximately 200 on day 50. 
# The recovered (R, green line) population increases and plateaus near 800.


## 2. Without household and regular network
# (alpha_h = alpha_c = 0, alpha_r =.04; beta varied) (left bottom)
# 
# The spread and resolution are faster than in the default scenario. The 
# peaks of states I and S are higher at approximately 200 and 150, respectively. 
# The peaks of infected and exposed also occur earlier, around days 25 and 
# 20 respectively, suggesting that the epidemic spreads and resolves 
# relatively quickly. The population of the state S drops sharply, reaching a 
# minimum of approximately 200 on day 40. The state's R population plateaus 
# at approximately 800. The high transmission through random mixing (alpha_r =.04) 
# likely drives this rapid spread.


## 3. Constant beta
# (alpha_h =.1, alpha_c = alpha_r =.01; beta constant) (right top)
#
#  The overall shape of the epidemic is very similar to the default parameters, 
# but the minimum value of state S's population is lower, approximately 150, 
# and the maximum value of R's population is higher, approximately 850. 
# This suggests that the homogeneity of beta slightly increased the spread.


## 4. Constant beta & random mixing
# (alpha_h = alpha_c = 0, alpha_r =.04; beta constant) (right bottom)
#
#  The peaks of the I and E states are higher than in scenarios 1 and 3, 
# at around 200 and 150 respectively. These peaks occur earlier, 
# on days 35 and 25 respectively.
# The minimum value of state S's population is lower, approximately 100,
# and the maximum value of R's population is higher, approximately 900. 
# Similar to 3rd scenario, using a constant beta over a varied beta 
# slightly increases the size of the epidemic.
