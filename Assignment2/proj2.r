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


