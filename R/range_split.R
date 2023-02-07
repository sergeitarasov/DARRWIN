'
Script developed by: Thomas Merrien
last update: 07/11/2022
File name: range_split.R

Description:
Testing models for splitting range when speciation happens
'

library(ggplot2)

#### 1D example ####

mean_st1 <- 5
mean_st2 <- 1
sd_st1 <- 0.5
sd_st2 <- 0.5

sigma_mean <- 1
sigma_sd <- 1

delta_t <- 2

brownian_motion <- function(x,sigma,time){
  x <- x + rnorm(1,0,sigma)*time
}

new_mean_st1 <- brownian_motion(mean_st1, sigma_mean, delta_t)
new_sd_st1 <- brownian_motion(sd_st1, sigma_sd, delta_t)

x_dir <- new_mean_st1 - mean_st1

ci95bot <- new_mean_st1 - 1.96*new_sd_st1
ci95up <- new_mean_st1 + 1.96*new_sd_st1

percentage_inherited <- 0.2
dist_interval <- ci95up-ci95bot
mean_daughter1 <- ci95bot + 0.2*dist_interval
mean_daughter2 <- new_mean_st1

ggplot(data = data.frame(x = c(-3, 12)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = mean_daughter1, sd = sd_st1), col = "blue") +
  stat_function(fun = dnorm, n = 101, args = list(mean = mean_st1, sd = sd_st1), col = "red") +
  stat_function(fun = dnorm, n = 101, args = list(mean = new_mean_st1, sd = new_sd_st1), col = "black") +
  ylab("") +
  scale_y_continuous(breaks = NULL)

#### 2D example ####

TwoDexample <- matrix(1:(9), ncol = 3, nrow = 3)
Aexample <- matrix(0,9,9)
Aexample[1,2] <- 1
Aexample[1,4] <- 1

Aexample[2,1] <- 1
Aexample[2,3] <- 1
Aexample[2,5] <- 1

Aexample[3,2] <- 1
Aexample[3,6] <- 1

Aexample[4,1] <- 1
Aexample[4,7] <- 1
Aexample[4,5] <- 1

Aexample[5,2] <- 1
Aexample[5,4] <- 1
Aexample[5,8] <- 1
Aexample[5,6] <- 1

Aexample[6,3] <- 1
Aexample[6,5] <- 1
Aexample[6,9] <- 1

Aexample[7,4] <- 1
Aexample[7,8] <- 1

Aexample[8,7] <- 1
Aexample[8,5] <- 1
Aexample[8,9] <- 1

Aexample[9,7] <- 1
Aexample[9,8] <- 1

Qexample <- Aexample
Qexample[Qexample==1] <- 0.1
for (i in 1:9){
  Qexample[i,i] <- -sum(Qexample[,i])
}

ancestor <- rep(0,9)
ancestor[1] <- 1
delta_t <- 2

#range evolution

at_speciation <- ancestor%*%expm(Qexample*delta_t)
at_speciation <- at_speciation/sum(at_speciation)

#split potential range
percentage_inherited <- 0.2

