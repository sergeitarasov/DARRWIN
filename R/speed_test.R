'
Script developed by: Thomas Merrien
creation: 05/05/2023
last update: 05/05/2023
File name: speed_test.R

Description:
Script to compare speed time of different function. Aim is to reduce computation
time when this is possible
'

my_time <- c()
R_time <- c()


for (i in 1:500){

  x <- rgamma(100000,2,2)

  startTime <- Sys.time()

  for (i in 1:100000){
    dens_gamma(x[i],2,2)
  }

  endTime <- Sys.time()

  my_time <- c(my_time, endTime-startTime)

  startTime <- Sys.time()

  for (i in 1:100000){
    dgamma(x[i],2,2)
  }

  endTime <- Sys.time()

  R_time <- c(R_time, endTime-startTime)

}

all_time <- data.frame(time=rep(0,2000), method=rep("R", 2000))
all_time$time[1:1000] <- my_time
all_time$method[1:1000] <- "dens_gamma"
all_time$time[1001:2000] <- R_time
all_time$method[1001:2000] <- "dgamma"

library(ggplot2)

p <- ggplot() + geom_boxplot(data=all_time, mapping = aes(y=time, x=method))
