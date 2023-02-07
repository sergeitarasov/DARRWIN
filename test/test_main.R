'
Script developed by: Thomas Merrien
last update: 26/04/2022
File name: test_main.R

Description:
Testing the different function using the data generated in the data_generation.R
file
'

# Call the different functions -------------------------------------------------
source("R/list_neighbor.R")
source("R/branching_time.R")
source("R/reps_per_period.R")
source("R/pruning_algo.R")
source("R/ancestral_state_probability.R")
library(expm)
library(corHMM)

# Using pruning algortihm to reconstruct ancestral ranges ----------------------
test_3cell <- exp_pruning_simple(brs,Y,neighbor, matQ, dt)

test_lemur <- matrix_pruning(brs, Y, neighbor, infQdt, dt)
test_lemur <- exp_pruning_simple(brs, Y, neighbor, matQ, dt)

test2 <- matrix_pruning(brs, Y, neighbor, infQdt, dt)
test3 <- exp_pruning(brs, Y, neighbor, matQ, dt)

Ptest <- ancestral_state_probability(n, brs, neighbor, test_lemur, infQdt, dt, init_proba = "equiproba")
Ptest <- ancestral_state_probability_exp_simple(n, brs, neighbor, test_lemur, matQ, dt, init_proba = "equiproba")

Ptest2 <- ancestral_state_probability_matrix(n, brs, neighbor, test3, matQ, dt, init_proba = "equiproba")
Ptest_3cell <- ancestral_state_probability_exp_simple(n, brs, neighbor, test_3cell, matQ, dt, init_proba = "equiproba")

data_test <- data.frame(matrix(NA, 3,2))
data_test[,1] <- tree$tip.label
data_test[1,2] <- 1
data_test[2,2] <- 2
data_test[3,2] <- 3

colnames(data_test) <- c("Genus_sp", "T1")
rownames(data_test) <- data_test$Genus_sp
MK_3state <- corHMM(tree,data_test, rate.cat = 1)

# # one way to get the parameters from your corHMM object in the correct order
p <- sapply(1:max(MK_3state$index.mat, na.rm = TRUE), function(x)
  na.omit(c(MK_3state$solution))[na.omit(c(MK_3state$index.mat) == x)][1])

# using custom params
states_1 <- ancRECON(phy = tree, data = MK_3state$data, p = c(0.1,0.0,0.1,0.1,0.0,0.1), method = "marginal",
                     rate.cat <- MK_3state$rate.cat, ntraits = 1, rate.mat = MK_3state$index.mat,
                     root.p = MK_3state$root.p)

# importing raster data --------------------------------------------------------
library(raster)
library(RcppArmadillo)
library(ggplot2)
library(sf)
library(ggmap) #to import and create maps
library(tmap) #to make static and interactive map

lemur_distrib <- st_read("test/data_0.shp")
#lemur_distrib <- shapefile("test/data_0.shp")

mada <- c(left = 43, bottom = -26, right = 51, top = -11.5)
map <- get_stamenmap(mada, zoom = 5, maptype = "toner-background")
ggmap(map)

#trying to plot some distrib
ggplot() +
  geom_sf(data = lemur_distrib[1,]) +
  coord_sf()

ggmap(map) + geom_raster(lemur_distrib[1,]) + coord_cartesian()

