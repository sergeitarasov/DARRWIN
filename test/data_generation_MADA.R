'
Script developed by: Thomas Merrien
last update: 1/06/2022
File name: data_generation_MADA.R

Description:
Using real Madagascar data needed to test the pruning algorithm of our
ancestral range reconstruction model
'

library(dplyr)
library(ape)
library(geiger)

# Importation of the Geographic data -------------------------------------------

all_species <- read.csv("test/all_species_in_grid.csv")
mada_grid <- read.csv("test/grid_forest.csv")

#SOI = Species of interest

all_species$BINOMIAL <- gsub(" ", "_", all_species$BINOMIAL)

SOI <- unique(all_species$BINOMIAL)

SOI <- SOI[grep("Microcebus",SOI)]

# Generation of a the phylogenetic tree -------------------------------------

tree <- read.nexus(file = "test/Phylo/fbd421agerange.tre")
small_tree <- drop.tip(tree, tree$tip.label[grepl("Microcebus",tree$tip.label)==FALSE])

#To drop
drop <- SOI[grep(FALSE, SOI %in% small_tree$tip.label)]
drop <- c(drop, small_tree$tip.label[grep(FALSE, small_tree$tip.label %in% SOI)])
small_tree <- drop.tip(small_tree, drop)
SOI <- SOI[grep(FALSE, SOI %in% drop)]

n <- length(small_tree$tip.label)

# Grid map parameters ----------------------------------------------------------

w = length(mada_grid$fid)
h = 1 #height: number of rows

cell_id <- sort(mada_grid$id)

# Generation of the biogeographic database -------------------------------------

Y <- list(matrix(0, nrow = 1, ncol = w),matrix(0, nrow = 1, ncol = w))
Y <- rep(Y,small_tree$Nnode+n)

for (i in 1:n){
  occup_cell <- c()
  occup_cell <- all_species$id[grep(TRUE, all_species$BINOMIAL==small_tree$tip.label[i])]
  mod_occup_cell <- rep(0, length(occup_cell))
  for (j in 1:length(occup_cell)){
    mod_occup_cell[j] <- grep(occup_cell[j], cell_id)
  }

  Y[[(i)]][mod_occup_cell] <- 1
  #Y[[(2*i)]] <- 1-Y[[(2*i)-1]]
}



# Generation of the grid map landscape structure -------------------------------

geo_types <- c("land", "forest")

R <- rep("land",w)
for (i in 1:length(R)){
  if (mada_grid$Forest[i]==1){
    R[i]<-"forest"
  }
}

# Generation of the transition matrix ------------------------------------------

#Area types
K = length(geo_types)
Q <- matrix(0, nrow = K, ncol = K)
colnames(Q) <- geo_types
rownames(Q) <- geo_types

Q[1]=0.2
Q[2]=0.2
Q[3]=0.200
Q[4]=0.2

# Delta time for the analysis --------------------------------------------------

dt <- 0.005

# Generation of the square occupation probability matrices ---------------------

Pi0 <- list(matrix(0.5, nrow = h, ncol = w),matrix(0, nrow = h, ncol = w))
Pi1 <- list(matrix(0.5, nrow = h, ncol = w),matrix(0, nrow = h, ncol = w))

#Adjacency matrix receive 1 if cells are neighbor, 0 otherwise
#A <- matrix(0, nrow=h*w, ncol = h*w)
A <- list.neighbor2(mada_grid)
diag(A) <- 0

# Create the list of neighbor --------------------------------------------------
#neighbor <- list.neighbor(map)

neighbor <- list()

for (i in 1:(h*w)){
  neighbor[[i]] <- c(grep(1, A[i,]),0,0,0)
}

# Get the phylogenetic tree structure ------------------------------------------
brs <- branching_time(small_tree)
brs <- data.table::as.data.table(brs)
colnames(brs) <- c("V1", "V2", "V3", "V4", "V5")
brs <- dplyr::arrange(brs, desc(V4))
#increase resolution (need to lower the R parameter if done)
brs[,3:5] <- brs[,3:5]*1

# Q matrix design --------------------------------------------------------------


#Rate matrix at the grid level
infQ <- matrix(0, nrow=h*w, ncol = h*w)
infQdt <- matrix(0, nrow=h*w, ncol = h*w)
matQ <- matrix(0, nrow=h*w, ncol = h*w)

for(i in 1:(h*w)){
  for (j in grep(1, A[i,])){
    matQ[i,j] <- Q[R[i],R[j]]
    infQ[i,j] <- Q[R[i],R[j]]
    infQdt[i,j] <- Q[R[i],R[j]]*dt
  }
}


for (i in 1:(h*w)){
  matQ[i,i] <- -sum(matQ[i,])
  infQ[i,i] <- 1 - sum(infQ[i,])
  infQdt[i,i] <- 1 - sum(infQdt[i,])
}

matQ <- t(matQ)
