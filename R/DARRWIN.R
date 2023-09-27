'
Script developed by: Thomas Merrien
last update: 31/08/2023
File name: DARRWIN.R

Description: wrapper to launch DARRWIN
'


library(expm)
library(corHMM)
library(nloptr)
library(diversitree)

source("R/list_neighbor.R")
source("R/branching_time.R")
source("R/reps_per_period.R")
source("R/pruning_algo.R")
source("R/ancestral_state_probability.R")
source("R/distribution_auxiliary.R")
source("R/parameter_optimization.R")
source("R/center_of_mass.R")
source("R/plot_result.R")
source("R/raster_discretization.R")

#DARRWIN(tree, "all_species_in_grid.csv", "grid_forest.csv", "test", c("Land", "Forest"))

#'
#'
#' @param envt_type list of the different environmental types
#' @param lbr lower bound of the elements of the rate matrix
#' @param ubr upper bound of the elements of the rate matrix
#' @param branch_rate_multiplier TRUE if you want to have specific branch rate multiplier, FALSE if you want to have the same branch rate everywhere


DARRWIN <- function(tree, sp_dist, grid, directory, envt_type, lbr=0, ubr=100, branch_rate_multiplier){

  #Step 1: Declinaison of the tree under appropriate format --------------------

  n <- length(tree$tip.label)
  brs <- branching_time(tree)
  brs <- data.table::as.data.table(brs)
  colnames(brs) <- c("V1", "V2", "V3", "V4", "V5")
  brs <- dplyr::arrange(brs, desc(V4))

  #Step 2: import the geographic data ------------------------------------------

  species_dist <- read.csv(paste(directory,sp_dist, sep="/"))
  grid_map <- read.csv(paste(directory,grid, sep="/"))

  w = length(grid_map$fid)
  h = 1 #height: number of rows

  cell_id <- sort(grid_map$id)

  #Step 3: Formating the distribution data for algo ----------------------------
  Y <- list(matrix(0, nrow = 1, ncol = w),matrix(0, nrow = 1, ncol = w))
  Y <- rep(Y,tree$Nnode+n)


  species_dist$BINOMIAL <- gsub(" ", "_", species_dist$BINOMIAL)

  for (i in 1:n){
    occup_cell <- c()
    occup_cell <- species_dist$id[grep(TRUE, species_dist$BINOMIAL==tree$tip.label[i])]
    mod_occup_cell <- rep(0, length(occup_cell))
    for (j in 1:length(occup_cell)){
      mod_occup_cell[j] <- grep(occup_cell[j], cell_id)
    }

    Y[[i]][mod_occup_cell] <- 1
  }

  for (i in 1:n){
    Y[[i]] <- Y[[i]]/sum(Y[[i]])
  }

  #Step 4: Adding the grid landscape structure ---------------------------------

  R <- rep(envt_type[1],w)

  for (i in envt_type[-1]){
    index <- grep(TRUE, i==colnames(grid_map)) #find the index of the column containing data about the environment of interest
    for (j in 1:length(R)){
      if (grid_map[j,index]==1){
        R[j]<-i
      }
    }
  }

  #Step 5: Creation of the adjacency matrix and of the neighboring list
  A <- list.neighbor2(grid_map)
  diag(A) <- 0

  neighbor <- list()


  for (i in 1:(h*w)){
    neighbor[[i]] <- c(grep(1, A[i,]),0,0,0)
  }

  #Step 6: Parameter optimization for the Q matrix
  # We will do a first search with a larger threshold for validation of results
  # to speed up the research.
  #Then we do a second thinner search for optimum around the optima found on the
  #first step.
  #These 2 steps are done for 4 different starts and we then select the best one
  #based on the likelihood

  para <- list()
  Like <- c()

  for (q in 1:10){

    x0 <- runif(length(envt_type)^2, lb, ub)

    opts <- list("algorithm" = "NLOPT_GN_DIRECT_L_RAND",
                 "xtol_rel" = 1e-2,
                 maxeval = 1000)

    params <- nloptr(x0 = x0,
                     lb = rep(lbr, length(envt_type)^2),
                     ub = rep(ubr, length(envt_type)^2),
                     eval_f = optimization_function,
                     opts = opts)


    opts <- list("algorithm" = "NLOPT_GN_DIRECT_L_RAND",
                 "xtol_rel" = 1e-5,
                 maxeval = 1000)

    params <- nloptr(x0 = params$solution,
                     lb = rep(lb, length(envt_type)^2),
                     ub = rep(ub, length(envt_type)^2),
                     eval_f = optimization_function,
                     opts = opts)

    para[[q]] <- params$solution
    Like[q] <- params$objective

  }

  select <- grep(min(Like), Like)[1] #index of the minimum negative log likelihood solution


  #Step 7: Construction of the Q matrix

  Q <- matrix(para[select], nrow = length(envt_type), ncol = length(envt_type))
  colnames(Q) <- envt_type
  rownames(Q) <- envt_type

  print(Q)

  matQ <- matrix(0, nrow=h*w, ncol = h*w)

  for(i in 1:(h*w)){
    for (j in grep(1, A[i,])){
      matQ[i,j] <- Q[R[i],R[j]]
    }
  }

  for (i in 1:(h*w)){
    matQ[i,i] <- -sum(matQ[i,])
  }

  matQ <- t(matQ)

  #Step 8: Running Likelihood analysis ------------------------------------------
  Darrwin_LNL <- tree_likelihood(Y, brs, matQ, rep(1,length(Y)))

  #Step 9: Running marginal reconstruction analysis -----------------------------
  reconstruction <- ancestral_state_probability(n, brs, neighbor, Darrwin_LNL$Y, matQ, Darrwin_LNL$Delta, init_proba = "equiproba")

  for (i in 1:length(reconstruction)){

    reconstruction[[i]] <- reconstruction[[i]]/max(reconstruction[[i]])

  }
  #Step 10: Displaying results

  plot_phylomap(grid_map, cell_id, tree, reconstruction)

  return(Q)

}
