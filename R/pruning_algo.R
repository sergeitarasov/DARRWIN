'
Script developed by: Thomas Merrien
last update: 09/05/2022
File name: pruning_algo.R
'

#' Likelihood function for one branch with a branch rate multiplier
#'
#' Function to compute the likelihood of one branch with a rate multiplier
#'
#' @param br State/Likelihood at the tip of the branch, given as a list of the probability of each state (list)
#' @param dt Length of the branch (float)
#' @param delta branch specific rate multiplier (float)
#' @param Q transition matrix between the states (matrix)
#'
#'
#' @return the likelihood of one branch as a list of the likelihood of each given states (list)
#'

branch_likelihood_rate <- function(br,dt,delta,Q){

  Ln_br <- br%*%t(expm(Q*dt*delta))

  return (Ln_br)

}


#' Likelihood function for one branch
#'
#' Function to compute the likelihood of one branch
#'
#' @param br State/Likelihood at the tip of the branch, given as a list of the probability of each state (list)
#' @param dt Length of the branch (float)
#' @param delta branch specific rate multiplier (float)
#' @param Q transition matrix between the states (matrix)
#'
#'
#' @return the likelihood of one branch as a list of the likelihood of each given states (list)
#'


branch_likelihood <- function(br,dt,delta,Q){

  Ln_br <- br%*%t(expm(Q*dt*delta))

  return (Ln_br)

}


#' Likelihood function at a node
#'
#' Function to compute the likelihood at a node
#'
#' @param Ln_br1 the likelihood of the first branch as a list of the likelihood of each given states (list)
#' @param Ln_br2 the likelihood of the second branch as a list of the likelihood of each given states (list)
#'
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list)
#'

node_likelihood <- function(Ln_br1, Ln_br2){

  Ln_node <- Ln_br1 * Ln_br2

  return (Ln_node)

}


#' Likelihood function for the whole phylogenetic tree
#'
#' Function to compute the likelihood of a phylogenetic tree and return the scaling parameters as
#'  well as the value of the likelihood for every different states of every nodes.
#'
#' @param Y Probability of presence in each cell across the tree for each species (list of arrays)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Q transition matrix between the states (matrix)
#' @param m branch specific rates vector in the same order as the brs file (array)
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list of arrays)
#'

tree_likelihood <- function(Y,brs,Q, Delta){

  LnL <- 0
  if (length(brs$V1)>50){
    scaling_time <- seq(50,length(brs$V1),50) #We will save the timing when to perform likelihood scaling to avoid underflow
  } else {
    scaling_time <- 50
  }
  Lm <- rep(1, max(brs$V1)) #Here we store the scaling parameters across the tree

  nb <- 0

  for (t in (unique(brs$V1))){ #We will go through each branches
    nb <- nb+1 #We check the number of iteration

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    Ln_br1 <- branch_likelihood(br = Y[[(brs[[index[1],2]])]], dt = brs$V3[index[1]], Q = Q, delta = m[index[1]])
    Ln_br2 <- branch_likelihood(br = Y[[(brs[[index[2],2]])]], dt = brs$V3[index[2]], Q = Q, delta = m[index[2]])

    Y[[t]] <- node_likelihood(Ln_br1 = Ln_br1, Ln_br2 = Ln_br2)

    if (nb %in% scaling_time){ #every 50 nodes we will perform a scaling event

      Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
      Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states

    }

  }

  t <- unique(brs$V1)[length(unique(brs$V1))]

  LnL <- log(prod(Y[[t]]*(1/length(Y[[t]])))) + sum(log(Lm))

  return(list (Y = Y, Lm = Lm , LnL = LnL, scaling_time = scaling_time, Delta = m))

}


#' Likelihood function for the whole phylogenetic tree with gamma distributed branch rate multiplier
#'
#' Function to compute the likelihood of a phylogenetic tree and return the scaling parameters as
#'  well as the value of the likelihood for every different states of every nodes.
#'
#' @param Y Probability of presence in each cell across the tree for each species (list of arrays)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Q transition matrix between the states (matrix)
#' @param m vector of branch rates multiplier organized in the same order as the branches in the brs file (array)
#' @param alpha alpha parameter of the gamma distribution (float)
#' @param beta beta parameter of the gamma distribution (float)
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list of arrays)
#'

tree_likelihood_gamma <- function(Y,brs,Q, m, alpha, beta){

  LnL <- 0
  if (length(brs$V1)>50){
    scaling_time <- seq(50,length(brs$V1),50) #We will save the timing when to perform likelihood scaling to avoid underflow
  } else {
    scaling_time <- 50
  }
  Lm <- rep(1, max(brs$V1)) #Here we store the scaling parameters across the tree

  nb <- 0

  for (t in (unique(brs$V1))){ #We will go through each branches
    nb <- nb+1 #We check the number of iteration

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    Ln_br1 <- branch_likelihood_rate(br = Y[[(brs[[index[1],2]])]], dt = brs$V3[index[1]], Q = Q, delta = m[index[1]])
    Ln_br2 <- branch_likelihood_rate(br = Y[[(brs[[index[2],2]])]], dt = brs$V3[index[2]], Q = Q, delta = m[index[2]])

    Y[[t]] <- node_likelihood(Ln_br1 = Ln_br1, Ln_br2 = Ln_br2)

    if (nb %in% scaling_time){ #every 50 nodes we will perform a scaling event

      Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
      Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states

    }

  }

  t <- unique(brs$V1)[length(unique(brs$V1))]

  LnL <- log(sum(Y[[t]]*(1/length(Y[[t]])))) + sum(log(Lm)) + sum(log(dens_gamma(m,alpha,beta)))

  return(list (Y = Y, Lm = Lm , LnL = LnL, scaling_time = scaling_time, Delta = m))

}



#' Likelihood function for the whole phylogenetic tree with branch rate multiplier
#'
#' Function to only compute the likelihood of a phylogenetic tree with branch rate multiplier
#'
#' @param Y Probability of presence in each cell across the tree for each species (list of arrays)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Q transition matrix between the states (matrix)
#' @param m vector of branch rates multiplier organized in the same order as the branches in the brs file (array)
#' @param alpha alpha parameter of the gamma distribution (float)
#' @param beta beta parameter of the gamma distribution (float)
#'
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list)
#'

tree_likelihood_only <- function(Y,brs,Q, m, alpha, beta){

  LnL <- 0
  if (length(brs$V1)>10){
    scaling_time <- seq(10,length(brs$V1),10) #We will save the timing when to perform likelihood scaling to avoid underflow
  } else {
    scaling_time <- 10
  }
  Lm <- rep(1, max(brs$V1)) #Here we store the scaling parameters across the tree

  nb <- 0

  for (t in (unique(brs$V1))){ #We will go through each nodes
    nb <- nb+1 #We check the number of iteration

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    Ln_br1 <- branch_likelihood_rate(br = Y[[(brs[[index[1],2]])]], dt = brs$V3[index[1]], Q = Q, delta = m[index[1]])
    Ln_br2 <- branch_likelihood_rate(br = Y[[(brs[[index[2],2]])]], dt = brs$V3[index[2]], Q = Q, delta = m[index[2]])

    Y[[t]] <- node_likelihood(Ln_br1 = Ln_br1, Ln_br2 = Ln_br2)

    if (nb %in% scaling_time){ #every 10 nodes we will perform a scaling event

      Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
      Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states

    }

    if (log(prod(Y[[t]]*(1/length(Y[[t]])))) == -Inf){ #Additional scaling if we rich scaling issue
      Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
      Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states
      scaling_time <- sort(c(scaling_time, nb))
    }

    #Lik <- sum(log(Y[[t]])) + log(Lm[t]) #We calculate here the likelihood of the branch

  }

  t <- unique(brs$V1)[length(unique(brs$V1))]

  LnL <- log(prod(Y[[t]]*(1/length(Y[[t]])))) + sum(log(Lm)) + sum(log(dens_gamma(m, alpha, beta)))

  if (LnL==-Inf){
    LnL <- NA
  }
  return(as.numeric(-LnL))

}


#' Likelihood function for the whole phylogenetic tree without branch rate multiplier
#'
#' Function to only compute the likelihood of a phylogenetic tree with branch rate multiplier
#'
#' @param Y Probability of presence in each cell across the tree for each species (list of arrays)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Q transition matrix between the states (matrix)
#' @param m vector of branch rates multiplier organized in the same order as the branches in the brs file (array)
#' @param alpha alpha parameter of the gamma distribution (float)
#' @param beta beta parameter of the gamma distribution (float)
#'
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list)
#'

tree_likelihood_only_no_rates <- function(Y,brs,Q){

  LnL <- 0
  if (length(brs$V1)>50){
    scaling_time <- seq(50,length(brs$V1),50) #We will save the timing when to perform likelihood scaling to avoid underflow
  } else {
    scaling_time <- 50
  }
  Lm <- rep(1, max(brs$V1)) #Here we store the scaling parameters across the tree

  nb <- 0

  for (t in (unique(brs$V1))){ #We will go through each nodes
    nb <- nb+1 #We check the number of iteration

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    Ln_br1 <- branch_likelihood_rate(br = Y[[(brs[[index[1],2]])]], dt = brs$V3[index[1]], Q = Q, delta = 1)
    Ln_br2 <- branch_likelihood_rate(br = Y[[(brs[[index[2],2]])]], dt = brs$V3[index[2]], Q = Q, delta = 1)

    Y[[t]] <- node_likelihood(Ln_br1 = Ln_br1, Ln_br2 = Ln_br2)

    if (nb %in% scaling_time){ #every 50 nodes we will perform a scaling event

      Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
      Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states

    }

    #Lik <- sum(log(Y[[t]])) + log(Lm[t]) #We calculate here the likelihood of the branch

  }

  t <- unique(brs$V1)[length(unique(brs$V1))]

  LnL <- log(sum(Y[[t]]*(1/length(Y[[t]])))) + sum(log(Lm))

  return(LnL)

}








#' Likelihood function for the whole phylogenetic tree
#'
#' Function to compute the likelihood of a phylogenetic tree only
#'
#' @param Y Probability of presence in each cell across the tree for each species (list of arrays)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Q transition matrix between the states (matrix)
#'
#' @return the likelihood at the node as a list of the likelihood of each given states (list)
#'

# tree_likelihood_only <- function(Y,brs,Q){
#
#   LnL <- 0
#   if (length(brs$V1)>50){
#     scaling_time <- seq(50,length(brs$V1),50) #We will save the timing when to perform likelihood scaling to avoid underflow
#   } else {
#     scaling_time <- 50
#   }
#   Lm <- rep(1, max(brs$V1)) #Here we store the scaling parameters across the tree
#
#   nb <- 0
#
#   for (t in (unique(brs$V1))){ #We will go through each nodes
#     nb <- nb+1 #We check the number of iteration
#
#     #Select the two species' distributions we are working with
#     index <- grep(TRUE, brs$V1==t)
#
#     Ln_br1 <- branch_likelihood(br = Y[[(brs[[index[1],2]])]], dt = brs$V3[index[1]], Q = Q)
#     Ln_br2 <- branch_likelihood(br = Y[[(brs[[index[2],2]])]], dt = brs$V3[index[2]], Q = Q)
#
#     Y[[t]] <- node_likelihood(Ln_br1 = Ln_br1, Ln_br2 = Ln_br2)
#
#     if (nb %in% scaling_time){ #every 50 nodes we will perform a scaling event
#
#       Lm[t] <- max(Y[[t]]) #here we store the value of the scaling parameter
#       Y[[t]] <- Y[[t]]/Lm[t] #here we scale across the different states
#
#     }
#
#     #Lik <- sum(log(Y[[t]])) + log(Lm[t]) #We calculate here the likelihood of the branch
#
#   }
#
#   t <- unique(brs$V1)[length(unique(brs$V1))]
#
#   LnL <- log(prod(Y[[t]]*(1/length(Y[[t]])))) + sum(log(Lm))
#
#   return(LnL)
#
# }



#'Description:
#'exp_pruning_simple_split(brs, Y, R, Q)
#'
#'Function to compute the likelihood of ancestral ranges of species using matrix exponentiation
#'We are going backward in time to compute the likelihood and saving the centroid of each likelihood
#'distribution
#'
#'Input:
#'- brs: table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#'- neighbor: list of the neighbor of each cell (list of list)
#'- Y: distribution of the species across the map with for each species the probability to be present and the probability to be absent(list of arrays)
#'- infQ: transition matrix (matrix)
#'-width: width of the map (float)
#'-height: of the map (float)
#'true_id: list of the occupied cells id (list)
#'


exp_pruning_simple_split <- function(n, brs, Y, neighbor, infQ, width, height, true_id){
  #In this function we will go backward in the tree and compute the likelihood of
  #species distribution as well as the centroid of all the distribution on the map
  #in order to use it in the next step for range split at speciation

  center <- rep(0, 2*length(Y))
  for (t in (unique(brs$V1))){
    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    col_start <- seq(1,width*height,height)
    row_start <- c(1:height)

    #For the first branch
    br1 <- Y[[(brs[[index[1],2]])]]
    laybis <- matrix(0, height, width)
    laybis[true_id] <- br1


    centers <- center_of_mass(laybis, width, height)
    center[(index[1]*2)-1] <- centers[1]
    center[(index[1]*2)] <- centers[2]

    #For the second branch
    br2 <- Y[[(brs[[index[2],2]])]]
    laybis <- matrix(0, height, width)
    laybis[true_id] <- br2

    centers <- center_of_mass(laybis, width, height)
    center[(index[2]*2)-1] <- centers[1]
    center[(index[2]*2)] <- centers[2]

    dt1 <- brs$V3[index[1]]
    dt2 <- brs$V3[index[2]]

    Y[[t]] <- br1%*%t(expm(infQ*dt1)) * br2%*%t(expm(infQ*dt2))

    if (t == n+1){
      laybis <- matrix(0, height, width)
      laybis[true_id] <- Y[[t]]

      centers <- center_of_mass(laybis, width, height)
      center[(t*2)-1] <- centers[1]
      center[(t*2)] <- centers[2]

    }

  }

  return(Y, center)

}

