'
Script developed by: Thomas Merrien
last update: 10/05/2022
File name: ancestral_state_probability.R
'

#' Ancestral state probability calculation across tree using matrix exponential
#'
#' Function to compute the probability of ancestral distribution of species using matrix exponentiation and felsenstein pruning algorithm (re-rooting tree)
#' We are going forward in time to compute the probability using the likelihood of each state.
#' User have to choose what are the initial probabilities of presence to be attributed to each cell.
#'
#' @param n number of species (float)
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array). This can be obtained with the branching_time function
#' @param neighbor list of the neighbor of each cells (list of list)
#' @param Y distribution of the species across the map and the time with for each species the likelihood to be present and to be absent (list of arrays). This is obtained through the different pruning functions
#' @param Q transition matrix (matrix)
#' @param Delta vector of branch rates multiplier organized in the same order as the branches in the brs file (array)
#' @param init_proba original a priori knowledge of the probability of cells occupancy by the ancestral species (string). This can be "equiproba" if you have no a priori knowledge of the distribution and that will give the same probability of occupancy to each cell
#'
#'
#' @return Returns a list of matrix with the probability of cell occupancy at each tip and node
#' @export
#'
#' @examples


ancestral_state_probability <- function(n, brs, neighbor, Y, Q, Delta, init_proba){

  neighbor2 <- neighbor

  for (i in 1:length(neighbor2)){

    neighbor2[[i]] <- neighbor2[[i]][neighbor2[[i]]!=0]

  }

  #Create intermediate Q matrix
  Qi <- list(rep(0, length(Y[[1]])))
  Qi <- rep(Qi,tree$Nnode+n)

  #Create matrices of probability at the nodes fill with 0 at the beginning
  Pi <- Qi

  #Create matrices of marginal probability
  P <- Pi

  #We will go forward in the tree using likelihood of presence/absence

  t_ini <- rev(unique(brs$V1))[1]
  t <- t_ini

  #Step 0: Initialization of the tree probabilities
  if (init_proba == "equiproba"){

    Pi[[t]] <- 1/length(Pi[[t]])

    P[[t]] <- Pi[[t]] * Y[[t]] / sum(Pi[[t]] * Y[[t]])

  }

  for (t in brs$V2[brs$V1 %in% brs$V1[brs$V2%in%(rev(unique(brs$V2[brs$V2>n])))]]){

    #Step 1: Creation of the transition matrix

    Qi[[t]] <- Y[[t]]%*%t(expm(Q*(brs$V3[grep(TRUE, brs$V2==t)] * Delta[grep(TRUE, brs$V2==t)])))

  }

  #Step 2: Computing the prior probability at each node

#  for (t in rev(brs$V2[brs$V1 %in% brs$V1[brs$V2%in%((unique(brs$V2[brs$V2>n])))]])){
  for (t in rev(brs$V2[brs$V2>n])){
    anc <- brs$V1[grep(TRUE, brs$V2==t)]
    sis_node <- brs$V2[brs$V1==anc & brs$V2!=t]

    #initialization
    Pi[[t]] <- Qi[[sis_node]] * Pi[[anc]]

    Pi[[t]] <- Pi[[t]]%*%t(expm(Q*(brs$V3[grep(TRUE, brs$V2==t)] * Delta[grep(TRUE, brs$V2==t)])))

    #Step 3: Computing the marginal probability at each node


    P[[t]] <- Pi[[t]]*Y[[t]]/sum(Pi[[t]]*Y[[t]])


  }

  for (j in 1:n){
    P[[j]] <- Y[[j]]
  }

  return(P)

}


ancestral_state_probability_exp_simple_split <- function(n, brs, neighbor, Y, infQ, dt, init_proba, width, height, true_id, center){

  neighbor2 <- neighbor
  x_graph_map <- t(matrix(1:width, width, height))
  y_graph_map <- matrix(1:height, height, width)

  for (i in 1:length(neighbor2)){

    neighbor2[[i]] <- neighbor2[[i]][neighbor2[[i]]!=0]

  }

  #Create intermediate Q matrix
  Qi <- list(matrix(0, nrow = dim(Y[[1]])[1], ncol = dim(Y[[1]])[2]))
  Qi <- rep(Qi,tree$Nnode+n)

  #Create matrices of probability at the nodes fill with 0 at the beginning
  Pi <- Qi

  #Create matrices of marginal probability
  P <- Pi

  #We will go forward in the tree using likelihood of presence/absence

  t_ini <- rev(unique(brs$V1))[1]
  t <- t_ini

  #Step 0: Initialization of the tree probabilities
  if (init_proba == "equiproba"){

    Pi[[t]] <- 1/length(Pi[[t]])

    P[[t]] <- Pi[[t]] * Y[[t]] / sum(Pi[[t]] * Y[[t]])

  }

  for (t in brs$V2[brs$V1 %in% brs$V1[brs$V2%in%(rev(unique(brs$V2[brs$V2>n])))]]){

    #Step 1: Creation the transition matrix

    Qi[[t]] <- Y[[t]]%*%expm(infQ*(brs$V3[grep(TRUE, brs$V2==t)]))

  }

  #Step 2: Computing the prior probability at each node

  for (t in rev(brs$V2[brs$V1 %in% brs$V1[brs$V2%in%((unique(brs$V2[brs$V2>n])))]])){

    anc <- brs$V1[grep(TRUE, brs$V2==t)]
    sis_node <- brs$V2[brs$V1==anc & brs$V2!=t]

    #Setting the range split
    #For the range split we consider that the ancestor species will occupy only
    #one half of its range and that will be what it will transmit to the daughter

    x_anc <- center[(2*anc)-1]
    y_anc <- center[2*anc]

    x_daughters <- center[(2*sis_node)-1]-center[(2*t)-1]
    y_daughters <- center[(2*sis_node)]-center[(2*t)]

    if (x_daughters!=0 & y_daughters!=0){
      y_star <- 0
      x_star <- (x_daughters*y_star + (y_anc*x_daughters-x_anc*y_daughters))/(y_daughters)
    } else {
      if (x_daughters==0){
        x_star <- x_anc
        y_star <- 0
      } else {
        y_star <- y_anc
        x_star <- 0
      }
    }

    #We now can have the split direction which is the altitude of the triangle defined
    #by the 3 center of mass of the ancestor and daughter species likelihoods

    x_alt <- x_star - x_anc
    y_alt <- y_star - y_anc

    #The altitude can be represented as a linear regression that will split the
    #map according to the following expression y = ax + b

    Pi_bis <- Pi[[anc]]

    if (x_alt == 0){
      a <- 1
      b <- 0

      if (center[(2*t)-1] <= x_anc) {

        for (i in 1:length(Pi_bis)){

          if (x_graph_map[true_id[i]] >= x_anc){

            Pi_bis[i] <- 0

          }

        }
      } else {

        for (i in 1:length(Pi_bis)){

          if (x_graph_map[true_id[i]] <= x_anc){

            Pi_bis[i] <- 0

          }

        }

      }

    } else {
      a <- y_alt/x_alt
      b <- y_anc - a*x_anc

      if (center[2*t] <= a*centroid[(2*t)-1]+b){

        for (i in 1:length(Pi_bis)){

          if (y_graph_map[true_id[i]] >= a*graph_map[true_id[i]]+b){

            Pi_bis[i] <- 0

          }

        }

      } else {

        for (i in 1:length(Pi_bis)){

          if (y_graph_map[true_id[i]] <= a*graph_map[true_id[i]]+b){

            Pi_bis[i] <- 0

          }

        }

      }

    }

    #initialization
    Pi[[t]] <- Qi[[sis_node]] * Pi_bis

    Pi[[t]] <- Pi[[t]]%*%expm(infQ*(brs$V3[grep(TRUE, brs$V2==t)]))

    #Step 3: Computing the marginal probability at each node


    P[[t]] <- Pi[[t]]*Y[[t]]/sum(Pi[[t]]*Y[[t]])


  }

  return(P)

}


