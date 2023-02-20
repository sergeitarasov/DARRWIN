'
Script developed by: Thomas Merrien
last update: 26/04/2022
File name: pruning_algo.R
'
#'Description:
#'exp_pruning_simple(brs, Y, R, Q)
#'
#'Function to compute the likelihood of ancestral ranges of species using matrix exponentiation
#'We are going backward in time to compute the likelihood
#'
#'Input:
#'- brs: table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#'- neighbor: list of the neighbor of each cell (list of list)
#'- Y: distribution of the species across the map with for each species the probability to be present and the probability to be absent(list of arrays)
#'- R: matrix of envrionment type that follow the grided map (matrix)
#'- Q: transition matrix (matrix)
#'- dt: small delta time used for the Euler approximation of the likelihood calculation (float)
#'

matrix_pruning_simple <- function(brs, Y, neighbor, infQ, dt){

  neighbor2 <- neighbor

  for (i in 1:length(neighbor2)){

    neighbor2[[i]] <- neighbor2[[i]][neighbor2[[i]]!=0]

  }
  #We will go backward in the tree and compute the potential distribution of
  #ancestors using matrix exponentiation
  for (t in (unique(brs$V1))){

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    lay1s1 <- Y[[(brs[[index[1],2]])]]
    lay2s1 <- Y[[(brs[[index[2],2]])]]


    #for the first branch

    rep <- reps_per_period(brs[[index[1],3]], dt)
    for (i in 1:rep){

      lay1st1 <- matrix(0, dim(lay1s1)[1], dim(lay1s1)[2])
      lay2st1 <- matrix(0, dim(lay2s1)[1], dim(lay2s1)[2])

      #for the first state: the cell is occupied

      occupancy <- grep(TRUE, lay1s1!=0)
      occupancy <- unique(c(occupancy, sapply(neighbor[occupancy],"[[",1),sapply(neighbor[occupancy],"[[",2),sapply(neighbor[occupancy],"[[",3),sapply(neighbor[occupancy],"[[",4)))
      occupancy <- occupancy[occupancy!=0]

      for (cell in occupancy){

        lay1st1[cell] <- lay1st1[cell] + lay1s1[cell]*infQ[cell,cell]

        for (nb in neighbor[[cell]][neighbor[[cell]]>0]){

          lay1st1[cell] <- min(1,lay1st1[cell] + (1-lay1s1[cell]) * lay1s1[nb]*(infQ[nb,cell]*(1/length(neighbor2[[cell]]))))

        }

      }

      lay1s1 <- lay1st1

    }

    #for the second branch

    rep <- reps_per_period(brs[[index[2],3]], dt)
    for (i in 1:length(rep)){

      occupancy <- grep(TRUE, lay2s1!=0)
      occupancy <- unique(c(occupancy, sapply(neighbor[occupancy],"[[",1),sapply(neighbor[occupancy],"[[",2),sapply(neighbor[occupancy],"[[",3),sapply(neighbor[occupancy],"[[",4)))
      occupancy <- occupancy[occupancy!=0]

      for (cell in occupancy){

        lay2st1[cell] <- lay2st1[cell] + lay2s1[cell]*infQ[cell,cell]

        for (nb in neighbor[[cell]][neighbor[[cell]]>0]){

          lay2st1[cell] <- min(1,lay2st1[cell] + lay2s1[nb] * (1-lay2s1[cell]) * (infQ[nb,cell]*(1/length(neighbor2[[cell]]))))

        }

      }

      lay2s1 <- lay2st1

    }

    Y[[t]] <- lay1s1 * lay2s1

  }

  return(Y)

}

#' Likelihood calculation across the tree
#'
#' Function to compute the likelihood of ancestral ranges of species using matrix exponentiation
#' We are going backward in time to compute the likelihood.
#'
#' @param brs table summarizing the phylogenetic tree with the branch length and the nodes that are linked(array)
#' @param Y distribution of the species across the tree with for each species the probability to be present and the probability to be absent(list of arrays)
#' @param neighbor list of the neighbor of each cell (list of list)
#' @param infQ transition matrix (matrix)
#' @param dt
#'
#' @return the likelihood at each nodes
#'


exp_pruning_simple <- function(brs, Y, neighbor, infQ){

  #We will go backward in the tree and compute the potential distribution of
  #ancestors using matrix exponentiation
  for (t in (unique(brs$V1))){

    #Select the two species' distributions we are working with
    index <- grep(TRUE, brs$V1==t)

    br1 <- Y[[(brs[[index[1],2]])]]
    br2 <- Y[[(brs[[index[2],2]])]]

    dt1 <- brs$V3[index[1]]
    dt2 <- brs$V3[index[2]]


    #for the first state
    Y[[t]] <- br1%*%t(expm(infQ*dt1)) * br2%*%t(expm(infQ*dt2))

  }

  return(Y)

}




#' Likelihood function for trees
#'
#' Function to compute the overall likelihood of the tree
#'
#' @param Y distribution of the likelihood of each states across the tree (list of arrays)
#'
#'
#' @return the overall likelihood of the tree
#'

tree_likelihood <- function(Y){

  Ln <- 0

  return (Ln)

}



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

