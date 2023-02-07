#'
#'Script developed by: Thomas Merrien
#'last update: 26/04/2022
#'File name: data_generation.R
#'
#'Description:
#'Generating the different data needed to test the pruning algorithm of our
#'ancestral range reconstruction model
#'
#'

# Grid map parameters ----------------------------------------------------------

w = 3 #width: number of columns
h = 1 #height: number of rows
landtype = c("land", "water")
landtype_proportion <- c(0.9, 0.1)

# Generation of the grid map landscape structure -------------------------------

#' Random map generation
#'
#' This function will create a rectangular grid of dimension h x w that will be
#' filled with different landscape types according to the proportion given by the
#' user
#'
#' @param h The height of the map in number of cells
#' @param w The width of the map in number of cells
#' @param landtype The list of landscape type you want in your map
#' @param prop The proportion of each landscape type in your map. Proportion are given in the same order as the landtype and should sum to one
#'
#' @return Returns a matrix of dimension h x w filled with the different landscape types
#' @export
#'
#' @examples
#' #Generating map with land and water of dimension 1 x 3 with 1/3 of water
#' landscape_generation(1,3,c("land","water"), c(2/3,1/3))

landscape_generation <- function(h,w,landtype, prop){

  R <- matrix(landtype[1],nrow = h, ncol = w)
  if (length(landtype)>1){
    for (i in 2:length(landtype)){
      R[sample(R[R==landtype[1]],round(prop[i]*w*h))] <- landtype[i]
    }
  }

  return (R)
}

# Generation of the indexed map ------------------------------------------------
map <- matrix(1:(h*w), ncol = w, nrow = h)

# Generation of the transition matrix ------------------------------------------

#Area types
geo_types <- c("land", "water")

K = length(geo_types)
Q <- matrix(0, nrow = K, ncol = K)
colnames(Q) <- geo_types
rownames(Q) <- geo_types

Q[1]=0.10
Q[2]=0.10
Q[3]=0.10
Q[4]=0.10

# Delta time for the analysis --------------------------------------------------

dt <- 0.05

# Generation of the square occupation probability matrices ---------------------

Pi0 <- list(matrix(0.5, nrow = h, ncol = w),matrix(0, nrow = h, ncol = w))
Pi1 <- list(matrix(0.5, nrow = h, ncol = w),matrix(0, nrow = h, ncol = w))

# Generation of a random phylogenetic tree -------------------------------------

n <- 3 #number of species

library(ape)
set.seed(78)
tree <- rcoal(n,
              rooted = TRUE)

# Create the list of neighbor --------------------------------------------------
neighbor <- list.neighbor(map) #Function neighbor need to be checked

# Get the phylogenetic tree structure ------------------------------------------
brs <- branching_time(tree)
brs <- data.table::as.data.table(brs)
colnames(brs) <- c("V1", "V2", "V3", "V4", "V5")
brs <- dplyr::arrange(brs, desc(V4))
#increase resolution (need to lower the R parameter if done)
brs[,3:5] <- brs[,3:5]


# Adjacency and Q matrix design ------------------------------------------------

#Adjacency matrix receive 1 if cells are neighbor, 0 otherwise
A <- matrix(0, nrow=h*w, ncol = h*w)
#Rate matrix at the grid level
infQ <- matrix(0, nrow=h*w, ncol = h*w)
matQ <- matrix(0, nrow=h*w, ncol = h*w)
infQdt <- matrix(0, nrow=h*w, ncol = h*w)

for(i in 1:(h*w)){
  for (j in neighbor[[i]]){
    A[i,j] <- 1
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

# Generation of the biogeographic dataset --------------------------------------

Y <- list(matrix(0, nrow = h, ncol = w),matrix(0, nrow = h, ncol = w))
Y <- rep(Y,tree$Nnode+n)

#creation of simulated biogeographic history
Y_true <- Y

Y_true[[(n+1)*2-1]][c(361,362,363,364,381,382,383,384)] <- 1
Y_true[[(n+1)*2]][] <- 1
Y_true[[(n+1)*2]][c(361,362,363,364,381,382,383,384)] <- 0

# random version ---------------------------------------------------------------

for (t in (rev(unique(brs$V2)))){

  rep <- reps_per_period(brs$V3[grep(TRUE, brs$V2==t)], dt)
  Y_true[[(2*t)-1]] <- Y_true[[(2*brs$V1[grep(TRUE, brs$V2==t)])-1]]

  for (i in 1:rep){

    for (cell in grep(TRUE,Y_true[[(t*2)-1]]>0)){

      random <- runif(1)
      prob <- cumsum(c(0.05*dt, infQdt[cell,neighbor[[cell]]]))
      prob <- sort(c(random,prob))

      if (grep(TRUE, prob==random)[1]==1){

        Y_true[[(2*t)-1]][cell] <- 0
        print(random)

      } else {

        if (grep(TRUE, prob==random)<(length(neighbor[[cell]])+2)) {

          Y_true[[(2*t)-1]][neighbor[[cell]][grep(TRUE, prob==random)[1]-1]] <- 1

        }
      }

    }

  }

}

for (i in 1:n){
  Y[[(2*i)-1]] <- Y_true[[(2*i)-1]]
  Y[[(2*i)]] <- 1-Y_true[[(2*i)-1]]
}

Y_v2 <- list(rep(0,h*w),rep(0,h*w))
Y_v2 <- rep(Y_v2,tree$Nnode+n)

for (i in 1:n){
  for (j in 1:h*w){
    Y_v2[[(2*i)-1]][j] <- Y[[(2*i)]][j]
    Y_v2[[(2*i)]][j] <- Y[[(2*i)]][j]
  }
}

# deterministic version --------------------------------------------------------
for (t in (rev(unique(brs$V2)))){

  rep <- reps_per_period(brs$V3[grep(TRUE, brs$V2==t)], dt)
  Y_true[[(2*t)-1]] <- Y_true[[(2*brs$V1[grep(TRUE, brs$V2==t)])-1]]

  for (i in 1:rep){

    for (cell in grep(TRUE,Y_true[[(t*2)-1]]>0)){

      for (nb in neighbor[[cell]][neighbor[[cell]]>0]){

        Y_true[[(2*t)-1]][nb] <- min(1, Y_true[[(2*t)-1]][nb]+(1-Y_true[[(2*t)-1]][nb])*Y_true[[(2*t)-1]][cell]*infQdt[cell,nb])

      }

      Y_true[[(2*t)-1]][cell] <- Y_true[[(2*t)-1]][cell]*infQdt[cell, cell]

    }

  }

}


# V1 ---------------------------------------------------------------------------
## For species 1 ##
number = 1
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 2 ##
number = 2
Y[[(number*2)-1]][c(81,82)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 3 ##
number = 3
Y[[(number*2)-1]][c(41,42,43,63, 44, 64, 45, 65, 61, 62, 82, 83, 84)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 4 ##
number = 4
Y[[(number*2)-1]][c(81,82,83,101,102,103,104)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 5 ##
number = 5
Y[[(number*2)-1]][c(63,64,83,84,103,104)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]


# V2 ---------------------------------------------------------------------------

## For species 1 ##
number = 1
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 2 ##
number = 2
Y[[(number*2)-1]][c(81,82)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 3 ##
number = 3
Y[[(number*2)-1]][c(41,42,43,63, 44, 64, 45, 65, 61, 62, 82, 83, 84)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 4 ##
number = 4
Y[[(number*2)-1]][c(88,108, 89, 109, 110)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 5 ##
number = 5
Y[[(number*2)-1]][c(67,87,88,108)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 6 ##
number = 6
Y[[(number*2)-1]][c(154,174,155,195,156,176,196)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 7 ##
number = 7
Y[[(number*2)-1]][c(135, 155, 156)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 8 ##
number = 8
Y[[(number*2)-1]][c(52,53,54,55,56,72,73,74,75,76,92,93,94,95,96)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 9 ##
number = 9
Y[[(number*2)-1]][c(14,16,34,35,36)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 10 ##
number = 10
Y[[(number*2)-1]][c(33,53,34,54)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]


# V3 ---------------------------------------------------------------------------
## For species 1 ##
number = 1
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 2 ##
number = 2
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 3 ##
number = 3
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 4 ##
number = 4
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]

## For species 5 ##
number = 5
Y[[(number*2)-1]][c(1,2,3,21,23,41,43,61,62,63)] <- 1
Y[[(number*2)]] <- 1-Y[[(number*2)-1]]


