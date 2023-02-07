'
Script developed by: Ignacio Quintero adapted to R by Thomas Merrien
last update: 22/04/2022
File name: branching_time.R
'

#'Description:
#'branching_time(tree)
#'
#'Function adapted from package TRIBE from Ignacio Quintero that estimate absolute branching times, with 0 at the present, time going backward
#'
#'Input:
#'-phylogenetic tree (list)



# Branching_time ---------------------------------------------------------------

branching_time  <- function(tree){
  n <- tree$Nnode + 1
  el_t <- grep(TRUE,tree$edge[,2]<=n)

  brs <- matrix(0, length(tree$edge.length),5)

  brs[,1:2] <- tree$edge
  brs[,3] <- tree$edge.length

  for (j in 1:length(tree$edge.length)){
    if (brs[j,1]==(n+1)){
      brs[j,4] <- 0.0
      brs[j,5] <- brs[j,4]-brs[j,3]
    } else {
      brs[j,4] <- brs[,5][grep(TRUE, brs[j,1]==brs[,2])][1]
      brs[j,5] <- brs[j,4]-brs[j,3]
    }
  }

  tree_depth <- brs[grep(TRUE, n==brs[,2]),5][1]

  for (j in 1:length(tree$edge.length)){
    brs[j,4] <- tree_depth - brs[j,4]
    brs[j,5] <- tree_depth - brs[j,5]
  }

  brs[el_t, 5] <- 0.0

  brs <- brs[order(brs[,5]),]

  return(brs)
}
