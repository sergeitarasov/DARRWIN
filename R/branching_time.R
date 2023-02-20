'
Script developed by: Ignacio Quintero adapted to R by Thomas Merrien
last update: 07/02/2023
File name: branching_time.R
'

#' Tree branching times summary
#'
#' Function adapted from package TRIBE from Ignacio Quintero that summarize the
#' tree in a matrix and estimate absolute branching times, with 0 at the present.
#'
#' @param tree a phylogenetic tree (tree)
#'
#' @return Returns a matrix with the first column of the table gives the internal
#' node where the branch start, the second column gives the node or tips where
#' the branch ends. The third column gives the length of the branch. The fourth
#' column gives the length to the root from the starting node (actually it is the
#'  negative length, as it is in fact the date of the branch internal node with
#'  0 being the present).Finally the fifth column gives the date of the external
#'  node/tip, 0 being present and the root being -1.(tree_length).
#'
#' @export
#'
#' @examples
#' #For a simulated tree with 10 tips
#' tree <- ape::rcoal(10)
#' branching_time(tree)



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
