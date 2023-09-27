#'
#'Script developed by: Thomas Merrien
#'last update: 03/07/2023
#'File name: plot_result.R
#'
#'Description:
#'Displaying the results of the ancestral ranges reconstruction.
#'
#'


#### LIBRARIES ####

library(phytools)

#### AUXILIARY FUNCTIONS ####
# asp and add_grid function are coming from Phytools blog and are written by Liam Revell. Source here: http://blog.phytools.org/2023/03/adding-grids-to-nodes-of-plotted-tree.html

# function to compute effective aspect ratio -----------------------------------
asp<-function(){
  dx<-diff(par()$usr[1:2])
  dy<-diff(par()$usr[3:4])
  asp<-(dy/dx)*(par()$pin[1]/par()$pin[2])
  asp
}

# function to add grid to a node -----------------------------------------------
add_grid<-function(x,y,nr,nc,cols,box_size=0.75, missing){
  box_size<-box_size/asp()
  xleft<-x-box_size/2
  ybottom<-y-box_size/2*asp()*(nr/nc)
  xright<-x+box_size/2
  ytop<-y+box_size/2*asp()*(nr/nc)
  nn<-max(c(nr,nc))
  m<-1
  for(k in 0:(nr-1)) for(j in 0:(nc-1)){
    xx<-c(xleft+k*box_size/nn,
          xleft+(k+1)*box_size/nn)
    yy<-c(ytop-j*box_size/nn*asp(),
          ytop-(j+1)*box_size/nn*asp())
    if (m %in% missing){
      rect(xx[1],yy[1],xx[2],yy[2],col="lightblue")
    } else {
      rect(xx[1],yy[1],xx[2],yy[2],col=cols[m])
    }

    m<-m+1
  }
}

#### MAIN FUNCTION ####

#' Plot phylogenetic tree with the ancestral range reconstruction at the node
#'
#' This function allow to display the phylogenetic tree with the range reconstruction of the ancestral species displayed on each nodes
#'
#' @param tree phylogenetic tree
#' @param P array of array of ancestral range probabilities at each nodes for each cell grid of your map
#' @param h The height of the map in number of cells
#' @param w The width of the map in number of cells
#'
#'
#' @return the likelihood of one branch as a list of the likelihood of each given states (list)
#'

plot_phylomap <- function(grid_map, cell_id, tree, P){

  plotTree(tree,ftype="i")
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)


  height = max(length(unique(grid_map$left)),length(unique(grid_map$right)))+1
  width = max(length(unique(grid_map$top)),length(unique(grid_map$bottom)))+1

  missing_id <- setdiff(c(1:(width*height)), cell_id)

  ## set colors
  #cols<-setNames(colorRampPalette(c("white","darkgreen"))(6),4:9)
  cols<-colorRampPalette(c("white","darkred"))(1000)

  ## pull labels
  labs<-c(tree$tip.label,1:tree$Nnode+Ntip(tree))

  ## add grids
  for(i in 1:(length(tree$tip.label) + tree$Nnode)){
    grid <- matrix(0, height, width)
    x<-pp$xx[i]
    y<-pp$yy[i]

    for (j in 1:length(P[[i]])){
      if (P[[i]][j]>0){
        grid[cell_id[j]] <- P[[i]][j]
      } else {
        grid[cell_id[j]] <- 0
      }
    }

    add_grid(x,y,height,width,cols[grid*999 + 1], missing = missing_id)
  }

  ## add legend
  #legend("bottomleft",names(cols),pch=22,pt.bg=cols,cex=1,
  #       pt.cex=2,bty="n")


}
