'
Script developed by: Thomas Merrien
last update: 23/11/2022
File name: center_of_mass.R

Description:
center_of_mass(Y, width, height)

Function to compute the center of mass of each cells

Input:
- Y: distribution of the species across the map with for each species the
probability to be present (list of arrays)
- width: the width of the initial rectangular grid map (int)
- height: the height of the initial rectangular grid map (int)
'

center_of_mass <- function(Y, width, height){

  center <- c(0,0)
  for (j in 1:height){
    center[1] <- center[1] + (sum(c(0.5:(width-0.5))*Y[seq(j,width*height,height)])/sum(Y[seq(j,width*height,height)]))
  }
  center[1] <- center[1]/height

  for (j in seq(1,width*height,height)){
    center[2] <- center[2] + (sum(c(0.5:(height-0.5))*Y[j:(height+(j-1))])/sum(Y[j:(height+(j-1))]))
  }
    center[2] <- center[2]/width

  return(center)

}

