'
Script developed by: Thomas Merrien
last update: 26/04/2022
File name: plot_main.R

Description:
Plotting the results of the main test section
'


# Libraries --------------------------------------------------------------------
library(plot.matrix)
library(plotly)
library(rgl)

# create a directory to which the images will be written
dir_out <- file.path(getwd(), "test/distrib")
dir.create(dir_out, recursive = TRUE)

## prepare data
name <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")

#for (i in 1:max(brs[,1])){
for (i in 1:20){
  gridbis <- Y_true[[(2*i)-1]]

  fp <- file.path(dir_out, paste0("distrib_V1_", name[i], ".png"))
  png(fp)

  p <- plot(gridbis, breaks = 50)

  dev.off()

}

for (i in 1:1){
  gridbis <- Ptest[[((n+1)*2)-1]]
  fp <- file.path(dir_out, paste0("distrib_Proba_", ".png"))
  png(fp)

  p <- plot(gridbis, breaks = 50)

  dev.off()
}

persp(x=c(1:20),y=c(1:20), z=Y_true[[41]])
persp(x=c(1:20),y=c(1:20), z=Ptest[[41]])


persp3d(x=c(1:20),y=c(1:20), z=Y_true[[41]])
persp3d(x=c(1:20),y=c(1:20), z=Ptest[[41]], col="green")

#Plot the data with lemurs -----------------------------------------------------

#Real grid
real_grid <- read.csv("test/full_grid.csv")

width = 10
height = 18

grid = matrix(0,height,width)
forest = matrix(NA,height,width)

true_id <- cell_id
missing_id <- setdiff(c(1:(width*height)), true_id)

for (i in 1:length(Ptest[[1]])){
  if (Ptest[[25]][i]==1){
    grid[true_id[i]] <- 0
  } else {
    grid[true_id[i]] <- Ptest[[35]][i]
  }

  if (R[i] == "land"){
    forest[true_id[i]] <- 0
  } else {
    forest[true_id[i]] <- 1
  }
}


grid[missing_id] <- NA

plot(grid, col = rev(heat.colors(10)))
persp(x=c(1:18),y=c(1:10), z=grid, col = t(forest))

grid = matrix(0,height,width)
for (i in 1:length(Y[[1]])){
  grid[true_id[i]] <- Y[[29]][i]
}
grid[missing_id] <- NA
plot(grid)







#Test plot nice map ------------------------------------------------------------
mat_test <- matrix(0,10,10)
neighbor_test <- list.neighbor(mat_test)
infQ_test <- matrix(0,100,100)
for (i in 1:100){
  for (j in neighbor_test[[i]]){
    infQ_test[i,j] <- 0.25
  }
  infQ_test[i,i] <- -sum(infQ_test[i,])
}

mat_test <- rep(0,100)
mat_test[55] <- 1
mat_test <- mat_test %*% matexpo(infQ_test*5)
mat_test <- matrix(mat_test, nrow = 10, ncol = 10)*10

mat_test2 <- rep(0,100)
mat_test2[95] <- 1
mat_test2 <- mat_test2 %*% matexpo(infQ_test*5)
mat_test2 <- matrix(mat_test2, nrow = 10, ncol = 10)*10


w2 <- (mat_test)/max(mat_test)
cr <- colorRamp(c('white','yellow','red'))
w3 <- cr(w2)
w4 <- rgb(w3[,1], w3[,2], w3[,3], maxColorValue=255)
dim(w4) <- c(10,10)

open3d()
surface3d(c(1:10), c(1:10), z=mat_test, col = w4)
m <- as.mesh3d()
close3d()

col.pal<-colorRampPalette(c('white','yellow','red'))
colors<-col.pal(100)
# height of facets
z.facet.center <- (mat_test[-1, -1] + mat_test[-1, -ncol(mat_test)] + mat_test[-nrow(mat_test), -1] + mat_test[-nrow(mat_test), -ncol(mat_test)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=c(1:10),y=c(1:10), z=mat_test, col = colors[z.facet.range], xlab = NA, ylab = NA, theta = 0, phi = 30, r = 2, d = 3, box = FALSE, scale = FALSE, ltheta = -120)

open3d()
surface3d(c(1:10), c(1:10), z=mat_test, col = w4)
m <- as.mesh3d()
close3d()


col.pal2<-colorRampPalette(c('white','lightblue','blue'))
colors2<-col.pal2(100)
# height of facets
z.facet.center2 <- (mat_test2[-1, -1] + mat_test2[-1, -ncol(mat_test2)] + mat_test2[-nrow(mat_test2), -1] + mat_test2[-nrow(mat_test2), -ncol(mat_test2)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range2<-cut(z.facet.center2, 100)

persp(x=c(1:10),y=c(1:10), z=mat_test2, col = colors2[z.facet.range2], xlab = NA, ylab = NA, theta = 0, phi = 30, r = 2, d = 3, box = FALSE, scale = FALSE, ltheta = -120)

#Trying to mix the two matrix to get one map
mat_test_mix <- (mat_test2+mat_test)/2
z.facet.center_mix <- (mat_test_mix[-1, -1] + mat_test_mix[-1, -ncol(mat_test_mix)] + mat_test_mix[-nrow(mat_test_mix), -1] + mat_test_mix[-nrow(mat_test_mix), -ncol(mat_test_mix)])/4
z.facet.range_mix<-cut(z.facet.center_mix, 100)

persp(x=c(1:10),y=c(1:10), z=mat_test_mix, col = colors[z.facet.range_mix], xlab = NA, ylab = NA, theta = 0, phi = 30, r = 2, d = 3, box = FALSE, scale = FALSE, ltheta = -120)

#### INTERACTIVE PLOT ####
ind = 1
fig <- plot_ly(
    z = (Ptest[[ind]]*1),
    frame = ind,
    type = 'heatmap')

fig
