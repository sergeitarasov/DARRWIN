'
Script developed by: Thomas Merrien
last update: 07/02/2023
File name: list_neighbor.R
'

#' Finding neighbors cells in matrix
#'
#' Function that for each pixel of the grid list its neighbouring cells
#'
#' @param map matrix filled with the number of each cells (matrix)
#'
#' @return a list of list which returns for each cell the list of the neighbouring cells.
#' each list ends with 2 zeros to have a minimum length of 2


list.neighbor <- function(map){

  neighbor <- list()

  r <- length(map)
  h <- dim(map)[1]
  l <- dim(map)[2]
  row_start <- seq(1, l)
  row_end <- seq(h*(l-1)+1, r)
  col_start <- seq(1,r,by=h)
  col_end <- seq(h,r,by=h)


  for (i in 1:r){
    if (length(grep(TRUE,row_start==i))==0 &
        length(grep(TRUE,row_end==i))==0 &
        length(grep(TRUE,col_start==i))==0 &
        length(grep(TRUE,col_end==i))==0){
      nb <- c(i-1, i+1, i-h, i+h)
      nb[nb<=0 | nb>r] <- NA
      nb <- nb[!is.na(nb)]
      neighbor[[i]] <- c(nb,0,0)
    } else {
      if (length(grep(TRUE,row_end==i))>0){
        nb <- c(i-1, i+1, i-h)
        nb[nb<=0 | nb>r] <- NA
        if (length(grep(TRUE,col_start==i))>0){
          nb[nb==(i-1)] <- NA
        }
        nb <- nb[!is.na(nb)]
        neighbor[[i]] <- c(nb,0,0)
      } else {
        if (length(grep(TRUE,col_end==i))>0){

          nb <- c(i-1, i-h, i+h)
          nb[nb<=0 | nb>r] <- NA
          nb <- nb[!is.na(nb)]
          neighbor[[i]] <- c(nb,0,0)

        } else {

          if (length(grep(TRUE, row_start==i))>0){

            nb <- c(i-1, i+1, i+h)
            nb[nb<=0 | nb>r] <- NA
            nb <- nb[!is.na(nb)]
            neighbor[[i]] <- c(nb,0,0)

          } else{

              nb <- c(i+1, i-h, i+h)
              nb[nb<=0 | nb>r] <- NA
              nb <- nb[!is.na(nb)]
              neighbor[[i]] <- c(nb,0,0)

          }
        }
      }
    }
  }

  return(neighbor)
}


#' Finding neighbors cells in matrix second version
#'
#' Function that for each pixel of the grid list gives its neighbouring cells according to the pixel edges
#'
#' @param map dataframe with the pixel id and edges
#'
#' @return a list of list which returns for each cell the list of the neighbouring cells.
#'

list.neighbor2 <- function(map){

  neighbor <- matrix(0, nrow = dim(map)[1], ncol = dim(map)[1])

  list_id <- sort(map$id)

  for (i in 1:dim(map)[1]){

    for (j in 1:dim(map)[1]){

      if (length(intersect(c(map[i,9:12]), c(map[j,9:12])))>2){

        neighbor[grep(TRUE, map$id[i] == list_id),grep(TRUE, map$id[j] == list_id)] <- 1

      }

    }

  }

  return(neighbor)

}



