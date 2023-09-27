'
Script developed by: Thomas Merrien
last update: 24/08/2023
File name: Stochastic_mapping.R

Description: functions for performing stochastic mapping of characters.
'

#' what is needed?
#' Estimate the extirpation parameter otherwise it is not gonna work
#' Draw the number of events from a distribution
#' Decide what is the starting point

stoch_map <- function(nsim, n, brs, neighbor, Y, Q, Delta, init_proba, in_state){

  #creation of object that will store the multiple phylogenies --> this part was taken from Phytools developped by Liam Revell
  multi_trees <- vector(mode="list", length=nsim)
  Y_list <- rep(Y, nsim)


  #Creation of the appropriate neighbouring data (without 0)
  neighbor2 <- neighbor

  for (i in 1:length(neighbor2)){

    neighbor2[[i]] <- neighbor2[[i]][neighbor2[[i]]!=0]

  }

  for (i in 1:nsim){

    #We will go forward in the tree using likelihood of presence/absence

    t_ini <- rev(unique(brs$V1))[1]
    t <- t_ini

    #Step 0: Initialization of the tree probabilities
    if (init_proba == "equiproba"){

      #set the probability of each cell to the same amount (equiprobability scenario)
      Pi[[t]] <- 1/length(Pi[[t]])

      #Draw the initial state
      #3 possibilities for that: either you specify the number of cells at the beginning
      #Or they will be drawn from a random number generator:
      #which will be gamma distributed if we suspect small areas of origin
      #and uniformly distributed if we have no idea

      if (is.numeric(in_state) == TRUE){

        in_state <- as.numeric(round(in_state))

        fr
        Pi[[t]] * Y[[t]] / sum(Pi[[t]] * Y[[t]])

      } else {

        if (in_state == "small"){


        } else {

          if (in_state == "unknown"){


          } else {

            return('Set a correct value for the "in_state" variable. "in_state" can take values that are either integer, or string that can be "small" if the original area is small, and "unknown" if you have no knowledge of the potential size of the original area occupied by the ancestral species.')

          }

        }

      }

      P[[t]] <- Pi[[t]] * Y[[t]] / sum(Pi[[t]] * Y[[t]])

    }

      #Simulation continues

  }


}
