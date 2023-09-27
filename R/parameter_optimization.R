'
Script developed by: Thomas Merrien
last update: 27/03/2023
File name: parameter_optimization.R

Description: functions for performing parameter and likelihood optimization
'

######################
#### REQUIREMENTS ####
######################

library(nloptr)

################
#### SCRIPT ####
################

# opts <- list("algorithm" = "NLOPT_GN_DIRECT",
#              "xtol_rel" = 1e-10,
#              maxeval = 1000)
#
# opts2 <- list("algorithm" = "NLOPT_GN_DIRECT_L_RAND",
#              "xtol_rel" = 1e-10,
#              maxeval = 1000)
#
# res <- nloptr(x0 = 0.5,
#               lb = 0,
#               ub = 1,
#               eval_f = optimization_function,
#               opts = opts)
#
# res2 <- nloptr(x0 = c(0.5,0.5),
#               lb = c(0,0),
#               ub = c(1,1),
#               eval_f = optimization_function,
#               opts = opts)

############################
#### AUXILIARY FUNCTION ####
############################

optimization_function_no_rate <- function(lambda){

  Q <- matrix(lambda, length(envt_type), length(envt_type))
  colnames(Q) <- envt_type
  rownames(Q) <- envt_type

  matQ <- matrix(0, nrow=length(neighbor), ncol = length(neighbor))
  for(i in 1:(length(neighbor))){
    for (j in neighbor[[i]]){

      matQ[i,j] <- Q[R[i],R[j]]
    }
  }

  for (i in 1:(length(neighbor))){
    matQ[i,i] <- -sum(matQ[i,])
  }

  Lik<-tree_likelihood_only(Y,brs,matQ, m=rep(1,length(Y)), alpha = 1, beta = 1)

  return(Lik)

}



optimization_function_rate <- function(lambda){

  nb_branch <- length(brs$V1)
  nb_param <- length(lambda)

  Q <- matrix(lambda[1:(nb_param-nb_branch-2)], length(envt_type), length(envt_type))
  colnames(Q) <- envt_type
  rownames(Q) <- envt_type

  matQ <- matrix(0, nrow=length(neighbor), ncol = length(neighbor))
  for(i in 1:(length(neighbor))){
    for (j in neighbor[[i]]){

      matQ[i,j] <- Q[R[i],R[j]]
    }
  }

  for (i in 1:(length(neighbor))){
    matQ[i,i] <- -sum(matQ[i,])
  }

  Lik<-tree_likelihood_only(Y,brs,matQ, m=lambda[(nb_param-nb_branch-1):(nb_param+nb_branch-2)], alpha = lambda[-2], beta = lambda[-1])

  return(Lik)

}


optimization_function2 <- function(lambda){

  Q <- matrix(c(lambda[1], lambda[1], lambda[2], lambda[2]), length(envt_type), length(envt_type))
  colnames(Q) <- envt_type
  rownames(Q) <- envt_type

  matQ <- matrix(0, nrow=length(neighbor), ncol = length(neighbor))
  for(i in 1:(length(neighbor))){
    for (j in neighbor[[i]]){

      matQ[i,j] <- Q[R[i],R[j]]
    }
  }

  for (i in 1:(length(neighbor))){
    matQ[i,i] <- -sum(matQ[i,])
  }

  Lik<-tree_likelihood_only(Y,brs,matQ, m=rep(1,length(Y)), alpha = 1, beta = 1)

  return(Lik)

}


# results <- c()
# for (i in 1:1000){
#
#   results <- c(results,optimization_function(i/1000))
#
# }


#### VERIFICATION OF THE PARAM OPTIMIZATION BY CREATING A LNL MAP AND THEN PERFORMING NLOPTR ####
# # Creation of a landscape map of likelihood
#
# lambda1 <- seq(10,100,5)
# lambda2 <- seq(10,100,5)
#
# #lambda1[1] <- 0.001
# #lambda2[1] <- 0.001
#
#
# likelihood <- data.frame(matrix(NA, length(lambda1), length(lambda2)))
#
# for (i in 1:length(lambda1)){
#   for (j in 1:length(lambda2)){
#     likelihood[j,i] <- optimization_function2(c(lambda1[i], lambda2[j]))
#   }
# }
#
#
#
#
# opts <- list("algorithm" = "NLOPT_GN_DIRECT_L_RAND",
#              "xtol_rel" = 1e-2,
#              maxeval = 500)
#
# params <- nloptr(x0 = c(0.001,10),
#                  lb = rep(0, 2),
#                  ub = rep(100, 2),
#                  eval_f = optimization_function2,
#                  opts = opts)
#
#
# opts <- list("algorithm" = "NLOPT_GN_DIRECT_L_RAND",
#              "xtol_rel" = 1e-5,
#              maxeval = 2000)
#
# params <- nloptr(x0 = params$solution,
#                  lb = rep(0, 2),
#                  ub = rep(100, 2),
#                  eval_f = optimization_function,
#                  opts = opts)
#
#
#
#
#
# library(plot3D)
# persp3D(x=seq(10,100,5), y=seq(10,100,5), z=as.matrix(likelihood), zlim = c(570, 625), theta=180, phi = 20)


