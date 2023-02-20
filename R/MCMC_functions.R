'
Script developed by: Thomas Merrien
Date of creation: 07/02/2023
last update: 08/02/2023
File name: MCMC_infer.R
'

#' Hastings ratio MCMC
#'
#' Function that calculate the Hastings ratio (acceptance ratio) for MCMC parameters
#' selection
#'
#' @param Lp Likelihood of the data with the new (proposed) parameter
#' @param Lc Likelihood of the data with the current parameter
#' @param Qp Probability of drawing in the parameter distribution the new (proposed) parameter knowing the current one
#' @param Qc Probability of drawing in the parameter distribution the current parameter knowing the proposed one
#'
#' @return the Hastings ratio (numeric)
#'

H_ratio <- function(Lp,Lc,Qp,Qc){
  return ((Lp*Qc)/(Lc*Qp))
}



#' Log of the Hastings ratio MCMC
#'
#' Function that calculate the log of the Hastings ratio (acceptance ratio) for MCMC parameters
#' selection
#'
#' @param LLp Log-Likelihood of the data with the new (proposed) parameter
#' @param LLc Log-Likelihood of the data with the current parameter
#' @param LQp Log-probability of drawing in the parameter distribution the new (proposed) parameter knowing the current one
#' @param LQc Log-probability of drawing in the parameter distribution the current parameter knowing the proposed one
#'
#' @return the log of the Hastings ratio (numeric)
#'

LH_ratio <- function(LLp,LLc,LQp,LQc){
  return (LLp-LLc + LQc-LQp)
}




#' Research of the tuning parameter for a MCMC with several parameters
#'
#' Function that will perform some MCMC with different tuning parameters to look for the optimal one in a context of multiple or single parameters MCMC
#'
#' @param tuning initial tuning parameter (numeric)
#' @param parameters list of the parameters to estimate (list)
#'
#' @return the right tuning parameter for the MCMC (numeric)
#'

single_tuning <- function(tun, parameters){

  d <- length(parameters) #number of parameters to estimate
  acc_rate <- 0 #set acceptance rate at 0 at the beginning
  while ((0.234 <= acc_rate) && (acc_rate <= 0.44)==FALSE){ #we want acceptance rate to fall within the 0.234-0.44 interval. 0.234 being optimal when we have around 5 parameters and 0.44 being optimal when we have 1 parameter

    tun_test <- rnorm(1, tun, 0.1) #propose a new tuning parameter
    while (tun_test <= 0){ #the tuning parameter has to be positive
      tun_test <- rnorm(1, tun, 0.1)
    }

    ratio <- c() #list of acceptance ratio
    LLc #Calculate likelihood of the model with current parameter values
    for (i in 1:(ceiling(100/d)*d)){ #let's see what is the acceptance rate we get with the new tuning parameter

      #Run MCMC iteration
      MCMC_single <- MCMC_SingleChain(parameters = parameters, tuning = tun, parameter_to_update = (i%%d))
      parameters_p <- MCMC_single[1]
      LLp <- MCMC_single[2] #Likelihood of the model with proposed parameter values
      LQp
      LQc

      ratio[i] <- LH_ratio(LLp,LLc,LQp,LQc) #store the Hastings ratio in a list

      if (min(0,ratio[i]) < -rexp(1)){
        ratio[i] <- 1
        parameters <- parameters_p
        LLc <- LLp
      } else {
        ratio[i] <- 0
      }


    }

    acc_rate <- sum(ratio)/(ceiling(100/d)*d)

  }
}


#' Realization of a single chain of MCMC
#'
#' Function that will perform a single chain of a MCMC
#'
#' @param parameters list of the parameters value of your process.
#' @param tuning tuning parameter (numeric)
#' @param parameter_to_update index in the parameter list of the parameter to update (numeric)
#'
#' @return the right tuning parameter for the MCMC (numeric)
#'

MCMC_SingleChain <- function(parameters, tuning, parameter_to_update){

  ö

}

