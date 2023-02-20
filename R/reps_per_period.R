'
Script developed by: Ignacio Quintero adapted to R by Thomas Merrien
last update: 07/02/2023
File name: reps_per_period.R
'


#' Number of time slots per period
#'
#' This function estimates the number of repetitions for each speciation waiting time. Function adapted
#' from TRIBE developed by Ignacio Quintero
#'
#' @param br_length branch length (numeric)
#' @param const_dt small and constant timestep (numeric)
#'
#' @return Returns the number of time period that will happen on one branch (numeric)
#' @export
#'
#' @examples
#' #For a branch of length 1 with a delta time of 0.05
#' reps_per_period(1, 0.05)


reps_per_period <- function(br_length, const_dt){

  return (round(br_length/const_dt))

}

