'
Script developed by: Ignacio Quintero adapted to R by Thomas Merrien
last update: 22/04/2022
File name: reps_per_period.R
'

#'Description:
#'reps_per_period(br_length, const_dt)
#'
#'Estimate number of repetitions for each speciation waiting time. Function adapted
#'from TRIBE developed by Ignacio Quintero
#'
#'Input:
#'-br_length: branch length (numeric)
#'-const_dt: small and constant timestep (numeric)
#'
#'
reps_per_period <- function(br_length, const_dt){

  return (round(br_length/const_dt))

}

