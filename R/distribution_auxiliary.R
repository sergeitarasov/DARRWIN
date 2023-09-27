'
Script developed by: Thomas Merrien
creation: 05/05/2023
last update: 05/05/2023
File name: distribution_auxiliary.R

Description:
Creating personal function to calculate density of various distribution in order
to speed up calculation
'


#' Density of gamma distribution
#'
#' Function that given the parameter of the gamma function and the point to estimate,
#' will return the density of the gamma function at this point. All the terms
#' should be strictly positives.
#' The formula of the distribution is the following one:
#' **f(x, \Alpha, \Beta) = frac{x^{\Alpha - 1} . e{- \Beta x} . \Beta^{\Alpha}{\Gamma (\Alpha)}}**
#' where * \Gamma (\Alpha) = (\Alpha - 1)! *
#'
#' @param x point(s) where to estimate the gamma density (float or vector)
#' @param alpha alpha parameter of the gamma distribution (float)
#' @param beta beta parameter of the gamma distribution (float)
#'
#' @return the value(s) of the gamma density estimated on x
#'

dens_gamma <- function(x, alpha, beta){

  return (x**(alpha -1) * exp(-beta * x) * beta**(alpha)) / factorial(alpha - 1)

}

