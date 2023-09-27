###################
#### LIBRARIES ####
###################

library(diversitree)
library(deSolve)

##############
#### DATA ####
##############

y <- c(0,1)
tt <- seq(0,5, length=101)
pars <- c(.5,1)

###################
#### FUNCTIONS ####
###################

derivs.mk2new <- function(t,y,pars){

  D0 <- y[1]
  D1 <- y[2]

  q01 <- pars[1]
  q10 <- pars[2]

  dDdt <- c(-q01 * D0 + q01 * D1,
            -q10 * D1 + q10 * D0)

  list(dDdt)

}


make.branches.dtlik <- diversitree:::make.branches.dtlik
check.loaded.symbol <- diversitree:::check.loaded.symbol
check.control.ode <- diversitree:::check.control.ode

out <- lsoda(y,tt,derivs.mk2new, pars)[,-1]
info <- list(name="mknew", np=2, ny=2, idx.d=1:2, derivs=derivs.mk2new)
branches.mk2new <- make.branches.dtlik(info, check.control.ode())

###################
#### GRAPHICAL ####
###################

matplot(tt, out, type="l", las=1)
legend("topright", c("State 0", "State 1"), col=1:2, lty=1:2)
abline(h=pars[1]/sum(pars), lty=3)

