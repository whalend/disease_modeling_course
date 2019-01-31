#' # Simulating a Multi-Population Disease
#' ### Generic code for simulation of a multi-population disease model
#' 
#' ## Initialize Modeling Environment
#+ load package
library(deSolve)
#' This package has general solvers for differential equations.
#'
#' Set up vectors of parameters, initial populations, and the time step.
#' 
#+ set parameters and initial values
parms <- c(beta1=1e-3, gamma1=1e-1, beta2=1e-4, gamma2=1e-1)
inits <- c(S1=499,I1=1 ,S2=1990, I2=10)# H for host, V for vector
dt <- seq(0,100,0.1)# time step
#' Note that there are separate beta and gamma values that correspond to the two different populations of susceptibles and infecteds.
#'
#' Write the simulation model function
#+ epidemic simulation model function
epi.sim <- function(t, x, parms){
      with(as.list(c(parms,x)),{
            
            dS1 <- -beta1*S1*I1
            dI1 <- beta1*S1*I1 - gamma1*I1
            dS2 <- -beta2*S2*I2
            dI2 <- beta2*S2*I2 - gamma2*I2
                  
            der <- c(dS1, dI1, dS2, dI2)
            list(der) # the output must be returned
      }) # end of ’with’
} # end of function definition
#'
#' Run a simulation using the initial parameters
#+ run simulation
sim1 <- as.data.frame(lsoda(inits, dt, epi.sim, parms = parms))
#'
#' Plot the results of the simulation run
#+ plot simulation run
plot.ts(sim1$I1,lwd=2,ylim=c(0,1000))
lines(c(1:1001),sim1$I2,col="blue",lty=2,lwd=2)
leg.txt<-c("Infected Group 1","Infected Group 2")
legend("topright",leg.txt,lty=c(1,2),col=c("black","blue"))
