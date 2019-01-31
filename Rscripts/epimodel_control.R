#' # Epidemiological Models for Disease Control
#' 
#' Generic code for simulating epidemiology with control
#' 
#+
library(deSolve)

parms <- c(beta = 1e-3, gamma = 1e-1)
inits <- c(S1 = 499, I1 = 1, R1 = 0)
dt <- seq(0, 100, 0.1)

epi.sim <- function(t, x, parms){
      with(as.list(c(parms, x)), {
            dS <- -(beta * S1 * I1)
            dI <- beta * S1 * I1
            dR <- gamma * I1
            der <- c(dS, dI, dR)      
            list(der)# return the output
      })# end 'with'
}# end function

sim.1 <- as.data.frame(lsoda(inits, dt, epi.sim, parms = parms))
plot.ts(sim.1$I, lwd = 2, ylim = c(0, 1000))


