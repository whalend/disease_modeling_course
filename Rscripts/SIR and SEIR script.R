
harmonic <- function(x){y=1/mean(1/x)
y
}

doggies <- function(x){y=(6*x)
y
}

#SIR Model
install.packages("deSolve")
library(deSolve)

sir.model <- function(t,x,params){
  S <- x[1] # Susceptible
  I <- x[2] # Infected
  R <- x[3] # Recovered
  
  beta <- parms[1] # Transmission rate
  gamma <- parms[2] # Recovery rate
  
  dS <- -beta*S*I
  dI <- beta*S*I - gamma*I
  dR <- gamma*I
  
  list(c(dS, dI, dR))
}

parms <- c(beta = 1e-3, gamma = 1e-1)
inits <- c(S = 499, I = 1, R = 0)
dt <- seq(0, 100, 0.1)

sim.2 <- as.data.frame(lsoda(inits, dt, sir.model, parms=parms))
plot.ts(sim.2$I, lwd=2, ylim=c(0,500), ylab="Number", main="SIR model")
lines(c(1:1001), sim.2$S, col="blue", lty=2, lwd=2)
lines(c(1:1001), sim.2$R,col="green", lty=2, lwd=2)
leg.txt<-c("Infected", "Susceptible", "Recovered")
legend("topright", leg.txt, lty = c(1,2,3), col = c("black", "blue", "green"))


# SEIR Model
seir.model <- function(t, x, params){
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  
  beta <- parms[1] # Transmission rate
  gamma <- parms[2] # Recovery rate
  infect <- parms[3]
  mu <- parms[4]
  
  
  dS <- -beta*S*I
  dE <- beta*S*I - infect*E - mu*E
  dI <- infect*E - gamma*I - mu*I
  dR <- gamma*I - mu*R
  
  list(c(dS, dE, dI, dR))
}

parms <- c(beta=1e-3, gamma=1e-1, infect=.25, mu=.05)
inits <- c(S=499, E=0, I=1, R=0)
dt <- seq(0,100,.1)

sim.3 <- as.data.frame(lsoda(inits, dt, seir.model, parms=parms))
plot.ts(sim.3$I,lwd=2,ylim=c(0,500), ylab="Number", main="SEIR model")
lines(c(1:1001), sim.2$S, col="blue", lty=2, lwd=2)
lines(c(1:1001), sim.2$E, col="red", lty=2, lwd=2)
lines(c(1:1001), sim.2$R,col="green", lty=2, lwd=2)
leg.txt<-c("Infected", "Susceptible","Exposed", "Recovered")
legend("topright", leg.txt, lty=c(1,2,3,4), col=c("black", "blue", "red", "green"))


# 1. Try to derive the steady-state (equilibrium) for the SIR model. Compare your result with the simulation figure, do they align well?

# 2. Change the parameters in the model and see how the result changes accordingly!
      
# 3. Modify the SIR code to accommodate the epidemic SEIR model 

# (HINT: use your homework’s solution as a guide for the codes!)
# You can try endemic SEIR model if you want some challenge!
      
