#' # CBS 810 Infectious Disease Modeling Homework
#' 
#' Whalen Dillon
#' 
#' 12/7/2015
#' 
#' ## Stochastic Simulation w/Disease-Induced Death Gillespie's Algorithm Sample Code for SIS Model
#' 
#' ### The Gillespie Function
#' 
#+
gillesp <- function(start, ratefun, trans, pars, times){
      t0 <- times[1] # set time to starting time
      ntimes <- length(times) # total time duration
      X <- start # set state to starting state
      res <- matrix(nrow = length(times), ncol = length(start), dimnames = list(times, names(start))) # matrix for results
      for(ctr in 1:(ntimes-1)){# loop over reporting times
            res[ctr,] <- X # record the current state
            while(t0 < times[ctr + 1]){
                  rates <- ratefun(X, pars, t0) # calculate current rates
                  if(all(rates == 0)) break # extinction
                  totrate <- sum(rates)
                  elapsed <- rexp(1, totrate) # sample elapsed time
                  which.trans <- sample(1:nrow(trans), size = 1, prob = rates) # pick the transition
                  t0 <- t0 + elapsed # update time
                  X <- X + trans[which.trans,] # add transition values to current state
            }
      }
      cbind(times, res)
}

#' Define starting conditions, 1 infected and 99 susceptibles
#+
start <- c(S = 99, I = 1, R = 0)
#'
#' Specify rate functions: beta*S*I & gamma*I for infection/transmission and recovery
#+
ratefun.SIR <- function(X, pars, time){
      vals = c(as.list(pars), as.list(X)) # attach state and pars as lists
      rates = with(vals, #allows reference to states and parameters by name
                   c(infection = beta * S * I, 
                     recovery = gamma * I, 
                     birth = mu * S, 
                     death = mu * I
                     )
                   )
}

statenames.SIR <- c("S", "I", "R") # state variable names
transnames.SIR <- c("infection", "recovery", "birth", "death") # transition names
#'
#' Define the transition matrix
#+
trans.SIR <- matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, 0, 0, -1, 0),
                    byrow = T, # defaults to column
                    ncol = 3, # equal to number of state variables
                    dimnames = list(transnames.SIR, statenames.SIR))
#'
#' Specify parameters
#+
pars.SIR <- c(beta = 0.07, gamma = 1, mu = 0.1)
#'
#' Specify simulation time
#+
times <- seq(0, 5, by = 0.05)
#'
#' ## Run simulation
#+
G.SIR.mult <- replicate(100, gillesp(start = start, times = times, ratefun = ratefun.SIR, trans = trans.SIR, pars = pars.SIR)[,"I"])

matplot(times, G.SIR.mult, type = "l", col = "gray", lty = 1, xlab = "Time", ylab = "Number infectious")
# add line showing mean value of infectious from all simulations
lines(times, rowMeans(G.SIR.mult), lwd = 2)
