#' # Gillespie's Algorithm Sample Code for SIS Model
#' 
#' ## The Gillespie Function
#' 
#+
gillesp <- function(start, ratefun, trans, pars, times){
      t0 <- times[1] # set time to starting time
      ntimes <- length(times) # totoal time duration
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
start <- c(S = 99, I = 1)
#'
#' Specify rate functions for infection/transmission and recovery
#+
ratefun.SIS <- function(X, pars, time){
      vals = c(as.list(pars), as.list(X)) # attache state and pars as lists
      rates = with(vals, #allows reference to states nd parameters by name
                   c(infection = beta * S * I, 
                     recovery = gamma * I))
}

statenames.SIS <- c("S", "I") # state variable names
transnames.SIS <- c("infection", "recovery") # transition names
#'
#' Define the transition matrix
#+
trans.SIS <- matrix(c(-1, 1, 1, -1),
                    byrow = T, # defaults to column
                    ncol = 2, # equal to number of state variables
                    dimnames = list(transnames.SIS, statenames.SIS))
#'
#' Specify parameters
#+
pars.SIS <- c(beta = 0.05, gamma = 1)
#'
#' Specify simulation time
#+
times <- seq(0, 5, by = 0.05)
#'
#' ## Run simulation
#+
G.SIS.mult <- replicate(100, gillesp(start = start, times = times, ratefun = ratefun.SIS, trans = trans.SIS, pars = pars.SIS)[,"I"])

matplot(times, G.SIS.mult, type = "l", col = "gray", lty = 1, xlab = "Time", ylab = "Number infectious")
# add line showing mean value of infectious from all simulations
lines(times, rowMeans(G.SIS.mult), lwd = 2)
#'
#' ### Why is the mean (solid black line) substantially lower than most of the simulations (grey lines)?
#' - it's being pulled down by the simulations that have very low values
#' 
#' 
#' # Gillespie's Algorithm Sample Code for SIR Model
#' 
#' ## The Gillespie Function
#' 
#+
gillesp <- function(start, ratefun, trans, pars, times){
      t0 <- times[1] # set time to starting time
      ntimes <- length(times) # totoal time duration
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
#'
#' Define starting conditions, 1 infected and 99 susceptibles
#+
start <- c(S = 99, I = 1)
#'
#' Specify rate functions for infection/transmission and recovery
#+
ratefun.SIS <- function(X, pars, time){
      vals = c(as.list(pars), as.list(X)) # attache state and pars as lists
      rates = with(vals, #allows reference to states nd parameters by name
                   c(infection = beta * S * I, recovery = gamma * I))
}

statenames.SIS <- c("S", "I") # state variable names
transnames.SIS <- c("infection", "recovery") # transition names
#'
#' Define the transition matrix
#+
trans.SIS <- matrix(c(-1, 1, 1, -1),
                    byrow = T, # defaults to column
                    ncol = 2, # equal to number of state variables
                    dimnames = list(transnames.SIS, statenames.SIS))
#'
#' Specify parameters
#+
pars.SIS <- c(beta = 0.05, gamma = 1)
#'
#' Specify simulation time
#+
times <- seq(0, 5, by = 0.05)
#'
#' ## Run simulation
#+
G.SIS.mult <- replicate(100, gillesp(start = start, times = times, ratefun = ratefun.SIS, trans = trans.SIS, pars = pars.SIS)[,"I"])

matplot(times, G.SIS.mult, type = "l", col = "gray", lty = 1, xlab = "Time", ylab = "Number infectious")
# add line showing mean value of infectious from all simulations
lines(times, rowMeans(G.SIS.mult), lwd = 2)

