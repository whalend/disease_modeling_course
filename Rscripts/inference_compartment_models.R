#' # Inference from Compartment Models
#' **CBS810 Lab 3**
#' Learning Objectives:
#'  - Undersand the relationship between "simulation" and "inference"
#'  - Know the basics of Sum of Square Error (SSE) and Least Square Estimator(LSE)
#'  - Estimate the LSE of compartment epidemic models
#'  
#' **Simulation:** given known paramters, get simulated observations
#' **Inference:** given actual observations, get optimized parameters
#' 
#' #### Example of inference
#' Linear Model: f(x) = ax + b
#' Residual Error: $r_i = y_i - f(x_i | beta(a,b))$
#' SSE: S = $sigma^n_i=1 r^2_i$
#' LSE: estimated parameter value of $beta(a,b)$ that minimizes S
#' 
#' **Linear Regression**
#+
x <- seq(1,10)
y <- x + rnorm(10, 0.1)
plot(x, y, type = "b", ylim = c(-2, 12))
fit1 <- lm(y ~ x)
abline(fit1, col = "blue")
fit1$coefficients

#' #### Inference from the real compartment model
#'  1. Specify model structure (e.g. SIS or SIR) to determine how many parameters need to be estimated
#'  2. Write (borrow) the existing simulation function
#'  3. Calcualte SSE
#'  4. Minimize SSE to get optimized value of LSE
#'
#' Additional Reference: "Fitting Epidemics in R" by John M. Drake  
#+
library(deSolve)# package for analyzing ordinary differential equations(ODE)
load("data/fluday.Rdat")# load sample observation data
summary(flu.day)
plot(sim.sis.I ~ day, data = flu.day, type = "b", xlab = "Day", ylab = "I(t)")

#' The SIS Model
#+
sis.model <- function(t, x, params){
      S <- x[1]
      I <- x[2]
      
      beta <- params[1]
      gamma <- params[2]
      
      dS <- -beta * S * I + gamma * I
      dI <- beta * S * I - gamma * I
      
      list(c(dS, dI))
}

#' The data and SSE function for the SIS model
#+
sse.sis <- function(params0,data){
      t <- data[, 1] # Do not confuse with x[1] and x[2]
      cases <- data[, 2]
      
      beta <- params0[1]# initial value for beta
      gamma <- params0[2]# initial value for gamma
      
      S0 <- 499 # initial number of S
      I0 <- 1   # initial number of I
      
      out <- as.data.frame(ode(y = c(S = S0, I = I0), times = t, sis.model, parms = c(beta, gamma)))
      sse <- sum((out$I - cases)^2) # This is the SSE
}


#' The actual minimization to get the optimized LSE 
params0 <- c(0.005, 0.15)
fit0 <- optim(params0, sse.sis, data = flu.day); fit0$par
fit1 <- optim(fit0$par, sse.sis, data = flu.day); fit1$par
