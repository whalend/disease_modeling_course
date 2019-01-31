#' If needed, set working directory:
# setwd("~")  ### YOU NEED YOU CHANGE TO YOUR OWN WORKING DIRECTORY

#' Read in mumps sero- data
#+
mumps <- read.csv("data/Seroprofilemumps.csv",header=T)
summary(mumps)

#' Create new variables in the dataframe mumps
#+
mumps$seroprevalence <- mumps$Seropositive/(mumps$Seropositive + mumps$Seronegative)
summary(mumps)

#' plot proportion of (ever) infected (i.e., seroprevalence) as function of age
#+
plot(mumps$Age, mumps$seroprevalence, xlab='Age', ylab='Seroprevalence')

#' fit the catalytic model to the seroprevalence data
#+
model <- nls(seroprevalence ~ 1 - exp(-lambda * Age), data = mumps, start = list(lambda = 0.40), trace = T)

#' to see the fit results 
#+
summary(model)

#' The average age to infection (A) ≈ 1/lambda = `r 1/0.204`
#' Why measles and mumps are called childhood diseases(?)

#' plot data against predicted model
#+
plot(mumps$Age, mumps$seroprevalence, xlab='Age', ylab='Seroprevalence')
lines(mumps$Age, 1-exp(-mumps$Age * summary(model)$coefficients[1]), 
      lwd=2, col="red")

#' The prediction from the model fits the data pretty closely.

#' Fitting age-dependent force of infections
#' explore whether the data show age-dependency in the force of infection (plot -log(sa))

#+
plot(mumps$Age, -log(1-mumps$seroprevalence), 
     xlab='Age', ylab='-log(Seronegatives)')

#' fit a linear model to the -log(Seronegatives) to estimate the two forces of infection
#+
fit <- lm(-log(1 - mumps$seroprevalence[1:15]) ~ mumps$Age[1:15])
abline(fit)

fit2 <- lm(-log(1 - mumps$seroprevalence[16:32]) ~ mumps$Age[16:32])
abline(fit2)

#' Note: I am using until de data point 32 to avoid log(0), some data points after 32 have 100 % seropositivity
#' 
#' ## Estimate the basic R-nought from seroprevalence data and model
#' 
#' R-nought = 1 + FOI(lambda) / mortality rate(mu)
#' For an assumed life-expectancy of 75 years:
#+
summary(model)
rnought_mumps = 1 + summary(model)$coefficients[1]/(1/75)

#' ## Rubella Data
#' 
#+
rubella <- read.csv("data/RubellaBangladesh.csv")
summary(rubella)
rubella

rubella$seroprevalence <- (rubella$Number.tested - rubella$Number.negative)/rubella$Number.tested
rubella

plot(seroprevalence ~ Age, data = rubella)

rub_model <- nls(seroprevalence ~ 1 - exp(-lambda * Age), data = rubella, start = list(lambda = 0.40), trace = T)
summary(rub_model)

lines(rubella$Age, 1 - exp(-rubella$Age * summary(rub_model)$coefficients[1]), 
      lwd=2, col="red")

#' Fewer observations, so the prediction from the model doesn't fit the obesrved values as well as they did for the mumps data
#' 
#' #### Estimate basic reproduction number
#' Assuming life expectancy of 68 years:
#' 
#+
rnought_rubella = 1 + summary(rub_model)$coefficients[1]/(1/68)

#' #### Is force of infection age dependent?
#'
#+
fit1 <- lm(-log(1 - seroprevalence[1:15]) ~ Age[1:15], data = rubella)
abline(fit1)

fit2 <- lm(-log(1 - seroprevalence[16:32]) ~ Age[16:32], data = rubella)
abline(fit2)
