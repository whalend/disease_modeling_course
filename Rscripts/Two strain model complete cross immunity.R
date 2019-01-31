
library (deSolve) ## you need to load the package deSolve

### create a function that contains the ODE model, the function
### takes as argument t (time), x (state variables), and parms (parameters)

  sir.twostrains=function(t,x,parms){
  s=x[1]; i1=x[2]; i2=x[3]; r=x[4];
  
  ds=parms["mho"]-parms["beta1"]*i1*s-parms["beta2"]*i2*s-parms["mho"]*s
  di1=parms["beta1"]*i1*s-parms["gamma1"]*i1-parms["mho"]*i1-parms["m1"]*i1
  di2=parms["beta2"]*i2*s-parms["gamma2"]*i2-parms["mho"]*i2-parms["m2"]*i2
  dr=parms["gamma1"]*i1-parms["gamma2"]*i2-parms["mho"]*r
  dY=c(ds,di1,di2,dr)
  return(list(dY))
  
  }

## vector with the time sequence

times=seq(0, 1000, 0.1)

## vector containing the parameters

parms=c(beta1=0.12, beta2=0.1,gamma1=0.04, gamma2=0.05, mho=0.001, m1=0.01, m2=0.02)

## vector containing the initial values
x0=c(s=0.96, i1=0.02, i2=0.02, r=0)

Rnot1 = parms["beta1"]/(parms["gamma1"]+parms["mho"]+parms["m1"])

Rnot2 = parms["beta2"]/(parms["gamma2"]+parms["mho"]+parms["m2"])

# out is a dataframe containing the model simulated trajectory
out=as.data.frame ( lsoda(x0, times, sir.twostrains, parms))

plot (times, out$s, type="l", xlab="time", ylab="number",  col="green", ylim=c(0,1), main = "2-strain Complete Cross-Immunity Model")
lines (times, out$i1, type="l", col="red")
lines (times, out$i2, type="l", col="blue")
legend (100, 0.95, c("s", "i1","i2"),lty=c( 2, 2, 2), col=c("green","red", "blue"))
