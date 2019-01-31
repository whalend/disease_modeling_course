# ODE vector-borne disease system 
library("deSolve")
parms <- c(beta1=1e-3, gamma=1e-1, beta2=1e-4, mu=1e-1)
inits <- c(SH=499,IH=1,RH=0,SV=1990,EV=0,IV=10) # H for host, V for vector
dt <- seq(0,100,0.1)

VBD <- function(t, x, parms){
      with(as.list(c(parms,x)),{
            dSH <-  -beta1*IV*SH
            dIH <-  +beta1*IV*SH - gamma*IH
            dRH <-  gamma*IH
            dSV <- -beta2*IH*SV
            dEV <-  beta2*IH*SV - mu*EV
            dIV <-  mu*EV
            der <- c(dSH, dIH, dRH, dSV, dEV, dIV)
            list(der) # the output must be returned
            }) # end of ??with??
} # end of function definition

sim.vbd <- as.data.frame(lsoda(inits, dt, VBD, parms=parms))
plot.ts(sim.vbd$IV,lwd=2,ylim=c(0,1000))
lines(c(1:1001),sim.vbd$IH,col="blue",lty=2,lwd=2)
leg.txt<-c("Infectious Vector","Infected Host")
legend("topright",leg.txt,lty=c(1,2),col=c("black","blue"))
