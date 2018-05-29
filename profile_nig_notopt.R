# Script for plotting profile likelihood (not in optimum) for NIG RV
# Berent Å. S. Lunde
# 18.05.2018

library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")
?ExactSPA::hello()

# BESSEL ####
loglikNIG=function(param,data){
    #This routine requires the observations to be stored in
    # a global variable "data"
    param[1]<-exp(param[1]) #Reparametrization
    param[2]<-exp(param[2])
    # param[1] : chi
    # param[2] : psi
    # param[3] : mu
    # param[4] : gamma
    # The square root expression
    y = sqrt((param[1] + (data - param[3])^2)*(param[2] + param[4]^2))
    # The log-likelihood
    loglik =
        -sum(
            0.5*log(param[1]) +
                log(param[2] + param[4]^2) +
                sqrt(param[1]*param[2]) -
                log(pi) +
                log(besselK(y, -1, expon.scaled = F)) + # - y + # - y if expon.scaled=T
                (data - param[3])*param[4] -
                log(y)
        )
    # Return the functional value
    return(loglik)
}

# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n = 200

# Simulate a dataset
set.seed(4321)
attach(par)
x.nig <- rNIG(n,  c(chi, psi, mu, gamma), seed=123)

nll_fun_nig(par, X=x.nig)
loglikNIG(par, x.nig)

# Look at distribution of parameters - from sim exper, to get an idea of parameters
par.eval <- c(lchi=log(0.01),
               lpsi = log(5000),
               mu = - 0.02,
               gamma = -30)

loglikNIG(par.start, x.nig)
nll_fun_nig(par.start, x.nig)

# Find initial area
gamma.val <- seq(-15000,15000, length.out=500)
nll_bessel <- numeric(length(gamma.val))
par["lpsi"] = log(5e7)
par["lchi"] = log(0.000001)
par["mu"] = -0.01
for(i in 1:length(gamma.val)){
    par["gamma"] = gamma.val[i]
    nll_bessel[i] <- loglikNIG(par, x.nig)
}
plot(gamma.val, nll_bessel)
# bessel gets problems around gamma > 10k

# Can espa handle it?
par["gamma"] <- 12000
nll_fun_nig(par, x.nig) # ok
par["gamma"] <- 15000
nll_fun_nig(par, x.nig) # great

gamma.val2 <- seq(1000, 15000, length.out=50)
nll_espa <- nll_bessel <- numeric(length(gamma.val2))
for(i in 1:length(gamma.val2)){
    cat("iter: ",i)
    par["gamma"] = gamma.val2[i]
    nll_espa[i] <- nll_fun_nig(par, x.nig)
    nll_bessel[i] <- loglikNIG(par, x.nig)
}
plot(gamma.val2, nll_espa, type="o",pch=3, lwd=2,
     ylab="Log-likelihood", xlab=expression(gamma))
lines(gamma.val2, nll_bessel, type="b", col="red", pch=2, lwd=2)
legend("topleft", c("Exact SPA", "Bessel"), col=c("black","red"), pch=c(3,2), lwd=c(2,2))
