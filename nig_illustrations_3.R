# Create renormalised spa

library(ExactSPA)
n <- 1000
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)


x <- sort(rNIG(n, c(chi, psi, mu, gamma), seed = 123))
m <- mean(x)
sd <- sd(x)
hist(x, freq=F)

# Plot density
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
        sum(
            0.5*log(param[1]) +
                log(param[2] + param[4]^2) +
                sqrt(param[1]*param[2]) -
                log(pi) +
                log(besselK(y, -1, expon.scaled = TRUE)) - y +
                (data - param[3])*param[4] -
                log(y)
        )
    # Return the functional value
    return(loglik)
}

f_bessel <- exp(sapply(x, loglikNIG, param=par))
lines(x,f_bessel,col="red")
plot(x,f_bessel, type="l", lwd=2)
hist(x,freq=F, add=T)

# espa 
f_espa <- exp(- sapply(x, nll_fun_nig, par=par, type="ExactSPA"))
lines(x, f_espa, lty=4, col="red", lwd=2)

# spa
f_spa <- exp(- sapply(x, nll_fun_nig, par=par, type="SPA"))
lines(x, f_spa, lty=2, col="blue", lwd=2)

# respa
f_respa <- exp(-sapply(x, nll_fun_nig, par=par, type="reSPA"))
lines(x, f_respa, lty=5, col="green", lwd=2)
