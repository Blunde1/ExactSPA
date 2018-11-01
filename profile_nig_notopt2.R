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
nll_fun_nig <- function(par, X, type="ExactSPA"){
    lchi = par[1]
    lpsi = par[2]
    mu = par[3]
    gamma = par[4]
    if(type=="ExactSPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll
    }else if(type=="reSPA"){
        n <- 200
        x_sim <- sort(rNIG(n, c(exp(lchi), exp(lpsi), mu, gamma), seed = 1234))
        c <- sum(sapply(2:n, function(i) exp(-nll_nig(x_sim[i], lchi, lpsi, mu, gamma, 100, 512, 1)$nll)*(x_sim[i]-x_sim[i-1])))
        cat("value of renormalisation: ", c)
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll + length(X)*log(c)
    }else if(type=="Simpson"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 150, 512, 3)$nll
        
    }
    return(nll)
}

# Set parameters
#chi = 3.0e-4; psi=2e2; mu=-3e-4; gamma=2
#chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n = 100 #n=200
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
sqrt(VarX)


# Simulate a dataset
set.seed(4321)
x.nig <- rNIG(n,  c(chi, psi, mu, gamma), seed=123)

nll_fun_nig(par, X=x.nig)
loglikNIG(par, x.nig)
nll_fun_nig(par, X=x.nig, type="SPA")
nll_fun_nig(par, X=x.nig, type="reSPA")


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
par["gamma"] <- 100
nll_fun_nig(par, x.nig) 
nll_fun_nig(par, x.nig, type="reSPA") 
loglikNIG(par,x.nig)
nll_fun_nig(par,x.nig, type="Simpson")

gamma.val2 <- exp(seq(log(100), log(15000), length.out=20)) # second try 
#gamma.val2 <- seq(1000, 12000, length.out=20)
nll_espa <- nll_bessel <- nll_spa <- nll_respa <- numeric(length(gamma.val2))
for(i in 1:length(gamma.val2)){
    cat("iter: ",i,"\n")
    par["gamma"] = gamma.val2[i]
    nll_espa[i] <- nll_fun_nig(par, x.nig)
    nll_spa[i] <- nll_fun_nig(par, x.nig, type="SPA")
    nll_respa[i] <- nll_fun_nig(par, x.nig, type="reSPA")
    nll_bessel[i] <- loglikNIG(par, x.nig)
}

setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_besselfail_tail.pdf", width=7, height=4+1/3)
plot(gamma.val2, nll_espa, type="l", lwd=2, lty=4, log="x",
     ylab="Negative log-likelihood", xlab=expression(gamma))
points(gamma.val2, nll_bessel, col="red", pch=2, lwd=3)
points(gamma.val2, nll_spa, pch=3, col="blue", lwd=2)
#lines(gamma.val2, nll_respa, lty=5, col="green", lwd=2)
legend("topright", c("Saddlepoint adjusted IFT", "Bessel implementation", "Saddlepoint approximation"), 
       col=c("black","red", "blue"), pch=c(NA,2,3), lwd=c(2,3,2), lty=c(4,NA,NA))
dev.off()

# add spa
for(i in 1:length(gamma.val2)){
    cat("iter: ",i)
    par["gamma"] = gamma.val2[i]
    nll_respa[i] <- nll_fun_nig(par, x.nig, type="reSPA")
    #    nll_bessel[i] <- loglikNIG(par, x.nig)
}
plot(gamma.val2, nll_respa, lty=5, col="green", lwd=2)


# start opt parameters
# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n = 100

# Simulate a dataset
set.seed(4321)
x.nig <- rNIG(n,  c(chi, psi, mu, gamma), seed=123)
nll_nig(x.nig, par[1],par[2],par[3],par[4], 400, 128, 3)$nll


# testing when spa renorm is close to 1
#par["lpsi"] = log(1000)
#par["lchi"] = log(0.000001)
#par["mu"] = -0.01
par["gamma"] <- 1
nll_fun_nig(par, x.nig, type="reSPA")
loglikNIG(par, x.nig)
nll_fun_nig(par, x.nig, type="Simpson")

# GAMMA
gamma.val3 <- seq(1,150, length.out=25)
nll_espa <- nll_bessel <- nll_spa <- nll_respa <- nll_simpson <- numeric(length(gamma.val3))
for(i in 1:length(gamma.val3)){
    cat("iter: ",i, "\n")
    par["gamma"] = gamma.val3[i]
    nll_espa[i] <- nll_fun_nig(par, x.nig)
    nll_spa[i] <- nll_fun_nig(par, x.nig, type="SPA")
    #nll_respa[i] <- nll_fun_nig(par, x.nig, type="reSPA")
    nll_bessel[i] <- loglikNIG(par, x.nig)
    nll_simpson[i] <- nll_fun_nig(par, x.nig, type="Simpson")

}

# MU
par["gamma"] <- 2
mu.val3 <- seq(-0.1, 0.1, length.out=25)
nll_espa_mu <- nll_bessel_mu <- nll_spa_mu <- nll_simpson_mu <- numeric(length(mu.val3))
for(i in 1:length(mu.val3)){
    cat("iter: ",i, "\n")
    par["mu"] = mu.val3[i]
    nll_espa_mu[i] <- nll_fun_nig(par, x.nig)
    nll_spa_mu[i] <- nll_fun_nig(par, x.nig, type="SPA")
    nll_bessel_mu[i] <- loglikNIG(par, x.nig)
    nll_simpson_mu[i] <- nll_fun_nig(par, x.nig, type="Simpson")
}


library(latex2exp)
setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_spafail_iftfail3.pdf", width=15, height=7)

par(mfrow=c(1,2))

plot(gamma.val3, nll_espa, type="l",lwd=4,
     ylim=range(nll_spa,nll_espa,nll_simpson), 
     ylab="Negative log-likelihood", xlab=expression(gamma))
lines(gamma.val3, nll_bessel, type="p",col=2, pch=2, lwd=2)
lines(gamma.val3, nll_spa, type="l", lty=1, col=3, lwd=2)
lines(gamma.val3, nll_simpson, type="l", lty=2, col=4, lwd=2)
legend("topleft", c("Saddlepoint adjusted IFT", "Bessel implementation", "Saddlepoint approximation", "Simpson IFT"),
       col=1:4, pch=c(NA,2,NA,NA), lwd=c(4,2,2,2), lty=c(1,NA,1,2))

plot(mu.val3, nll_espa_mu, type="l",pch=3, lwd=4,
     ylab="Negative log-likelihood", xlab=TeX("$\\mu$"))
lines(mu.val3, nll_bessel_mu, type="p",col=2, pch=2, lwd=2)
lines(mu.val3, nll_spa_mu, type="l", lty=1, col=3, lwd=2)
lines(mu.val3, nll_simpson_mu, type="l", lty=2, col=4, lwd=2)
legend("bottomleft", c("Saddlepoint adjusted IFT", "Bessel implementation", "Saddlepoint approximation", "Simpson IFT"),
       col=1:4, pch=c(NA,2,NA,NA), lwd=c(4,2,2,2), lty=c(1,NA,1,2))

dev.off()

