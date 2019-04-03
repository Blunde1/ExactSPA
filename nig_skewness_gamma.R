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

chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n = 100 #n=200
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
sqrt(VarX)

par(mfrow=c(1,2))
library(latex2exp)
par["gamma"] <- 2
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
x <- seq(EX-8*sqrt(VarX), EX+8*sqrt(VarX), length.out=50)
y <- exp(-sapply(x, loglikNIG, param=par))
plot(x,y, main=TeX("$\\gamma = 2$"))

par["gamma"] <- 150
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
x <- seq(EX-8*sqrt(VarX), EX+8*sqrt(VarX), length.out=50)
y <- exp(-sapply(x, loglikNIG, param=par))
plot(x,y, main=TeX("$\\gamma = 150$"))


# IG variance = 10
par(mfrow=c(1,1))
vw <- 0.1
chi <- psi <- 1/vw
mu <- gamma <- 0
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
x <- seq(EX-8*sqrt(VarX), EX+12*sqrt(VarX), length.out=300)
y <- exp(-sapply(x, loglikNIG, param=par))
plot(x,y, main=TeX("$\\theta = 10$"))
