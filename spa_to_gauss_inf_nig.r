# illustrate spa --> Gaussian under optimisation
# INFINITE SAMPLING
# illustrate with NIG
# Berent Lunde
# 27.08.2018

# # CONDITIONS # #
# EW = 1 --> sqrt(chi) = sqrt(psi)
# mu = gamma = 0
# VarW = sqrt(chi) / psi^(3/2) = 1/ psi = 1/chi

# # TO DO:
# Function that takes in $\theta
# Outputs negative int log(spa(f(theta); x)) p_X(x) dx

# Minimise the function over different theta

# Call library
library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")
ExactSPA::hello()


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

c(100, 512) # IG: 0.1, 1
c(300, 1024) # IG: 3, 5
c(1000, 2048) # IG: 10
# Likelihood wrappers
nll_nig_locked <- function(ltheta, X, type="ExactSPA"){
    VarW <- exp(ltheta)
    lchi = -log(VarW)
    lpsi = -log(VarW)
    mu = 0
    gamma = 0
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
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 1500, 252, 3)$nll
        
    }
    return(nll)
}

# normal for reference to spa --> gaussian under optimisation
nll_normal <- function(par, X){
    # par[1]: mu
    # par[2]: log(sigma)
    mu <- par[1]
    sigma <- exp(par[2])
    nll <- -sum(dnorm(X, mean=mu, sd=sigma, log=T))
    return(nll)
}

# Integrand
nll_integrand <- function(x, ltheta_est, ltheta_true, type="ExactSPA"){
    
    if(type=="StandardGauss"){
        res <- nll_normal(c(0,0), X=x) * exp(-nll_nig_locked(ltheta = ltheta_true, X=x, type="ExactSPA"))
    }else{
        res <- nll_nig_locked(ltheta=ltheta_est, X=x, type=type) * exp(-nll_nig_locked(ltheta = ltheta_true, X=x, type="ExactSPA"))   
    }
    
    return(res)

}
nll_integrand(0, log(2), log(2), "ExactSPA")


nll_expected <- function(ltheta_est, ltheta_true, type="ExactSPA"){
    
    # mu +- range
    vw <- exp(ltheta_true)
    chi <- psi <- 1/vw
    mu <- gamma <- 0
    par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
    EX <- mu + sqrt(chi/psi)*gamma
    VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
    lower <- EX - 12*sqrt(VarX)
    upper <- EX + 12*sqrt(VarX)
    
    # # Simulate to find range
    # nsim <- 10000
    # vw <- exp(ltheta_true)
    # chi <- psi <- 1/vw
    # mu <- gamma <- 0
    # par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
    # x.nig <- rNIG(nsim,  c(chi, psi, mu, gamma), seed=123)
    # lower <- min(x.nig)
    # upper <- max(x.nig)
    
    # Calculate integrand over range
    m <- 200
    # dx <- (upper-lower)/m
    # val <- numeric(m)
    # for(i in 0:(m-1)){
    #     val[i+1] <- nll_integrand(x=lower + i*dx,
    #                               ltheta_est = ltheta_est,
    #                               ltheta_true=ltheta_true,
    #                               type=type) *
    #         dx
    # }
    
    h <- (upper-lower)/m
    xj <- seq.int(lower, upper, length.out = m + 1)
    xj <- xj[-1]
    xj <- xj[-length(xj)]
    val <- numeric(m+1)
    
    val[1] <- nll_integrand(x=lower, ltheta_est = ltheta_est, ltheta_true=ltheta_true, type=type)
    for(i in seq.int(2, length(xj), 2)){
        val[i+1] <- 2*nll_integrand(x=xj[i], 
                                           ltheta_est = ltheta_est, 
                                           ltheta_true=ltheta_true, 
                                           type=type)
    }
    for(i in seq.int(1, length(xj), 2)){
        val[i+1] <- 4*nll_integrand(x=xj[i], 
                                           ltheta_est = ltheta_est, 
                                           ltheta_true=ltheta_true, 
                                           type=type)
    }
    val[m+1] <- nll_integrand(x=upper, ltheta_est = ltheta_est, ltheta_true=ltheta_true, type=type)
    val <- val * h/3

    return(sum(val, na.rm=T))
    
}

nll_expected(log(2), log(2), "ExactSPA")
nll_expected(log(2), log(2), "SPA")
nll_expected(log(2), log(2), "StandardGauss")

# Check optimisation
opt_exact <- nlminb(log(1), nll_expected, ltheta_true=log(1), control=list(trace=1))
opt_spa <- nlminb(log(2), nll_expected, ltheta_true=log(2), type="SPA", control=list(trace=1))

# Parameters
opt_exact$par
format(opt_exact$par, scientific = F, digits=3)

exp(opt_exact$par)
format(exp(opt_exact$par), scientific = F, digits=3)
exp(opt_spa$par)

# Objective
opt_exact$objective
opt_spa$objective
nll_expected(NULL, log(2), "StandardGauss")

#Plot and see
# IG variance = 10
par(mfrow=c(1,1))
vw <- 1
chi <- psi <- 1/vw
mu <- gamma <- 0
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
x <- seq(EX-12*sqrt(VarX), EX+12*sqrt(VarX), length.out=300)
y <- exp(-sapply(x, loglikNIG, param=par))
plot(x,y, main=TeX("$\\theta = 10$"))


# For different theta = varw

# Show that the situation above happens for multiple simulations of VarW
#theta.grid <- seq(0.1,20, length.out = 30)
theta.grid <- c(0.1, 0.2, 0.5, 1)
nll.exact <- nll.spa <- nll.gauss <- numeric(length(theta.grid))
nll.exact.theta <- nll.spa.theta <- nll.exact
for(i in 1:length(theta.grid)){
    cat("iter: ", i, "\n")
    
    # Optimisation
    theta <- theta.grid[i]
    opt.exact <- nlminb(log(theta), nll_expected, ltheta_true=log(theta), control=list(trace=1))
    opt.spa <- nlminb(log(theta), nll_expected, ltheta_true=log(theta), type="SPA", control=list(trace=1))
    
    # Store solution objective
    nll.exact[i] <- opt.exact$objective
    nll.spa[i] <- opt.spa$objective
    nll.gauss[i] <- nll_expected(NULL, log(theta), "StandardGauss")
    
    # Store parameters
    nll.exact.theta[i] <- exp(opt.exact$par)
    nll.spa.theta[i] <- exp(opt.spa$par)
}


# Create a table
theta.grid
nll.exact.theta
format(nll.exact.theta, scientific = F, digits=5)
format(nll.exact, scientific = F, digits=5)
format(nll.spa.theta, scientific = T, digits=3)
format(nll.spa, scientific = F, digits=4)
format(nll.gauss, scientific = F, digits=4)


#idea: add likelihood values


par(mar=c(5, 4, 4, 2) + 0.1)

library(latex2exp)
# Plot results parameters
if(FALSE){
    setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
    pdf("spa_to_gauss_nig_par.pdf", width=8, height=5)
}
plot(theta.grid, nll.exact.theta, ylim=range(nll.exact.theta, nll.spa.theta), 
     xlab=TeX("$\\theta_0$"), ylab="",#ylab=TeX("$\\hat{\\theta}(\\theta_0)$"), 
     type="l", col=1, pch=1, lwd=4, lty=1)
mtext(2, text=TeX("$\\hat{\\theta}(\\theta_0)$"), line=2.5)
lines(theta.grid, nll.spa.theta, type="l", col=3, pch=NA, lwd=2, lty=1)
legend("topleft", c("Saddlepoint adjusted IFT", "Saddlepoint approximation"),
       lty=c(1,1), col=c(1,3), pch=c(NA,NA), lwd=c(4,2))
if(FALSE){
    dev.off()
}

# Plot results likelihood
if(FALSE){
    setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
    pdf("spa_to_gauss_nig2.pdf", width=8, height=5)
}
plot(theta.grid, nll.exact, ylim=range(nll.exact, nll.spa, nll.gauss), 
     #main="NIG likelihood VS IG variance",
     xlab="Inverse Gaussian variance", ylab="Negative log-likelihood", 
     type="b", col=1, pch=1, lwd=3, lty=1)
points(theta.grid, nll.spa, type="l", col=2, pch=2, lwd=3, lty=2)
points(theta.grid, nll.gauss, type="p", col=3, pch=3, lwd=3, lty=3)
legend("bottomleft", c("Saddlepoint adjusted IFT", "Saddlepoint approximation", "Standard Gaussian"),
       lty=c(1,2,NA), col=1:3, pch=c(1,NA,3), lwd=c(3,3,3))
if(FALSE){
    dev.off()
}



