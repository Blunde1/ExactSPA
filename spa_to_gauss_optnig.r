# illustrate spa --> Gaussian under optimisation
# illustrate with NIG
# Berent Lunde
# 23.08.2018

# # CONDITIONS # #
# EW = 1 --> sqrt(chi) = sqrt(psi)
# mu = gamma = 0
# VarW = sqrt(chi) / psi^(3/2) = 1/ psi = 1/chi

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

# Likelihood wrappers
nll_nig_locked <- function(lVarW, X, type="ExactSPA"){
    VarW <- exp(lVarW)
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

# Locked parameters
mu <- gamma <- 0

# Simulate nig (under conditions)
n = 100 #n=200
vw <- 2
chi <- psi <- 1/vw
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
set.seed(4321)
x.nig <- rNIG(n,  c(chi, psi, mu, gamma), seed=123)

# Check valules of likelihoods
nll_fun_nig(par, X=x.nig)
loglikNIG(par, x.nig)
nll_fun_nig(par, X=x.nig, type="SPA")
nll_fun_nig(par, X=x.nig, type="reSPA")
nll_nig_locked(log(vw), X=x.nig)
# OK

# Estimate parameters (psi = chi)
opt.exact <- nlminb(log(vw), nll_nig_locked, X=x.nig, type="ExactSPA", control=list(trace=1))
exp(opt.exact$par)
opt.spa <- nlminb(log(vw), nll_nig_locked, X=x.nig, type="SPA", control=list(trace=1))
exp(opt.spa$par)
opt.respa <- nlminb(log(vw), nll_nig_locked, X=x.nig, type="reSPA", control=list(trace=1))
exp(opt.respa$par)
opt.gauss <- nlminb(c(0,0), nll_normal, X=x.nig, control=list(trace=1))

nll_normal(opt.gauss$par, x.nig)
nll_nig_locked(opt.spa$par, X=x.nig, type="SPA")
nll_nig_locked(opt.exact$par, X=x.nig, type="ExactSPA")

# Under conditions, N(0,1) is a special case, when VarW-->0
# Indeed, spa seeks this solution
# --> Do not need the reference gauss opt
# --> Only N(0,1) nll


# Show that the situation above happens for multiple simulations of VarW
VarW.grid <- seq(0.1,4, length.out = 30)
nll.exact <- nll.spa <- nll.gauss <- numeric(length(VarW.grid))
set.seed(52321)
for(i in 1:length(VarW.grid)){
    cat("iter: ", i, "\n")
    
    # simulate
    vw <- VarW.grid[i]
    chi <- psi <- 1/vw
    x.nig <- rNIG(n,  c(chi, psi, mu, gamma), seed = 123)
    
    # Estimate
    opt.exact <- nlminb(log(vw), nll_nig_locked, X=x.nig, type="ExactSPA", control=list(trace=1))
    opt.spa <- nlminb(log(vw), nll_nig_locked, X=x.nig, type="SPA", control=list(trace=1))
    
    # Store solutions
    nll.exact[i] <- opt.exact$objective
    nll.spa[i] <- opt.spa$objective
    nll.gauss[i] <- nll_normal(c(0,0), x.nig)
}

# Plot results
if(FALSE){
    setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
    pdf("spa_to_gauss_nig.pdf", width=8, height=5)
}
plot(VarW.grid, nll.exact, ylim=range(nll.exact, nll.spa, nll.gauss), 
     main="NIG likelihood VS IG variance",
     xlab="Inverse Gaussian variance", ylab="Negative log-likelihood", 
     type="b", col=1, pch=1, lwd=3, lty=1)
points(VarW.grid, nll.spa, type="l", col=2, pch=2, lwd=3, lty=2)
points(VarW.grid, nll.gauss, type="p", col=3, pch=3, lwd=3, lty=3)
legend("bottomleft", c("Saddlepoint adjusted IFT", "Saddlepoint approximation", "Standard Gaussian"),
       lty=c(1,2,NA), col=1:3, pch=c(1,NA,3), lwd=c(3,3,3))
if(FALSE){
    dev.off()
}

#