# R script for simulating random variables

# Merton Jump Diffusion ####

mjd_increment <- function(dt, x0, r, sigma, jump_intensity, mu, nu){
    lx0 <- log(x0)
    number_of_jumps <- rpois(1,dt*jump_intensity)
    compounded_poisson <- sum(rnorm(number_of_jumps, mu, nu))
    k = exp(mu+0.5*nu^2)-1
    gbm_increment <- (r-0.5*sigma^2 - jump_intensity*k)*dt + sigma*rnorm(1,0,sqrt(dt))
    exp(lx0 + gbm_increment + compounded_poisson)
}

rMJD <- function(N, Time, x0, r, sigma, jump_intensity, mu, nu, seed){
    set.seed(seed)
    x <- c(x0)
    dt<-Time/N
    for (i in 2:(N+1)) {
        x[i]  <-  mjd_increment(dt, x[i-1], r, sigma, jump_intensity, mu, nu)
    }
    return(x);
}
#test<-rMJD(N=20*250, Time=20, x0=1, 0.1, 0.1, 100, -0.001, 0.015, seed=123)
#plot(test,type="l")

# Non-central chi-square ####
### rchisq(n,df,ncp)

# Normal Inverse Gaussian ####
require(SuppDists)
rNIG=function(n,param, seed){
    set.seed(seed)
    # param[1] : chi
    # param[2] : psi
    # param[3] : mu
    # param[4] : gamma
    #Alternative parametrization used in SuppDists
    nu<-sqrt(param[1]/param[2])
    #Simulate n inverse gaussian observations
    IGsim<-rinvGauss(n,nu,param[1])
    #Constructing n NIG observations based upon IGsim
    NIGsim<-param[3]+param[4]*IGsim+sqrt(IGsim)*rnorm(n)
    #returning the simulated NIG observations
    return(NIGsim)
}
#x <- sort(rNIG(n, c(chi, psi, mu, gamma)))