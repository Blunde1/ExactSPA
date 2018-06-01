

nig_cgf <- function(s, par){
    chi <- exp(par[1])
    psi <- exp(par[2])
    mu = par[3]
    gamma <- par[4]
    return(s*mu + sqrt(chi*psi) - sqrt(chi*psi - chi*s*s - 2.0*s*chi*gamma))
}
nig_dcgf <- function(s,par){
    chi <- exp(par[1])
    psi <- exp(par[2])
    mu = par[3]
    gamma <- par[4]
    mu - (sqrt(chi)*(-2.0*(gamma + s))) / (2.0*sqrt(psi-2.0*gamma*s-s*s))
}
nig_ddcgf <- function(s,par){
    chi <- exp(par[1])
    psi <- exp(par[2])
    mu = par[3]
    gamma <- par[4]
    as.matrix(sqrt(chi)*(psi+gamma*gamma)/( (psi-s*(2.0*gamma+s))^1.5 ) )
}
sp_fun <- function(x,par){
    nlminb(0,
           function(s,par){nig_cgf(s,par) - s*x},
           function(s,par){nig_dcgf(s,par) - x},
           nig_ddcgf, par=par#, 
           #control=list(trace=1)
           )$par
}
sp_fun(x.val, par)

r_nig_lspa <- function(x, par){
    sp <- sp_fun(x,par)
    exp(nig_cgf(sp,par) - sp*x - 0.5*log(2*pi+nig_ddcgf(sp,par)))
}
r_nig_lspa(x.val,par)

x_sim <- sort(rNIG(n, c(exp(par[1]), exp(par[2]), par[3], par[4]), seed = 12324))
hist(x_sim, freq=F)
