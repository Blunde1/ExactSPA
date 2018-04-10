# Normal inverse Gaussian - profile likelihood
library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")
pnll_fun_nig <- function(par, X, map=NULL, type="ExactSPA"){
    par <- c(par,map)
    lchi = par["lchi"]
    lpsi = par["lpsi"]
    mu = par["mu"]
    gamma = par["gamma"]
    if(type=="ExactSPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll
    }
    return(nll)
}

# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n=200
x <-rNIG(n, c(chi, psi, mu, gamma), seed=123)

map <- c(lpsi=log(100))
par2 <- par[!names(par)%in%names(map)]
opt2 <- nlminb(par2, pnll_fun_nig, map=map, X=x, control=list(trace=1), type="ExactSPA")

lpsi2 <- seq(4,6, length.out=30)
opt.nig.profile <- opt.nig.spa.profile <- list()
for(i in 1:length(lpsi2)){
    cat("iter:",i,"\n")
    map <- c(lpsi=lpsi2[i])
    opt.nig.profile[[i]] <- nlminb(par2, pnll_fun_nig, map=map, X=x, control=list(trace=1))
    opt.nig.spa.profile[[i]] <- nlminb(par2, pnll_fun_nig, map=map, X=x, type="SPA", control=list(trace=1))
    save(opt.nig.profile, file="opt.nig.profile.RData")
    save(opt.nig.spa.profile, file="opt.nig.spa.profile.RData")
}
load("opt.nig.profile.RData")
load("opt.nig.spa.profile.RData")
nll.val <- sapply(1:length(opt.nig.profile), function(i){opt.nig.profile[[i]]$objective})
nll.val.spa <- sapply(1:length(opt.nig.spa.profile), function(i){opt.nig.spa.profile[[i]]$objective})

plot(lpsi2, -nll.val, type="b",ylim=c(min(c(-nll.val,-nll.val.spa)),
                                     max(c(-nll.val,-nll.val.spa))))
lines(lpsi2, -nll.val.spa, col="red", type="b")

