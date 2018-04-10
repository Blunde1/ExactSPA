# Non-central chi - profile nll
library(ExactSPA)

pnll_fun_nchisq <- function(par, X, map=NULL, type="ExactSPA"){
    par <- c(par,map)
    ldf <- par["ldf"]
    lncp <- par["lncp"]
    if(type=="ExactSPA"){
        nll <- nll_nchisq(X, ldf, lncp, 12, 64, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nchisq(X, ldf, lncp, 12, 64, 1)$nll
    }
    return(nll)
}

set.seed(1234)
df = 100
ncp = 40
n = 500
x <- sort(rchisq(n, df, ncp))

par <- c(ldf=log(df), lncp = log(ncp))
map <- c(ldf = log(df))
par2 <- par[!names(par)%in%names(map)]
opt2 <- nlminb(par2, pnll_fun_nchisq, map=map, X=x, control=list(trace=1), type="SPA")

ldf2 <- seq(2,6, length.out=40)
opt.nchisq.profile <- opt.nchisq.spa.profile <- list()
for(i in 1:length(ldf2)){
    cat("iter:",i,"\n")
    map <- c(ldf=ldf2[i])
    opt.nchisq.profile[[i]] <- nlminb(par2, pnll_fun_nchisq, map=map, X=x, control=list(trace=1))
    opt.nchisq.spa.profile[[i]] <- nlminb(par2, pnll_fun_nchisq, map=map, X=x, type="SPA", control=list(trace=1))
    save(opt.nchisq.profile, file="opt.nchisq.profile.RData")
    save(opt.nchisq.spa.profile, file="opt.nchisq.spa.profile.RData")
}

load("opt.nchisq.profile.RData")
load("opt.nchisq.spa.profile.RData")
nll.val <- sapply(1:length(opt.nchisq.profile), function(i){opt.nchisq.profile[[i]]$objective})
nll.val.spa <- sapply(1:length(opt.nchisq.spa.profile), function(i){opt.nchisq.spa.profile[[i]]$objective})

plot(ldf2[1:30], -nll.val[1:30], type="b")#,ylim=c(min(c(-nll.val,-nll.val.spa)),
#                                               max(c(-nll.val,-nll.val.spa))))
lines(ldf2[1:30], -nll.val.spa[1:30], col="red", type="b")


