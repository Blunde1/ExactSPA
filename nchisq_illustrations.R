# NCHISQ ILLUSTRATIONS
# Berent Lunde
# 05.05.2018

# Poisson central chi mix
nll_nchisq_r <- function(par, data){
    -sum(log(dchisq(data, exp(par[1]), exp(par[2]), log=F)))
}

# wrapper - part of lib
nll_fun_nchisq <- function(par, X, type="ExactSPA"){
    ldf = par[1]
    lncp = par[2]
    if(type=="ExactSPA"){
        nll <- nll_nchisq(X, ldf, lncp, 2000, 5000, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nchisq(X, ldf, lncp, 12, 64, 1)$nll
    }
    return(nll)
}

df = 1
ncp = 0.5
par <- c(ldf=log(df), lncp = log(ncp))
n = 1

nll_nchisq_r(par, x.nchisq)
nll_fun_nchisq(par, x.nchisq, "ExactSPA")
nll_fun_nchisq(par, x.nchisq, "SPA")

# Plot transition density
x <- sort(rchisq(1000, df=df, ncp=ncp))
range(x)
hist(x, freq = F)
f.r <- dchisq(x, df=df, ncp=ncp, log=F)
lines(x,f.r)





plot(x,f.mix)
