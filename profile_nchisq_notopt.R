# Script for plotting profile likelihood (not in optimum) for NCHISQ RV
# Berent Å. S. Lunde
# 18.05.2018

library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")

# Poisson central chi mix
nll_nchisq_r <- function(par, data){
    -sum(log(dchisq(data, exp(par[1]), exp(par[2]), log=F)))
}

# Simulate datasets
df = 100
ncp = 40
par <- c(ldf=log(df), lncp = log(ncp))
n <- 500

# NEW PAR
df = 1
ncp = 0.5
par <- c(ldf=log(df), lncp = log(ncp))
n = 1

set.seed(123)
x.nchisq <-rchisq(n, df, ncp)

nll_nchisq_r(par, x.nchisq)
nll_fun_nchisq(par, x.nchisq, "ExactSPA")
nll_fun_nchisq(par, x.nchisq, "SPA")


par["ldf"] <- log(10)
ncp.val <- seq(1,10000, length.out=100)
nll_r <- numeric(length(ncp.val))
for(i in 1:length(ncp.val)){
    par["lncp"] <- log(ncp.val[i])
    nll_r[i] <- nll_nchisq_r(par, x.nchisq)
}
plot(ncp.val, nll_r)

# Can espa handle it?
par["lncp"] <- log(4000)
nll_fun_nchisq(par, x.nchisq) # ok
par["lncp"] <- log(6000)
nll_fun_nchisq(par, x.nchisq) # great

ncp.val2 <- seq(1, 5000, length.out=50)
nll_espa <- nll_bessel <- numeric(length(ncp.val2))
for(i in 1:length(ncp.val2)){
    cat("iter: ",i)
    par["lncp"] = log(ncp.val2[i])
    nll_espa[i] <- nll_fun_nchisq(par, x.nchisq)
    nll_bessel[i] <- nll_nchisq_r(par, x.nchisq)
}
plot(ncp.val2, nll_espa, type="o",pch=3, lwd=2,
     ylab="Log-likelihood", xlab="Noncentrality")
lines(ncp.val2, nll_bessel, type="b", col="red", pch=2, lwd=2)
legend("topleft", c("Exact SPA", "Poisson - Chi-squared mixture"), col=c("black","red"), pch=c(3,2), lwd=c(2,2))


