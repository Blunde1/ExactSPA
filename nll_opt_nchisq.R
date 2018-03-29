library(ExactSPA)

df = 100
ncp = 40
n = 2000
x <- sort(rchisq(n, df, ncp))
hist(x, freq=F)
curve(dchisq(x, df, ncp), add=TRUE)

y <- rep(0, length(x))
for(i in 1:length(x)){
    y[i] <- exp(-nll_nchisq(as.vector(x[i]),log(df),log(ncp), 12.0, 64)$nll)
}
lines(x,y,col="blue")

-sum(dchisq(x, df, ncp, log=TRUE))
nll_nchisq(x, log(df), log(ncp), 20.0, 128)


par <- c(ldf=log(df), lncp = log(ncp))
opt <- nlminb(par, nll_fun_nchisq, nll_grad_nchisq, X=x, control=list(trace=1))
opt
opt$par
exp(opt$par)
