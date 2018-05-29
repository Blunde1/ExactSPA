library(ExactSPA)

par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015)) #jump intensity usually 100
r = par[1]
lsigma = par[2]
llambda = par[3]
mu = par[4]
lnu = par[5]
x0 <- 10
time <- 4
N <- time*250
dt <- time/N
X<-rMJD(N,time,x0=x0,
        r=par["r"], sigma=exp(par["lsigma"]), jump_intensity = exp(par["llambda"]),
        mu=par["mu"], nu=exp(par["lnu"]),
        seed=123)
Xt <- log(X)
range(Xt)

xt <- seq(2.1,2.5,length.out=100)
x0 <- log(10)
log_ift <- log_espa <- rep(0, length(xt))

for(i in 1:length(xt)){
    cat("iter: ",i, "\n")
    #log_ift[i] <- - nll_fun_mjd(par,X=c(x0,xt[i]) , dt=dt, type = "Simpson")
    log_ift[i] <- -nll_mjd(c(x0,xt[i]), dt, r, lsigma, llambda, mu, lnu, 700, 64, 3)$nll
}
plot(xt, log_ift, type="b")
for(i in 1:length(xt)){
    cat("iter: ",i, "\n")
    #log_espa[i] <- - nll_fun_mjd(par,X=c(x0,xt[i]) , dt=dt, type = "ExactSPA")
    log_espa[i] <- - nll_mjd(c(x0,xt[i]), dt, r, lsigma, llambda, mu, lnu, 15, 128, 2)$nll
}
lines(xt, log_espa, type="l",col="red")
nll <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 700, 64, 3)$nll
nll <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64, 2)$nll

# Standard Normal distribution
cf <- function(s){
    exp(- s^2 / 2)
}
ift <- function(s, x){Re(cf(s)*exp(-1i*s*x))/(2*pi)}
fx.ift <- function(x){
    #s <- seq(-6,6,length.out=256)
    #abs(sum(sapply(s, function(s){Re(cf(s)*exp(-1i*s*x))/(2*pi)}))*(s[2]-s[1]))
    abs(integrate(ift, -Inf, Inf, x=x)$value) #not guaranteed to stay positive
}
fx.ift(0)
dnorm(0)
cgf <- function(s){
    s^2 / 2
}
cgf.1 <- function(s){
    s
}
cgf.2 <- function(s){
    1
}
sp <- function(x){
    x
}
fx.spa <- function(x){
    s.hat <- sp(x)
    exp(cgf(s.hat)- s.hat*x)/sqrt(2*pi*cgf.2(s.hat))
}
fx.ift(-10)
fx.spa(-10)
dnorm(-10)

x <- seq(-15,15,length.out=100)
y1 <- log(sapply(x, fx.ift))
y2 <- log(sapply(x, fx.spa))
plot(x,y2, type="l")
lines(x,y1, type="l", col="red")
legend("topleft", c("SPA adjusted", "Direct integration"), lwd=c(1,1), col=c("black","red"))
