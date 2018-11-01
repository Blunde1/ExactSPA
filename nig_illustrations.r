# NIG Plots

setwd("C:/Users/Berent/Projects/it-ift/implementation v5/ExactSPA")

library(Rcpp)
Rcpp::sourceCpp('nig_illustrations.cpp')
library(ExactSPA)

# # # DEPRECATED - SEE nig_HighLowDens_ReCF.r # # # 

# NIG RE CF STANDEXPTILT ####
n <- 10000
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)


x <- rNIG(n, c(chi, psi, mu, gamma), seed = 123)
m <- mean(x)
sd <- sd(x)
hist(x)

s_val <- seq(0,40, by=0.1)
s_val2 <- seq(0,40, by=0.1)

# Mean and tail together
par(mfrow=c(1,1))
par(mar=c(5,6, 3, 3))
y <- sapply(s_val, reStandCFNIG, x=min(x), lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_stand_cf.pdf", width=7, height=4+1/3)
plot(s_val,y, type="l", lwd=2, lty=1, pch=2,
     main="",
     ylab="Value",#expression(paste(Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
     xlab=expression(s))
# N(0,1)
points(s_val2,sqrt(2*pi)*dnorm(s_val2), type="l", lty = 2, col="red", lwd=3)
# mean
y2 <- sapply(s_val, reStandCFNIG, x=mean(x), lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
points(s_val, y2, type="l", lty=3, lwd=3, col="blue")
legend("topright", 
       legend=c(expression(paste("High density ", Re,"[", varphi[bar(x)(tilde(tau))](s),"]")),
                expression(paste("Low density ", Re,"[", varphi[bar(x)(tilde(tau))](s),"]")),
                "Standard normal"),
       col=c("blue","black","red"), pch=c(NA,NA,NA), lty=c(3,1,2), lwd=c(3,2,3)
)
dev.off()






# tail
par(mfrow=c(1,1))
par(mar=c(5,6, 3, 3))
y <- sapply(s_val, reStandCFNIG, x=min(x), lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_stand_cf_tail_x.pdf", width=7, height=4+1/3)
plot(s_val,y, type="l", lwd=2, lty=1, pch=2,
     main="",
     ylab="Value",#expression(paste(Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
     xlab=expression(s))
# N(0,1)
points(s_val2,dnorm(s_val2), type="l", lty = 2, col="red", lwd=3)
# mean
y2 <- sapply(s_val, reStandCFNIG, x=mean(x), lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
points(s_val, y2, type="l", lty=3, lwd=2, col="blue")
legend("topright", 
       legend=c(expression(paste("High density ", Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
                expression(paste("Low density ", Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
                "Standard normal"),
       col=c("black","blue","red"), pch=c(NA,NA,NA), lty=c(3,1,2), lwd=c(2,2,3)
       )
dev.off()

# mean
pdf("nig_stand_cf_mean_x.pdf", width=7, height=4+1/3)
s_val <- seq(-10,10, by=0.1)
s_val2 <- seq(-10,10, by=0.5)
y <- sapply(s_val, reStandCFNIG, x=mean(x), lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
plot(s_val,y, type="l", lwd=2, lty=4,
     main="", ylab="Value", xlab=expression(s))
points(s_val2,dnorm(s_val2), pch=18, col="red", lwd=2)
legend("topright", 
       legend=c(expression(paste(Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
                "Standard normal"),
       col=c("black","red"), pch=c(NA,18), lty=c(4,NA), lwd=c(2,2)
)
dev.off()


# NIG SPA VS f VS p(0) vs sqrt(2*pi) ####
par(mfrow=c(1,1))
range(x)
x_min <- -0.5#min(x)
x_max <- 0.5#max(x)
x <- seq(x_min, x_max, length.out=1000)

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
        sum(
            0.5*log(param[1]) +
                log(param[2] + param[4]^2) +
                sqrt(param[1]*param[2]) -
                log(pi) +
                log(besselK(y, -1, expon.scaled = TRUE)) - y +
                (data - param[3])*param[4] -
                log(y)
        )
    # Return the functional value
    return(loglik)
}
loglikNIG(par, x)
nll_fun_nig(par, x, type = "SPA")

y <- exp(sapply(x, loglikNIG, param=par))
y2 <- exp(-sapply(x, nll_fun_nig, par=par))
yspa <- exp(-sapply(x, nll_fun_nig, par=par, type="SPA"))


# double plot
setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_exact_vs_spa.pdf", width=7, height=4+1/3)
par(mar=c(5,6, 3, 6))
plot(x,y, type="l", lwd=3, lty=1,xaxt="n", yaxt="n",
     main="", xlab=expression(x), ylab="Density")
points(x,yspa, type="l", lty=2, lwd=3, col="red")
axis(2, ylim=range(y), col="black",lwd=1,las=1)
axis(1)
legend("topleft", legend=c("Exact", "SPA"), 
       lwd=c(3,3), lty=c(1,2), col=c("black","red"))
par(new=T)
p0 <- y / (yspa * sqrt(2*pi))
plot(x,p0, type="l", lwd=2, lty=3,ylim=c(0,1),#ylim=c(min(1/sqrt(2*pi), p0), max(1/sqrt(2*pi), p0)),
     col="blue",axes=F, main="", xlab="", ylab="")
lines(x,rep(1/sqrt(2*pi),length(x)), lty=4, col="green", lwd=2)
axis(4, ylim=c(0,1))#ylim=c(min(p0), max(p0)))
mtext(4, text=expression(paste(
    p[bar(X)(hat(tau))](0)
)), line=3)
legend("topright", legend=c("Exact", "SPA"),
       col=c("blue","green"),lwd=c(2,2), lty=c(3,4))
dev.off()


# new version 29.10.2018
# two plots
setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_exact_vs_spa_2plot.pdf", width=9, height=6)
par(mar=c(5, 4, 4, 2) + 0.1) # default
par(mfrow=c(1,2))
plot(x,y, type="l", lwd=3, lty=1,xaxt="n", yaxt="n",
     main="", xlab=TeX("$x$"), ylab="Density")
points(x,yspa, type="l", lty=2, lwd=3, col="red")
axis(2, ylim=range(y), col="black",lwd=1,las=1)
axis(1)
legend("topleft", legend=c("Exact", "SPA"), 
       lwd=c(3,3), lty=c(1,2), col=c("black","red"))
p0 <- y / (yspa * sqrt(2*pi))
plot(x,p0, type="l", lwd=3, lty=1,ylim=c(min(1/sqrt(2*pi), p0), max(1/sqrt(2*pi), p0)),
     col=1,axes=T, 
     main="", 
     xlab= TeX("x"),
     ylab= TeX("$p_{\\bar{X}(\\hat{\\tau})}(0)$")
)
points(x,rep(1/sqrt(2*pi),length(x)), 
       type="l", lty=2, lwd=3, col="red")
legend("topright", legend=c("Exact", "SPA"),
       col=c("black","red"),lwd=c(3,3), lty=c(1,2))
dev.off()







pdf("nig_exact_vs_spa.pdf", width=7, height=4+1/3)
plot(x,y, type="l", lwd=2, lty=4,
     main="", xlab=expression(x), ylab="Density")
points(x,yspa, pch=18, lwd=1, col="red")
lines(x,y, type="l", lwd=2, lty=4) # Add one more to get black over red
legend("topright", legend=c("Exact", "NIG SPA"), 
       pch=c(NA,18), lwd=c(2,1), lty=c(4,NA), col=c("black","red"))
dev.off()
#lines(x,y2, lty=3, lwd=2, col="blue")

p0 <- y / (yspa * sqrt(2*pi))
pdf("nig_exp_tilt_p0.pdf", width=7, height=4+1/3)
plot(x,p0, type="l", lwd=2, lty=4,ylim=c(min(1/sqrt(2*pi), p0), max(1/sqrt(2*pi), p0)),
     main="Standardised characteristic function IFT", xlab=expression(x), ylab="Value")
lines(x,rep(1/sqrt(2*pi),length(x)), lty=3, col="red", lwd=2)
legend("topright", legend=c("Exact", "SPA"),
       col=c("black","red"),lwd=c(2,2), lty=c(4,3))
dev.off()
#abline(h=1/sqrt(2*pi), col="red", lty=2, lwd=2)


# tail likelihood
# Try to get non-gaussianity at one point and gaussianity at other, where ift and Bessel fails


# CF NIG
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
x <- rNIG(n, c(chi, psi, mu, gamma), seed = 123)

par["lpsi"] = log(5e7)
par["lchi"] = log(0.000001)
par["mu"] = -0.01
par["gamma"] <- 100


s_val <- seq(-1000,1000, length.out=1000)
y <- sapply(s_val, reCFNIG, x=mean(x), lchi=par[1], lpsi=par[2], mu=par[3], gamma=par[4])
plot(s_val, y) # independent of x


