# NIG High dens and Low dens Re CF
# Berent Lunde
# 30.08.2018


setwd("C:/Users/Berent/Projects/it-ift/implementation v5/ExactSPA")

library(Rcpp)
Rcpp::sourceCpp('nig_illustrations.cpp')
library(ExactSPA)
ExactSPA::hello()

# NIG RE CF STANDEXPTILT ####
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2
hDensX <- EX # EX
lDensX <- EX - 8*sqrt(VarX) #EX - 5 SDX

s_val <- seq(0,40, by=0.1)
s_val2 <- seq(0,40, by=0.1)

# Mean and tail together
par(mfrow=c(1,1))
par(mar=c(5,6, 3, 3))
y <- sapply(s_val, reStandCFNIG, x=lDensX, lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

library(latex2exp)
setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("nig_stand_cf3.pdf", width=7, height=4+1/3)
plot(s_val,y, type="l", lwd=2, lty=1, pch=2,
     main="",
     ylab=expression(paste(Re,"[", varphi[bar(x)(hat(tau))](s),"]")),
     xlab=expression(s))
# N(0,1)
points(s_val2,sqrt(2*pi)*dnorm(s_val2), type="l", lty = 2, col="red", lwd=3)
# mean
y2 <- sapply(s_val, reStandCFNIG, x=hDensX, lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
points(s_val, y2, type="l", lty=3, lwd=3, col="blue")
# legend("topright", 
#        legend=c(expression(paste("High density ", Re,"[", varphi[bar(x)(tilde(tau))](s),"]")),
#                 expression(paste("Low density ", Re,"[", varphi[bar(x)(tilde(tau))](s),"]")),
#                 "Standard normal"),
#        col=c("blue","black","red"), pch=c(NA,NA,NA), lty=c(3,1,2), lwd=c(3,2,3)
# )
legend("topright", 
       legend=c(TeX("High density $x_0$"), TeX("Low density $x_0$"), TeX("$exp(-0.5s^2)$")),
       col=c("blue","black","red"), pch=c(NA,NA,NA), lty=c(3,1,2), lwd=c(3,2,3)
)
dev.off()
