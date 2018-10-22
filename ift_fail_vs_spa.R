# Plotting log spa and log ift
# Illustrating how direct ift fails
# Berent Lunde
setwd("C:/Users/Berent/Projects/it-ift/implementation v5/ExactSPA")

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
plot(x,y2, type="l", lty=2, lwd=3, ylab="Log-density")
lines(x,y1, type="l", col="red", lwd=2)
lines(x, rep(log(1e-14),length(x)), lty=3, lwd=2, col="blue")
legend("bottom", c("Exact", "Direct IFT", "Machine precision"), lwd=c(3,2,2), 
       col=c("black","red","blue"), lty=c(2,1,3))

setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("ift_direct_fail.pdf", width=7, height=4+1/3)
plot(x,y2, type="l", lty=2, lwd=3, ylab="Log-density")
lines(x,y1, type="l", col="red", lwd=2)
lines(x, rep(log(1e-14),length(x)), lty=3, lwd=2, col="blue")
legend("bottom", c("Exact", "Direct IFT", "Machine precision"), lwd=c(3,2,2), 
       col=c("black","red","blue"), lty=c(2,1,3))
dev.off()
