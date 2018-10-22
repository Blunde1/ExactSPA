setwd("C:/Users/Berent/Projects/it-ift/implementation v5/ExactSPA")

K.tweedie <- function(s, x, xi, mu, phi){
    mu^(2-xi)/(phi*(2-xi)) * ((1-s*phi*(xi-1)*mu^(xi-1))^((2-xi)/(1-xi)) - 1)
}

dK.tweedie <- function(s, x, xi, mu, phi){
    mu * (1-s*phi*(xi-1)*mu^(xi-1))^(1/(1-xi)) 
}

ddK.tweedie <- function(s, x, xi, mu, phi){
    phi*mu^xi * (1-s*phi*(xi-1)*mu^(xi-1))^(xi/(1-xi)) 
}

tweedie.spa <- function(x, xi, mu, phi){
    # find saddlepoint
    s <- nlminb(0, function(s){K.tweedie(s,x,xi,mu,phi) - s*x}, function(s){dK.tweedie(s,x,xi,mu,phi) - x})$par
    lspa <- K.tweedie(s, x, xi, mu, phi) - s*x - 0.5*(log(2*pi) + log(ddK.tweedie(s,x,xi,mu,phi)))
}


p <- 5
mu <- 1
phi <- 1
y <- seq(0, 10, length=100)
fy <- dtweedie( y=y, power=p, mu=mu, phi=phi)
plot(y, fy, type="l")
# Compare to the saddlepoint density
f.saddle <- dtweedie.saddle( y=y, power=p, mu=mu, phi=phi)
lines( y, f.saddle, col=2 )

f.spa <- sapply(y[-1], tweedie.spa, xi=p, mu=mu, phi=phi)
lines(y[-1], exp(f.spa),col=3)
