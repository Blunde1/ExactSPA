

# Plot tweedie cf
phi.tweeide <- function(s,mu=2, a=1, p=2.5){
    exp(
        mu^(2-p)/(a*(2-p))*
            (1+1i*s*a*(1-p)*mu^(p-1))^((2-p)/(1-p))-1
    )
}

phi.tweeide(-1)

integrand <- function(s,x,mu, a, p){
    Re(phi.tweeide(s=s, mu=mu, a=a, p=p)*exp(-1i*s*x))/(pi)
}
integrand(s=-1,2.4, 2, 1, 2.5)

dens.tweedie <- function(x,a=1,mu=5,p=2.5){
    integrate(integrand, lower=0, upper=Inf, x=x, mu=mu, a=a, p=p)$value
}
dens.tweedie(2.4, 1, 2, 2.5)
dens.tweedie(2.4, p=3.2)

x <- seq(0,10, length.out=100)
fx <- sapply(x, dens.tweedie)
plot(x, fx)
fx.2 <- sapply(x, dens.tweedie, p=3.5)
plot(x, fx.2, col="red")
sum(fx.2)*(x[2]-x[1])
