

x=-0.0
p.0 <- function(x){
    u <- seq(-25,25,0.005)
    y <- sapply(u, function(u){Re(phi.stand(u, x))/(2*pi)})
    sum(y)*(u[2]-u[1])
}
p.0(x)

# GH:
install.packages("gaussquad")
library("gaussquad")
order <- 50
gh <- ghermite.h.quadrature.rules(n=order, mu=0)[[order]]
#gh
sum( gh$w * Re( sapply(gh$x, function(u)phi.stand(sqrt(2)*u,x) * exp(u^2) )) ) / (sqrt(2)*pi)

order <- 100
gl <- glaguerre.quadrature.rules(order, alpha=0)[[order]]
#gl
sum( gl$w * Re( sapply( gl$x, function(u) phi.stand( sqrt(2*u),x )*exp(u)/sqrt(2*u)  ))) / pi


?gauss.hermite

integrand <- function(u, x) Re(phi.stand(u,x))/(2*pi)

library(microbenchmark)
almost.exact <- integrate(integrand, lower=-Inf, upper=Inf, x=x)

composite.simpson <- function(f, a, b, n) {
    if (is.function(f) == FALSE) {
        stop('f must be a function with one parameter (variable)')
    }
    
    h <- (b - a) / n
    
    xj <- seq.int(a, b, length.out = n + 1)
    xj <- xj[-1]
    xj <- xj[-length(xj)]
    
    approx <- (h / 3) * (f(a) + 2 * sum(f(xj[seq.int(2, length(xj), 2)])) + 4 * sum(f(xj[seq.int(1, length(xj), 2)])) + f(b))
    
    return(approx)
    
}
f.simpson <- function(u) Re(phi.stand(u,x=x))/(2*pi)
f.simpson(0)
composite.simpson(f.simpson, -20,20,64)
f.simpson.half <- function(u) Re( phi.stand(u,x=x)) / pi
f.simpson.half(0)
composite.simpson(f.simpson.half, 0, 20, 32)

# check order for gh
max.order <- 200
step.length <- 10
int.eval.simpson <- int.eval.gh <- int.eval.gl <- rep(0, max.order/step.length)
for(i in 1:(max.order/step.length)){
    order <- i*step.length
    gh <- ghermite.h.quadrature.rules(n=order, mu=0)[[order]]
    gl <- glaguerre.quadrature.rules(order, alpha=0)[[order]]
    int.eval.gh[i] <- sum( gh$w * Re( sapply(gh$x, function(u)phi.stand(sqrt(2)*u,x) * exp(u^2) )) ) / (sqrt(2)*pi)
    int.eval.gl[i] <- sum( gl$w * Re( sapply( gl$x, function(u) phi.stand( sqrt(2*u),x )*exp(u)/sqrt(2*u)  ))) / pi
    int.eval.simpson[i] <- composite.simpson(f.simpson, -20, 20, i*step.length)
}
plot(log(int.eval.gh), ylab="Log Probability", xlab="Order", ylim=c(-1.2,-0.85))
points(log(int.eval.gl),col="red")
points(log(int.eval.simpson), col="purple")
abline(h=log(almost.exact$value),col="blue")
legend("bottomright", c("Gauss-Hermite", "Gauss-Laguerre","Simpson","kronrod"), col=c("black","red","green","blue"), 
       lty=c(1,1,1,1))

# Integration on direct characteristic function
x=-0.1
integrand.originalcf <- function(s,x){
    Re(exp(mjd.cgf(1i*s) -1i*s*x )) / (pi)
}
int <- integrate(integrand.originalcf, -Inf, Inf, x=x)
int$subdivisions

s <- seq(-1000,1000,1)
cf.eval <- sapply(s, integrand.originalcf, x=x)
plot(s, cf.eval, type="l")
