# MJD OPT and SD
# Berent Lunde
# 22.09.2018

setwd("C:/Users/Berent/Projects/it-ift/implementation v5")

# PACKAGES
library(ExactSPA)
ExactSPA::hello()
library(Quandl)

# Quandl API
qapi <- "teWRSvY5Kygi_etMye4C"
Quandl.api_key(qapi)

# DATA
start_date <- "2000-01-01"; end_training <- "2018-01-01";
DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
DJIA <- DJIA[rev(rownames(DJIA)),]
plot(DJIA,type="l")
Xt=log(DJIA$Value)

# PARAMETERS
par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015))
dt <- 1/252

# WRAPPERS
nll_mjd_wrap <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]

    n <- length(X)
    nll_contributions <- numeric(n-1)
    for(i in 2:n){
        nll_contributions[i-1] <- nll_mjd(X[(i-1):i], dt, r, lsigma, llambda, mu, lnu, 16, 128, 2)$nll
    }

    # return individual likelihood contributions
    return(sum(nll_contributions))
}

dnll_mjd_wrap <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]

    n <- length(X)
    grad <- rep(0,5)
    for(i in 2:n){
        res <- nll_mjd(X[(i-1):i], dt, r, lsigma, llambda, mu, lnu, 16, 128, 2)
        grad <- grad +
            c(
                res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
            )
    }


    # return individual likelihood contributions
    return(grad)
}

# Optimisation
opt <- nlminb(par, nll_mjd_wrap, dnll_mjd_wrap, X=Xt, dt=dt, control=list(trace=1))
opt
cat(
    format(c(opt$par[1],
             exp(opt$par[2:3]),
             opt$par[4],
             exp(opt$par[5])),
           scientific = F, digits=4),
    sep = " & "
)
par_orig <- c(opt$par[1],
          exp(opt$par[2:3]),
          opt$par[4],
          exp(opt$par[5]))

# SD
# build Hessian

dfun_trans <- function(par_orig, Xt, dt){
    # takes in actual parameters
    par_trans <- c(par_orig[1], log(par_orig[2]), log(par_orig[3]), par_orig[4], log(par_orig[5]))
    grad <- dnll_mjd_wrap(par_trans, Xt, dt)
    dtrans <- c(1, 1/par_orig[2], 1/par_orig[3], 1, 1/par_orig[5])
    return(dtrans * grad)
}
dfun_trans(par_orig, Xt, dt)

fun <- dfun_trans
x = par_orig

jacobian <- function(fun, x, dt=dt, Xt=Xt){
    n <- length(x)
    J <- matrix(nrow=n, ncol=n)
    fx <- fun(par=x, X=Xt, dt=dt)
    eps <- 1e-8
    xperturb = x
    for(i in 1:n){
        xperturb[i] <- xperturb[i] + eps
        J[,i] <- (fun(par=xperturb, X=Xt, dt=dt) - fx) / eps
        xperturb[i] <- x[i]
    }
    return(J)
}
Hess <- J
Hess <- jacobian(fun=dfun_trans, x=opt$par)
solve(Hess)
SD <- sqrt(diag(solve(Hess)))


tab_vec <- c()
for(i in 1:length(par)){
    tab_vec <- paste0(tab_vec,
                      cat(format(par_orig[i], scientific = F, digits=3)),
                      cat(" (", format(SD[i], scientific = F, digits=3), ") & ", sep="")
                      , sep ="")
}
tab_vec


# TRUNCATED
nll_mjd_trun <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    nll <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 700, 64, 4)$nll
    return(nll)
}
dnll_mjd_trun <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    res <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 700, 64, 4)
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )
}


nll_mjd_trun <- function(par, X, dt){

    r = par[1]
    lsigma = par[2]
    sigma <- exp(lsigma)
    llambda = par[3]
    lambda <- exp(llambda)
    mu = par[4]
    lnu = par[5]
    nu <- exp(lnu)

    k <- exp(mu + 0.5*nu^2) - 1

    # Calculate weights until P(jump) < 1e-12
    # But maximum 5 jumps
    n <- length(X)
    maxjumps <- 20
    nll <- 0
    for(j in 2:n){
        prob <- 0
        for(i in 0:maxjumps){
            w <- dpois(i, lambda*dt)
            if(w < 1e-14) break
            #cat(w, " - ", prob, "\n")
            prob <- prob + w * dnorm(X[j],
                                     X[j-1] + (r - sigma^2/2 - lambda*k)*dt + i*mu,
                                     sqrt(dt*sigma^2 + i*nu^2))
        }
        nll <- nll - log(prob)
    }
    return(nll)
}
nll_mjd_trun(opt$par, Xt, dt)
nll_mjd_wrap(opt$par, Xt, dt)
opt$objective

#grad
dnll_mjd_trun <- function(par, X, dt){
    m <- length(par)
    parpartub <- par
    eps <- 1e-6
    G <- rep(0,m)
    fx <- nll_mjd_trun(par, X, dt)

    for(i in 1:m){
        parpartub[i] <- parpartub[i] + eps
        fxeps <- nll_mjd_trun(parpartub, X, dt)
        G[i] <- (fxeps - fx) / eps
        parpartub[i] <- par[i]
    }

    return(G)
}
dnll_mjd_wrap(opt$par, Xt, dt)
dnll_mjd_trun(opt$par, Xt, dt)

