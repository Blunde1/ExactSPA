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

#jacobian <- function(fun, x, dt=dt, Xt=Xt){
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
#    return(J)
#}
Hess <- J
#Hess <- jacobian(fun=dfun_trans, x=opt$par)
solve(Hess)
SD <- sqrt(diag(solve(Hess)))


tab_vec <- c("")
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

    # n <- length(X)
    # grad <- rep(0,5)
    # for(i in 2:n){
    #     res <- nll_mjd(X[(i-1):i], dt, r, lsigma, llambda, mu, lnu, 16, 128, 4)
    #     grad <- grad +
    #         c(
    #             res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    #         )
    # }
    #
    #
    # # return individual likelihood contributions
    # return(grad)



    res <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 700, 64, 4)
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )
}
nll_mjd_trun(opt$par, Xt, dt)
nll_mjd_wrap(opt$par, Xt, dt)
opt$objective

# GRAD
dnll_mjd_wrap(opt$par, Xt, dt)
dnll_mjd_trun(opt$par, Xt, dt)

ind <- floor(seq(0, length(Xt), length.out = 20))
d1 <- d2 <- list()
for(i in 2:length(ind)){
    cat("iter: ", i-1, "\n")
    d1[[i-1]] <- dnll_mjd_wrap(opt$par, Xt[ind[i-1]:ind[i]], dt)
    d2[[i-1]] <- dnll_mjd_trun(opt$par, Xt[ind[i-1]:ind[i]], dt)
}
dr1 <- dlsigma1 <- dllambda1 <- dmu1 <- dlnu1 <- dr2 <- dlsigma2 <- dllambda2 <- dmu2 <- dlnu2 <- length(d1)
for(i in 2:length(ind)){
    dr1[i-1] <- d1[[i-1]][1]
    dlsigma1[i-1] <- d1[[i-1]][2]
    dllambda1[i-1] <- d1[[i-1]][3]
    dmu1[i-1] <- d1[[i-1]][4]
    dlnu1[i-1] <- d1[[i-1]][5]
    dr2[i-1] <- d2[[i-1]][1]
    dlsigma2[i-1] <- d2[[i-1]][2]
    dllambda2[i-1] <- d2[[i-1]][3]
    dmu2[i-1] <- d2[[i-1]][4]
    dlnu2[i-1] <- d2[[i-1]][5]
}
# dr
plot(ind[-1], dr1)
points(ind[-1], dr2, col=2, pch=2)
abline(h= sum(dr1), col=1)
abline(h=sum(dr2), col=2)
# dlsigma
plot(ind[-1], dlsigma1)
points(ind[-1], dlsigma2, col=2, pch=2)
# dllambda
plot(ind[-1], dllambda1)
points(ind[-1], dllambda2, col=2, pch=2)
# dmu
plot(ind[-1], dmu1)
points(ind[-1], dmu2, col=2, pch=2)
# dlnu
plot(ind[-1], dlnu1)
points(ind[-1], dlnu2, col=2, pch=2)

# different but that is okay --- for the time being

# Optimisation
opt_trun <- nlminb(par, nll_mjd_trun, X=Xt, dt=dt, control=list(trace=1))
# gradient does not work...
opt_trun
cat(
    format(c(opt_trun$par[1],
             exp(opt_trun$par[2:3]),
             opt_trun$par[4],
             exp(opt_trun$par[5])),
           scientific = F, digits=3),
    sep = " & "
)
par_trun_orig <- c(opt_trun$par[1],
              exp(opt_trun$par[2:3]),
              opt_trun$par[4],
              exp(opt_trun$par[5]))

# SD
# build Hessian

dfun_trans <- function(par_orig, Xt, dt){
    # takes in actual parameters
    par_trans <- c(par_orig[1], log(par_orig[2]), log(par_orig[3]), par_orig[4], log(par_orig[5]))
    grad <- dnll_mjd_trun(par_trans, Xt, dt)
    dtrans <- c(1, 1/par_orig[2], 1/par_orig[3], 1, 1/par_orig[5])
    return(dtrans * grad)
}
dfun_trans(par_trun_orig, Xt, dt)

fun <- dfun_trans
x = par_trun_orig

#jacobian <- function(fun, x, dt=dt, Xt=Xt){
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
#    return(J)
#}
Hess <- J
#Hess <- jacobian(fun=dfun_trans, x=opt$par)
solve(Hess)
SD <- sqrt(diag(solve(Hess)))
SD

tab_vec <- c("")
for(i in 1:length(par)){
    tab_vec <- paste0(tab_vec,
                      cat(format(par_trun_orig[i], scientific = F, digits=3)),
                      cat(" (", format(SD[i], scientific = F, digits=3), ") & ", sep="")
                      , sep ="")
}
tab_vec
# concusion: truncated likelihood is not accurate enough in the case of derivative information.
# - thus meaningful SD cannot be obtained? ( But nlminb is able to do it!?)
# - grad also does not work during optimisation - only very accurate gradients are useful? ask Tore
# - explain method with AD grad and FD Hess - not accurate to get positive definite Hess

# SADDLEPOINT APPROXIMATION
nll_mjd_spa <- function(par, X, dt){
    nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 12, 1)$nll
    return(nll)
}
dnll_mjd_spa <- function(par, X, dt){
    res <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 12, 1)
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )
}
nll_mjd_spa(opt$par, Xt, dt)
opt$objective
dnll_mjd_spa(opt$par, Xt, dt)

# Optimisation
opt_spa <- nlminb(par, nll_mjd_spa, dnll_mjd_spa, X=Xt, dt=dt, control=list(trace=1))
opt_spa2 <- nlminb(par, nll_mjd_spa, X=Xt, dt=dt, control=list(trace=1))

par_spa_orig <- c(opt_spa$par[1],
                   exp(opt_spa$par[2:3]),
                   opt_spa$par[4],
                   exp(opt_spa$par[5]))


dfun_trans <- function(par_orig, Xt, dt){
    # takes in actual parameters
    par_trans <- c(par_orig[1], log(par_orig[2]), log(par_orig[3]), par_orig[4], log(par_orig[5]))
    grad <- dnll_mjd_spa(par_trans, Xt, dt)
    dtrans <- c(1, 1/par_orig[2], 1/par_orig[3], 1, 1/par_orig[5])
    return(dtrans * grad)
}
dfun_trans(par_spa_orig, Xt, dt)

fun <- dfun_trans
x = par_spa_orig

#jacobian <- function(fun, x, dt=dt, Xt=Xt){
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
#    return(J)
#}
Hess <- J
#Hess <- jacobian(fun=dfun_trans, x=opt$par)
solve(Hess)
SD <- sqrt(diag(solve(Hess)))
SD

tab_vec <- c("")
for(i in 1:length(par)){
    tab_vec <- paste0(tab_vec,
                      cat(format(par_spa_orig[i], scientific = F, digits=3)),
                      cat(" (", format(SD[i], scientific = F, digits=3), ") & ", sep="")
                      , sep ="")
}
tab_vec
