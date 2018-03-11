# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Welcome to ExactSPA functionality!")
}

nll_fun_mjd <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64)$nll
}

nll_grad_mjd <- function(par, X, dt){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    res <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64)
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )
}

nll_fun_nchisq <- function(par, X){
    ldf = par[1]
    lncp = par[2]
    nll_nchisq(X, ldf, lncp, 12, 64)$nll
}

nll_grad_nchisq <- function(par, X){
    ldf = par[1]
    lncp = par[2]
    res <- nll_nchisq(X, ldf, lncp, 12, 64)
    c(
        res$ldf_grad, res$lncp_grad
    )
}
