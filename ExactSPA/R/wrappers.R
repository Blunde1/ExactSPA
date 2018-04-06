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

nll_fun_mjd <- function(par, X, dt, type="ExactSPA"){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    if(type=="ExactSPA"){
        nll <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64, 1)$nll
    }
    return(nll)
}

nll_grad_mjd <- function(par, X, dt, type="ExactSPA"){
    r = par[1]
    lsigma = par[2]
    llambda = par[3]
    mu = par[4]
    lnu = par[5]
    if(type=="ExactSPA"){
        res <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64, 2)
    }else if(type=="SPA"){
        res <- nll_mjd(X, dt, r, lsigma, llambda, mu, lnu, 12, 64, 1)
    }
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )
}

nll_fun_nchisq <- function(par, X, type="ExactSPA"){
    ldf = par[1]
    lncp = par[2]
    if(type=="ExactSPA"){
        nll <- nll_nchisq(X, ldf, lncp, 12, 64, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nchisq(X, ldf, lncp, 12, 64, 1)$nll
    }
    return(nll)
}

nll_grad_nchisq <- function(par, X, type="ExactSPA"){
    ldf = par[1]
    lncp = par[2]
    if(type=="ExactSPA"){
        res <- nll_nchisq(X, ldf, lncp, 12, 64, 2)
    }else if(type=="SPA"){
        res <- nll_nchisq(X, ldf, lncp, 12, 64, 1)
    }
    c(
        res$ldf_grad, res$lncp_grad
    )
}


nll_fun_nig <- function(par, X, type="ExactSPA"){
    lchi = par[1]
    lpsi = par[2]
    mu = par[3]
    gamma = par[4]
    if(type=="ExactSPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll
    }
    return(nll)

}

nll_grad_nig <- function(par, X, type="ExactSPA"){
    lchi = par[1]
    lpsi = par[2]
    mu = par[3]
    gamma = par[4]
    if(type=="ExactSPA"){
        res <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 2)
    }else if(type=="SPA"){
        res <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)
    }
    c(
        res$lchi_grad, res$lpsi_grad, res$mu_grad, res$gamma_grad
    )
}

