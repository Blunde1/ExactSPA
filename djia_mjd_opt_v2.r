# djia mjd opt 2
# new version that manually sets the gradient and likelihood in R
# 11.09.2018

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

pnll_fun_mjd <- function(par, X, map=NULL, dt){
    par <- c(par,map)
    par2 <- c(par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"])
    nll <- nll_mjd_wrap(par2, X, dt)
    return(nll)
}
pnll_fun_mjd_smpsn <- function(par, X, map=NULL, dt){
    par < c(par,map)
    nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 40, 252, 3)$nll
    return(nll)
}

pnll_fun_mjd_spa <- function(par, X, map=NULL, dt){
    par <- c(par,map)
    nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 64, 1)$nll
    return(nll)
}
pnll_fun_mjd_trun <- function(par, X, map=NULL, dt){
    par <- c(par,map)
    par2 <- c(par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"])
    nll <- nll_mjd_trun(par2, X, dt)
    return(nll)
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
        res <- nll_mjd(X[(i-1):i], dt, r, lsigma, llambda, mu, lnu, 12, 64, 2)
        grad <- grad +
            c(
                res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
            )
    }


    # return individual likelihood contributions
    return(grad)
}


# TEST THE FUNCTIONS
# Test
nll_fun_mjd(par, Xt, dt)
nll_mjd_wrap(par, Xt, dt)
pnll_fun_mjd(par, Xt, dt=dt)
nll_mjd_trun(par, Xt, dt)
pnll_fun_mjd_trun(par, Xt, dt=dt)
pnll_fun_mjd_smpsn(par, Xt, dt=dt)
pnll_fun_mjd_spa(par, Xt, dt=dt)
nll_grad_mjd(par, Xt, dt)
dnll_mjd_wrap(par, Xt, dt)

system.time(nll_fun_mjd(par, Xt, dt))
system.time(nll_mjd_wrap(par, Xt, dt))
system.time(nll_grad_mjd(par, Xt, dt))
system.time(dnll_mjd_wrap(par, Xt, dt))

# TEST OPTIMISATION
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

map <- c(opt$par["llambda"])
par2 <- par[!names(par)%in%names(map)]
opt2 <- nlminb(par2, pnll_fun_mjd, map=map, dt=dt, X=Xt, control=list(trace=1))
opt$par
opt2$par
map


# PROFILE LIKELIHOOD
# exact spa
llambda <- seq(2,8.5,length.out=30)
map <- c("llambda"=0)
par2 <- par[!names(par)%in%names(map)]
opt.list.profile <- list()
for(i in 1:length(llambda)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            map <- c("llambda"=llambda[i])
            opt.list.profile[[i]] <- nlminb(par2, pnll_fun_mjd, map=map,
                                            X=Xt, dt=dt, control=list(trace=1))
            save(opt.list.profile, file="optlist_mjd_profile_18092018.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}
# spa
opt.list.profile.spa <- list()
for(i in 2:length(llambda)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            map <- c("llambda"=llambda[i])
            opt.list.profile.spa[[i]] <- nlminb(par2, pnll_fun_mjd_spa, map=map,
                                                X=Xt, dt=dt, control=list(trace=1))
            save(opt.list.profile.spa, file="optlist_mjd_profile_spa20092018.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

# exact truncated
opt.list.profile.trun <- list()
#llambda <- seq(4.65, 4.91, length.out = 5)
for(i in 1:length(llambda)){
    tryCatch(
        {
            cat("iter: ", i, "\n")
            map <- c("llambda"=llambda[i])
            opt.list.profile.trun[[i]] <- nlminb(par2, pnll_fun_mjd_trun, map=map,
                                                 X=Xt, dt=dt, control = list(trace=1))
            save(opt.list.profile.trun, file="optlist_mjd_profile_trun17092018.RData")
            #par2 <- opt.list.profile.trun[[i]]$par
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}


# simpson
opt.list.profile.smpsn <- list()
for(i in 2:length(llambda)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            map <- c("llambda"=llambda[i])
            opt.list.profile.smpsn[[i]] <- nlminb(par2, pnll_fun_mjd_smpsn, map=map,
                                                X=Xt, dt=dt, control=list(trace=1))
            save(opt.list.profile.smpsn, file="optlist_mjd_profile_smpsn18102018.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}


# gauss
nll.gbm <- function(par, X, dt){
    mu <- par[1]; sigma <- par[2]
    nobs <- length(X)
    nll <- 0
    for(i in 2:nobs){
        nll <- nll - dnorm(X[i],
                           X[i-1]+(mu-0.5*sigma^2)*dt,
                           sqrt(dt)*sigma, log=TRUE)
    }
    return(nll)
}
par.gbm <- c(0.1,0.1)
opt.gbm <- nlminb(par.gbm, nll.gbm, dt=dt, X=Xt, control=list(trace=1))


setwd("C:/Users/Berent/Projects/it-ift/implementation v5")

load("optlist_mjd_profile_18092018.RData")
load("optlist_mjd_profile_spa20092018.RData")
load("optlist_mjd_profile_trun17092018.RData")
#load("optlist_mjd_profile_smpsn18102018.RData")
llambda <- seq(2,8.5,length.out=30)

nll.val <- sapply(1:length(opt.list.profile),function(i){opt.list.profile[[i]]$objective})
#nll.val.smpsn <- sapply(1:length(opt.list.profile.smpsn),function(i){opt.list.profile.smpsn[[i]]$objective})
nll.val.spa <- sapply(1:length(opt.list.profile.spa), function(i){opt.list.profile.spa[[i]]$objective})
nll.val.trun <- sapply(1:length(opt.list.profile.trun), function(i){opt.list.profile.trun[[i]]$objective})

plot(llambda, nll.val)
llambda2 <- llambda[1:40]
plot(llambda2,nll.val[1:40])


setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("mjd_pnll7.pdf", width=9, height=6)
plot(llambda[-c(1,2)], nll.val[-c(1,2)],
     ylim=c(min(c(nll.val,nll.val.spa, opt.gbm$objective)-1), max(c(nll.val,nll.val.spa, opt.gbm$objective)+1)),
     main="", xlab = expression(paste(log,(lambda))), ylab="Negative log-likelihood",
     type="l", lwd=4, lty=1)
lines(llambda[-c(1,2)], nll.val.trun[-c(1,2)], type="p", col=5, pch=3,lwd=2 )
lines(llambda[-c(1,2)], nll.val.spa[-c(1,2)], type="l", col=3, lwd=2, lty=1 )
lines(llambda[-c(1,2)], rep(opt.gbm$objective,length(llambda))[-c(1,2)], col=6, type="l", lwd=2, lty=6)
#lines(llambda[-c(1,2)], nll.val.smpsn[-c(1,2)], type="l", col=5, pch=5,lwd=3, lty=2 )

legend("bottomright", c("Saddlepoint adjusted IFT", "Truncated", "Saddlepoint approximation", "Gaussian GBM"),
       col=c(1,5,3,6), pch=c(NA,3,NA,NA), lwd=c(4,2,2,2), lty=c(1,NA,1,6))
dev.off()




