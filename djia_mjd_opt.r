library(ExactSPA)
library(Quandl)

# Load data
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")
start_date <- "2003-01-01"; end_training <- "2005-01-01";
#DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
#DJIA <- DJIA[rev(rownames(DJIA)),]
#plot(DJIA,type="l")
#save(DJIA, file="djia_01012003_01012005.RData")
load("djia_01012003_01012005.RData")


# new data 10.09.2018
start_date <- "2000-01-01"; end_training <- "2008-01-01";
DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
DJIA <- DJIA[rev(rownames(DJIA)),]
plot(DJIA,type="l")

Xt=log(DJIA$Value)

# Start params
par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015))
dt <- 1/252


# Profile likelihood - MJD
pnll_fun_mjd <- function(par, X, map=NULL, dt, type="ExactSPA"){
    par <- c(par,map)
    if(type=="ExactSPA"){
        nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 64, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 64, 1)$nll
    }else if(type=="Simpson"){
        nll <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 700, 64, 3)$nll
    }
    return(nll)
}

# Using the gradient does not work well with optimisation
pnll_grad_mjd <- function(par, X, map=NULL, dt, type="ExactSPA"){
    par <- c(par,map)
    map_indices <- which(par %in% map)
    if(type=="ExactSPA"){
        res <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 64, 2)
    }else if(type=="SPA"){
        res- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 12, 64, 1)
    }else if(type=="Simpson"){
        res <- nll_mjd(X, dt, par["r"], par["lsigma"], par["llambda"], par["mu"], par["lnu"], 700, 64, 3)
    }
    c(
        res$r_grad, res$sigma_grad, res$lambda_grad, res$mu_grad, res$nu_grad
    )[-c(map_indices)]
}



# Optimisation over all parameters - different methods
# Test
nll_fun_mjd(par, Xt, dt)
nll_grad_mjd(par, Xt, dt)

# Estimate parameters
opt <- nlminb(par, nll_fun_mjd, nll_grad_mjd, X=Xt, dt=dt, control=list(trace=1))
opt
opt$par
exp(opt$par)

# map for lambda - matches with optimisation over all param
map <- c(opt$par["llambda"])
par2 <- par[!names(par)%in%names(map)]
opt2 <- nlminb(par2, pnll_fun_mjd, map=map, dt=dt, X=Xt, control=list(trace=1))
opt$par
opt2$par
map

llambda <- seq(2,8.5,length.out=50)
map <- c("llambda"=0)
par2 <- par[!names(par)%in%names(map)]
opt.list.profile <- list()
for(i in 1:length(llambda)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            map <- c("llambda"=llambda[i])
            opt.list.profile[[i]] <- nlminb(par2, pnll_fun_mjd, map=map, type="ExactSPA",
                                    X=Xt, dt=dt, control=list(trace=1))
            save(opt.list.profile, file="optlist_mjd_profile.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

nll.val <- sapply(1:length(opt.list.profile),function(i){opt.list.profile[[i]]$objective})
plot(llambda, nll.val)
llambda2 <- llambda[1:40]
plot(llambda2,nll.val[1:40])

# Notices that we previously has found a local opt, not global
par.opt <- opt.list.profile[[which.min(nll.val)]]$par
par.opt <- c(par.opt[1:2],llambda[which.min(nll.val)], par.opt[3:4])
opt2 <- nlminb(par.opt, nll_fun_mjd, nll_grad_mjd, X=Xt, dt=dt, control=list(trace=1))
opt2
opt2$par
exp(opt2$par)

# Finds spa values for llambda2
opt.list.profile.spa <- list()
for(i in 1:length(llambda2)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            map <- c("llambda"=llambda2[i])
            opt.list.profile.spa[[i]] <- nlminb(par2, pnll_fun_mjd, map=map, type="SPA",
                                            X=Xt, dt=dt, control=list(trace=1))
            save(opt.list.profile.spa, file="optlist_mjd_profile_spa.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}
nll.val.spa <- sapply(1:length(opt.list.profile.spa), function(i){opt.list.profile.spa[[i]]$objective})
load("optlist_mjd_profile_spa.RData")
plot(llambda2, -nll.val[1:40], type="b")
lines(llambda2, -nll.val.spa, col="red")

# Add a normal distribution
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

load("optlist_mjd_profile.RData")
load("optlist_mjd_profile_spa.RData")

setwd("C:/Users/Berent/Projects/it-ift/implementation/plotting/test_plots")
pdf("mjd_pnll.pdf", width=7, height=4+1/3)
plot(llambda2, nll.val[1:40], type="l",ylim=c(min(c(nll.val[1:40],nll.val.spa,opt.gbm$objective)-1),
                                               max(c(nll.val[1:40],nll.val.spa,opt.gbm$objective)+1)),
     main="", xlab = expression(paste(log,(lambda))), ylab="Negative log-likelihood",
     lwd=3, lty=2)
lines(llambda2, nll.val.spa, type="l", col="red", pch=4,lwd=3, lty=1 )
lines(llambda2, rep(opt.gbm$objective,length(llambda2)), col="blue", type="b",pch=2, lwd=2)
legend("bottomleft", c("Exact","SPA","GBM"), lty=c(2,1,NA), col=c("black","red","blue"),
       lwd=c(3,3,2), pch=c(NA,NA,2))
dev.off()

# Formatted estimated parameters (for latex)
# Global
cat(
    format(c(opt2$par[1],
             exp(opt2$par[2:3]),
             opt2$par[4],
             exp(opt2$par[5])),
           scientific = F, digits=4),
    sep = " & "
)

# Local
cat(
    format(c(opt$par[1],
             exp(opt$par[2:3]),
             opt$par[4],
             exp(opt$par[5])),
           scientific = F, digits=4),
    sep = " & "
)
