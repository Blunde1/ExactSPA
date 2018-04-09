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

llambda <- seq(2,10,length.out=50)
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
