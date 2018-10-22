# Clocking SPA, ESPA512, IFT512
# Berent Lunde
# 30.08.2018

# Call library
library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")
ExactSPA::hello()

# Make sure that nll is standardized
nll_fun_nig <- function(par, X, type="ExactSPA"){
    lchi = par[1]
    lpsi = par[2]
    mu = par[3]
    gamma = par[4]
    if(type=="ExactSPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 2)$nll
    }else if(type=="SPA"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll
    }else if(type=="reSPA"){
        n <- 200
        x_sim <- sort(rNIG(n, c(exp(lchi), exp(lpsi), mu, gamma), seed = 1234))
        c <- sum(sapply(2:n, function(i) exp(-nll_nig(x_sim[i], lchi, lpsi, mu, gamma, 100, 512, 1)$nll)*(x_sim[i]-x_sim[i-1])))
        cat("value of renormalisation: ", c)
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 100, 512, 1)$nll + length(X)*log(c)
    }else if(type=="Simpson"){
        nll <- nll_nig(X, lchi, lpsi, mu, gamma, 1500, 512, 3)$nll
        
    }
    return(nll)
}

library(microbenchmark)
clock <- microbenchmark(dnorm(0), times = 10000)



# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

EX <- mu + sqrt(chi/psi)*gamma
VarX <- sqrt(chi/psi) + sqrt(chi/psi^3)*gamma^2

nll_fun_nig(par, EX, "ExactSPA")
nll_fun_nig(par, EX, "SPA")
nll_fun_nig(par, EX, "Simpson")

# Clock and convert from nanoseconds to microseconds
# At mean
espa1 <- microbenchmark(nll_fun_nig(par, EX, "ExactSPA"), times=1000, unit = "ns")$time / 1e9
spa1 <- microbenchmark(nll_fun_nig(par, EX, "SPA"), times=1000, unit = "ns")$time / 1e9
ift1 <- microbenchmark(nll_fun_nig(par, EX, "Simpson"), times=1000, unit = "ns")$time / 1e9
# 1 sd away from mean
espa2 <- microbenchmark(nll_fun_nig(par, EX + sqrt(VarX), "ExactSPA"), times=1000, unit = "ns")$time / 1e9
spa2 <- microbenchmark(nll_fun_nig(par, EX + sqrt(VarX), "SPA"), times=1000, unit = "ns")$time / 1e9
ift2 <- microbenchmark(nll_fun_nig(par, EX + sqrt(VarX), "Simpson"), times=1000, unit = "ns")$time / 1e9

# 2 sd away from mean
espa3 <- microbenchmark(nll_fun_nig(par, EX + 2*sqrt(VarX), "ExactSPA"), times=1000, unit = "ns")$time / 1e9
spa3 <- microbenchmark(nll_fun_nig(par, EX + 2*sqrt(VarX), "SPA"), times=1000, unit = "ns")$time / 1e9
ift3 <- microbenchmark(nll_fun_nig(par, EX + 2*sqrt(VarX), "Simpson"), times=1000, unit = "ns")$time / 1e9


# Plotting
espa <- data.frame(espa1 = espa1,
                   espa2 = espa2,
                   espa3 = espa3)
spa <- data.frame(spa1 = spa1,
                  spa2 = spa2,
                  spa3 = spa3)
ift <- data.frame(ift1 = ift1,
                  ift2 = ift2,
                  ift3 = ift3)
stargazer(espa, summary.stat = c("mean","sd"))


library(stargazer)
tf_stargazer <- function(cvec){
    library(readr)
    a <- cvec[12:14] # relevant information
    b <- strsplit(a, "&")
    c <- numeric(length(b))
    for(i in 1:length(b)){
        c[i] <- paste0( format(parse_number(b[[i]][2]), scientific = T, digits = 3), 
                        paste0("(", format(parse_number(b[[i]][3]), scientific = T, digits = 3), ")"), sep="")
    }
    return(c)
}
# format of table
clocking_tab <- data.frame(evaluation=c("$\\E[X]$", "$\\E[X] + \\sqrt{\\Var[X]}$", "$\\E[X] + 2\\sqrt{\\Var[X]}$"),
                           espa = tf_stargazer(stargazer(espa, summary.stat = c("mean","sd"))),
                           spa = tf_stargazer(stargazer(spa, summary.stat = c("mean","sd"))),
                           ift = tf_stargazer(stargazer(ift, summary.stat = c("mean","sd")))
                           )
for(i in 1:nrow(clocking_tab)){
    for(j in 1:ncol(clocking_tab)){
        cat(as.character(clocking_tab[i,j]))
        if(j < ncol(clocking_tab)){
            cat(" & ")   
        }else{
            cat("\\ \n")
        }
    }
}
