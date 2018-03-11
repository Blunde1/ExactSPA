library(ExactSPA)
library(Quandl)

# Load data
start_date <- "2003-01-01"; end_training <- "2005-01-01"; 
DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
#DJIA <- Quandl("BCB/UDJIAD1")
DJIA <- DJIA[rev(rownames(DJIA)),]
plot(DJIA,type="l")
Xt=log(DJIA$Value)

# Start params
par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015))
dt <- 1/252

# Test
nll_fun_mjd(par, Xt, dt)
nll_grad_mjd(par, Xt, dt)

# Estimate parameters
opt <- nlminb(par, nll_fun_mjd, nll_grad_mjd, X=Xt, dt=dt, control=list(trace=1))
opt
opt$par
exp(opt$par)
