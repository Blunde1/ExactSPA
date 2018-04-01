
nll_mjd(c(0,runif(10000,-0.05,0.05)), 1/250, .08, log(.1), log(100), -.001, log(.015), 12, 64)


par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015), x0=10) #jump intensity usually 100
time <- 4
N <- time*250
dt <- time/N


# Simulate mjd ####
#setwd("C:/Users/Berent/Projects/it-ift/implementation/simulation_experiment/mjd")
#source("../../simulation/Simulation_MJD_3.R")
library(ExactSPA)
X<-rMJD(N,time,x0=par["x0"],
                           r=par["r"], sigma=exp(par["lsigma"]), jump_intensity = exp(par["llambda"]),
                           mu=par["mu"], nu=exp(par["lnu"]),
                           seed=123)

plot(X)
nll_fun_mjd(par, log(X), dt)
nll_grad_mjd(par, log(X), dt)

opt <- nlminb(par[1:5], nll_fun_mjd, nll_grad_mjd, X=log(X), dt=time/N, control=list(trace=1))
opt$par
exp(opt$par)
opt

# Simulation experiment
set.seed(123)
seeds <- sample(10000,1000)
mjd.list <- list()
time <- 4
N <- time*250
par <- c(r=0.08, lsigma=log(0.1), llambda=log(100), mu=-0.001, lnu=log(0.015), x0=10) #jump intensity usually 100
for(i in 1:length(seeds)){
    mjd.list[[i]]<-rMJD(N,time,x0=par["x0"],
                               r=par["r"], sigma=exp(par["lsigma"]), jump_intensity = exp(par["llambda"]),
                               mu=par["mu"], nu=exp(par["lnu"]),
                               seed=seeds[i])
}
sum(which(sapply(1:length(seeds), function(i)any(is.na(mjd.list[[i]])))))

# # Problematic ones
# set.seed(4321)
# seeds.2 <- sample(10001:20000,100)
# mjd.list[[110]] <- mjd_process(N,time,x0=par["x0"],
#                                r=par["r"], sigma=exp(par["lsigma"]), jump_intensity = exp(par["llambda"]),
#                                mu=par["mu"], nu=exp(par["lnu"]),
#                                seed=seeds.2[1])


# Estimate parameters ####
opt.list <- list() #par.est.list <- list()
load("optimum_object_list_rcpp_4.RData")
for(i in 1:length(mjd.list)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            opt.list[[i]] <- nlminb(par[-c(6,7)], nll_fun_mjd, nll_grad_mjd, X=log(mjd.list[[i]]), dt=time/N, control=list(trace=1))
            save(opt.list, file="optimum_object_list_rcpp_4.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

# Plot parameters ####
r.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[1])
sigma.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[2]))
jump_intensity.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[3]))
mu.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[4])
nu.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[5]))

par(mfrow=c(3,2))

plot(density(r.ests),col="blue", main="r estimates")
hist(r.ests, freq = FALSE, add=T)
abline(v=par["r"],col="red",lwd=3)
abline(v=mean(r.ests),col="green",lwd=2)
legend("topright",legend = c("Histogram","Kernel", "True par", "Avg. par"),
       col=c("black","blue","red","green"),
       lwd=c(1,1,3,2))

plot(density(sigma.ests),col="blue",main="sigma estimates")
hist(sigma.ests, freq = FALSE,add=T)
abline(v=exp(par["lsigma"]),col="red",lwd=3)
abline(v=mean(sigma.ests), col="green", lwd=2)
legend("topright",legend = c("Histogram","Kernel", "True par", "Avg. par"),
       col=c("black","blue","red","green"),
       lwd=c(1,1,3,2))

plot(density(jump_intensity.ests),col="blue", main="lambda estimates")
hist(jump_intensity.ests, freq = FALSE,add=T)
abline(v=exp(par["llambda"]),col="red", lwd=3)
abline(v=mean(jump_intensity.ests), col="green", lwd=2)
legend("topright",legend = c("Histogram","Kernel", "True par", "Avg. par"),
       col=c("black","blue","red","green"),
       lwd=c(1,1,3,2))

plot(density(mu.ests),col="blue", main="mu estimates")
hist(mu.ests, freq = FALSE, add=T)
abline(v=par["mu"],col="red", lwd=3)
abline(v=mean(mu.ests), col="green", lwd=2)
legend("topright",legend = c("Histogram","Kernel", "True par", "Avg. par"),
       col=c("black","blue","red","green"),
       lwd=c(1,1,3,2))


plot(density(nu.ests),col="blue", main="nu estimates")
hist(nu.ests, freq = FALSE, add=T)
abline(v=exp(par["lnu"]),col="red", lwd=3)
abline(v=mean(nu.ests), col="green", lwd=2)
legend("topright",legend = c("Histogram","Kernel", "True par", "Avg. par"),
       col=c("black","blue","red","green"),
       lwd=c(1,1,3,2))

par(mfrow=c(1,1))

mean(kappa.ests)
sd(kappa.ests)
mean(alpha.ests)
sd(alpha.ests)
mean(exp(sigma.ests))

par(mfrow=c(1,1))
