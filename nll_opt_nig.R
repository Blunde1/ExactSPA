# Simulation experiment for exact SPA methodology
# Distribution: Normal inverse Gaussian
# 01.04.2018

# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)
n = 1000

# Simulate datasets
set.seed(123)
seeds <- sample(10000,1000)
nig.list <- list()
n=200
for(i in 1:length(seeds)){
    nig.list[[i]]<-rNIG(n, c(chi, psi, mu, gamma), seed=seeds[i])
}
sum(which(sapply(1:length(seeds), function(i)any(is.na(nig.list[[i]])))))

# Estimate parameters
opt.list <- list() #par.est.list <- list()
load("optlist_nig_v1.RData")
for(i in 1:length(nig.list)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            opt.list[[i]] <- nlminb(par, nll_fun_nig, nll_grad_nig, X=nig.list[[i]], control=list(trace=1))
            save(opt.list, file="optlist_nig_v1.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

# Plot parameters
# Plot parameters ####
chi.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[1]))
psi.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[2]))
mu.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[3])
gamma.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[4])

par(mfrow=c(2,2))

plot(density(chi.ests),col="blue", main="r estimates")
hist(chi.ests, freq = FALSE, add=T)
abline(v=exp(par["lchi"]),col="red",lwd=3)
abline(v=meanexp((r.ests)),col="green",lwd=2)
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