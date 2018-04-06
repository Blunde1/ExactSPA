# Simulation experiment for exact SPA methodology
# Distribution: Normal inverse Gaussian
# 01.04.2018

library(ExactSPA)
setwd("C:/Users/Berent/Projects/it-ift/implementation v5")

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
save(nig.list, file="nig_data_sim.RData")

# Estimate parameters
opt.list <- list() #par.est.list <- list()
load("optlist_nig_v1.RData")
for(i in 10:length(nig.list)){
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

plot(density(chi.ests),col="blue", main="chi estimates")
hist(chi.ests, freq = FALSE, add=T)
abline(v=exp(par["lchi"]),col="red",lwd=3)
abline(v=mean(chi.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))

plot(density(psi.ests),col="blue", main="psi estimates")
hist(psi.ests, freq = FALSE, add=T)
abline(v=exp(par["lpsi"]),col="red",lwd=3)
abline(v=mean(psi.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))

plot(density(mu.ests),col="blue", main="mu estimates")
hist(mu.ests, freq = FALSE, add=T)
abline(v=par["mu"],col="red",lwd=3)
abline(v=mean(mu.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))

plot(density(gamma.ests),col="blue", main="gamma estimates")
hist(gamma.ests, freq = FALSE, add=T)
abline(v=exp(par["gamma"]),col="red",lwd=3)
abline(v=mean(gamma.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))


par(mfrow=c(1,1))

# Output format for latex table

cat(format(exp(par[1:2]), scientific=F, digits = 2),
    format(par[3:4], scientific=F, digits=2),
    sep= " & ")

cat(
    paste(
        as.numeric(format(mean(chi.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(chi.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(psi.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(psi.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(mu.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(mu.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(gamma.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(gamma.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    sep = " & "
)

as.numeric(format(mean(chi.ests), scientific = F, digits=2))
format(srep, scientific=F, digits=2)


# SPA ####
load("nig_data_sim.RData")
# Set parameters
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
par <- c(lchi=log(chi), lpsi=log(psi), mu=mu, gamma=gamma)

# Estimate parameters
opt.list <- list() #par.est.list <- list()
load("optlist_nig_spa_v1.RData")
for(i in 1:length(nig.list)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            opt.list[[i]] <- nlminb(par, nll_fun_nig, nll_grad_nig, type="SPA",
                                    X=nig.list[[i]], control=list(trace=1))
            save(opt.list, file="optlist_nig_spa_v1.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}
chi.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[1]))
psi.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[2]))
mu.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[3])
gamma.ests <- sapply(1:length(opt.list), function(i)opt.list[[i]]$par[4])


# Output format for latex table
cat(
    paste(
        as.numeric(format(mean(chi.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(chi.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(psi.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(psi.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(mu.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(mu.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(gamma.ests), scientific = F, digits=2)),
        paste("(",as.numeric(format(sd(gamma.ests), scientific = F, digits=2)),")", sep=""),
        sep = " "
    ),
    sep = " & "
)
