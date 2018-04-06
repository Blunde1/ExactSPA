library(ExactSPA)

df = 100
ncp = 40
n = 4000
x <- sort(rchisq(n, df, ncp))
hist(x, freq=F)
curve(dchisq(x, df, ncp), add=TRUE)

y <- rep(0, length(x))
for(i in 1:length(x)){
    y[i] <- exp(-nll_nchisq(as.vector(x[i]),log(df),log(ncp), 12.0, 64)$nll)
}
lines(x,y,col="blue")

-sum(dchisq(x, df, ncp, log=TRUE))
nll_nchisq(x, log(df), log(ncp), 20.0, 128)


par <- c(ldf=log(df), lncp = log(ncp))
opt <- nlminb(par, nll_fun_nchisq, nll_grad_nchisq, X=x, control=list(trace=1))
opt
opt$par
exp(opt$par)


# Simulate datasets
df = 100
ncp = 40
par <- c(ldf=log(df), lncp = log(ncp))

set.seed(123)
seeds <- sample(10000,1000)
nchisq.list <- list()
n=500
for(i in 1:length(seeds)){
    set.seed(seeds[i])
    nchisq.list[[i]]<-rchisq(n, df, ncp)
}
sum(which(sapply(1:length(seeds), function(i)any(is.na(nchisq.list[[i]])))))
save(nchisq.list, file="nchisq_data_sim.RData")
load("nchisq_data_sim.RData")

# Estimate parameters
opt.list <- list() #par.est.list <- list()
load("optlist_nchisq_v1.RData")
for(i in 1:length(nchisq.list)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            opt.list[[i]] <- nlminb(par, nll_fun_nchisq, nll_grad_nchisq, X=nchisq.list[[i]], control=list(trace=1))
            save(opt.list, file="optlist_nchisq_v1.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

# Plot parameters ####
df.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[1]))
ncp.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[2]))

par(mfrow=c(1,2))

plot(density(df.ests),col="blue", main="df estimates")
hist(df.ests, freq = FALSE, add=T)
abline(v=exp(par["ldf"]),col="red",lwd=3)
abline(v=mean(df.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))

plot(density(ncp.ests),col="blue", main="ncp estimates")
hist(ncp.ests, freq = FALSE, add=T)
abline(v=exp(par["lncp"]),col="red",lwd=3)
abline(v=mean(ncp.ests),col="green",lwd=2)
legend("topright",legend = c("Kernel","Histogram", "True val", "Avg. est."),
       col=c("blue", "black","red","green"),
       lwd=c(1,1,3,2))

par(mfrow=c(1,1))

# Output format for latex table

cat(format(exp(par[1:2]), scientific=F, digits = 2),
        sep= " & ")

cat(
    paste(
        as.numeric(format(mean(df.ests), scientific = F, digits=4)),
        paste("(",as.numeric(format(sd(df.ests), scientific = F, digits=4)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(ncp.ests), scientific = F, digits=4)),
        paste("(",as.numeric(format(sd(ncp.ests), scientific = F, digits=4)),")", sep=""),
        sep = " "
    ),
    sep = " & "
)

# SPA ####
load("nchisq_data_sim.RData")
df = 100
ncp = 40

# Estimate parameters
opt.list <- list() #par.est.list <- list()
load("optlist_nchisq_spa_v1.RData")
for(i in 1:length(nchisq.list)){
    tryCatch(
        {
            cat("iter:",i,"\n")
            opt.list[[i]] <- nlminb(par, nll_spa_fun_nchisq, nll_spa_grad_nchisq, type="SPA", X=nchisq.list[[i]], control=list(trace=1))
            save(opt.list, file="optlist_nchisq_spa_v1.RData")
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
        }
    )
}

# Output format for latex table
df.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[1]))
ncp.ests <- sapply(1:length(opt.list), function(i)exp(opt.list[[i]]$par[2]))

cat(format(exp(par[1:2]), scientific=F, digits = 2),
    sep= " & ")

cat(
    paste(
        as.numeric(format(mean(df.ests), scientific = F, digits=4)),
        paste("(",as.numeric(format(sd(df.ests), scientific = F, digits=4)),")", sep=""),
        sep = " "
    ),
    paste(
        as.numeric(format(mean(ncp.ests), scientific = F, digits=4)),
        paste("(",as.numeric(format(sd(ncp.ests), scientific = F, digits=4)),")", sep=""),
        sep = " "
    ),
    sep = " & "
)
