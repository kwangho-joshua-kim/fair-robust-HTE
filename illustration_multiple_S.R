library(MASS)
library(latex2exp)
library(matrixcalc)
library(mbend)
library(Matrix)
library(Rmosek)
library(ggplot2)
library(ranger)
library(quadprog)
expit <- function(x){ exp(x)/(1+exp(x))}
source("utils.R")

####################################################################################
# Toy dataset for multiple sensitive features (S1, S2)
####################################################################################
set.seed(2023)
n = 1000
S1 <- rbinom(n = n, size = 1, prob = 0.5)
S2 <- rbinom(n = n, size = 1, prob = 0.5)
X11 <- rnorm(sum(S1), 1, 1/1.5)
X10 <- rnorm(n-sum(S1), -1, 1/1.5)
X21 <- rnorm(sum(S2), 1, 1/1.5)
X20 <- rnorm(n-sum(S2), -1, 1/1.5)
X3 <- rnorm(n, 0, 1)
X <- matrix(NA, nrow=n, ncol=3)
X[S1==1,1] <- X11
X[S1==0,1] <- X10
X[S2==1,2] <- X21
X[S2==0,2] <- X20
X[,3] <- X3
pa <- expit(X3)
A <- rbinom(n = n, size = 1, prob = pa)
e <- rnorm(n,0,0.1)
Y1 = X[,1]^3/4 + X[,2]/2 + X[,3] + e
Y0 = X[,3] + e
Y = A*Y1 + (1-A)*Y0
dat <- as.data.frame(cbind(X,S1,S2,A,Y))
colnames(dat) <- c("X1", "X2", "X3", "S1", "S2", "A", "Y")
tau <- X[,1]^3/4 + X[,2]/2


####################################################################################
# some plots
####################################################################################
# density plot: CATE vs sensitive features
pal1 <- c("#CC6677", "#88CCEE")
pal2 <- c("#DDCC77", "#44AA99")
d1.s0.inf <- density(tau[S1==0])
d1.s1.inf <- density(tau[S1==1])
d2.s0.inf <- density(tau[S2==0])
d2.s1.inf <- density(tau[S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y,d2.s0.inf$y,d2.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y,d2.s0.inf$y,d2.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x,d2.s0.inf$x,d2.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x,d2.s0.inf$x,d2.s1.inf$x))
den <- density(tau)
par(mfrow=c(2,1))
# CATE density and S1 
plot(den, xlab = TeX(r'(${\tau}$)'), ylab="", main=TeX(r'($\tau$ vs $S_1$)'), cex.main=2, cex.lab=2)
tau1 <- tau[S1==1]
tau0 <- tau[S1==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal1[1],0.5), pch=19)
points(tau0, rep(-0.006,length(tau0)), col=alpha(pal1[2],0.5), pch=19)
legend("topright", legend=c(TeX(r'($S_1=0$)'), TeX(r'($S_1=1$)')),
       col=pal1, cex=2, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal1[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal1[1], xlim=c(x.min,x.max))
# S2
plot(den, xlab = TeX(r'(${\tau}$)'), ylab="", main=TeX(r'($\tau$ vs $S_2$)'), cex.main=2, cex.lab=2)
tau1 <- tau[S2==1]
tau0 <- tau[S2==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(-0.006,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topright", legend=c(TeX(r'($S_2=0$)'), TeX(r'($S_2=1$)')),
       col=pal2, cex=2, pch=19)
abline(v=weighted.mean(d2.s0.inf$x, d2.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d2.s1.inf$x, d2.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
par(mfrow=c(1,1))

# bar plot: OTR vs sensitive features
par(mfrow=c(2,1))
par(mar = c(3.5, 4.5, 3, 7))
# S1
ts1 = sum(S1==1 & tau > 0)
ts0 = sum(S1==0 & tau > 0)
uts1 = sum(S1==1 & tau < 0)
uts0 = sum(S1==0 & tau < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal1 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\tau>0$)'),TeX(r'($\tau<0$)')), main = TeX(r'(OTR w.r.t.$S_1$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S_1=1$)'), TeX(r'($S_1=0$)')),
        col=pal1, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)
# S2
ts1 = sum(S2==1 & tau > 0)
ts0 = sum(S2==0 & tau > 0)
uts1 = sum(S2==1 & tau < 0)
uts0 = sum(S2==0 & tau < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal2 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\tau>0$)'),TeX(r'($\tau<0$)')), main = TeX(r'(OTR w.r.t.$S_2$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S_2=1$)'), TeX(r'($S_2=0$)')),
        col=pal2, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)
par(mfrow=c(1,1))


####################################################################################
# fair and robust CATE estimation
####################################################################################
nm.s.features <- c("S1", "S2")
fr.cate <- fr_cate(dat, nm.s.features,
                   nm.l.factors = c(NA,NA),
                   solver="mosek",
                   nuisance.est='rf',
                   delta=0.005, sqr=TRUE, interactions=TRUE)
b.mat <- fr.cate$b.mat
beta.hat <- fr.cate$beta.hat
tau.hat <- b.mat %*% beta.hat
# RMSE
sqrt(mean((tau-tau.hat)^2))


####################################################################################
# diagnostics 
####################################################################################
tau <- tau.hat
# density plot: CATE vs sensitive features
pal1 <- c("#CC6677", "#88CCEE")
pal2 <- c("#DDCC77", "#44AA99")
d1.s0.inf <- density(tau[S1==0])
d1.s1.inf <- density(tau[S1==1])
d2.s0.inf <- density(tau[S2==0])
d2.s1.inf <- density(tau[S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y,d2.s0.inf$y,d2.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y,d2.s0.inf$y,d2.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x,d2.s0.inf$x,d2.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x,d2.s0.inf$x,d2.s1.inf$x))
den <- density(tau)
par(mfrow=c(2,1))
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\tau}}$)'), ylab="", main=TeX(r'($\hat{\tau}$ vs $S_1$)'), cex.main=2, cex.lab=2)
tau1 <- tau[S1==1]
tau0 <- tau[S1==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal1[1],0.5), pch=19)
points(tau0, rep(-0.006,length(tau0)), col=alpha(pal1[2],0.5), pch=19)
legend("topright", legend=c(TeX(r'($S_1=0$)'), TeX(r'($S_1=1$)')),
       col=pal1, cex=2, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal1[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal1[1], xlim=c(x.min,x.max))
# S2
plot(den, xlab = TeX(r'(${\hat{\tau}}$)'), ylab="", main=TeX(r'($\hat{\tau}$ vs $S_2$)'), cex.main=2, cex.lab=2)
tau1 <- tau[S2==1]
tau0 <- tau[S2==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(-0.006,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topright", legend=c(TeX(r'($S_2=0$)'), TeX(r'($S_2=1$)')),
       col=pal2, cex=2, pch=19)
abline(v=weighted.mean(d2.s0.inf$x, d2.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d2.s1.inf$x, d2.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
par(mfrow=c(1,1))

# bar plot: OTR vs sensitive features
par(mfrow=c(2,1))
par(mar = c(3.5, 4.5, 3, 7))
# S1
ts1 = sum(S1==1 & tau > 0)
ts0 = sum(S1==0 & tau > 0)
uts1 = sum(S1==1 & tau < 0)
uts0 = sum(S1==0 & tau < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal1 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{\tau}>0$)'),TeX(r'($\hat{\tau}<0$)')), main = TeX(r'(OTR w.r.t.$S_1$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S_1=1$)'), TeX(r'($S_1=0$)')),
        col=pal1, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)
# S2
ts1 = sum(S2==1 & tau > 0)
ts0 = sum(S2==0 & tau > 0)
uts1 = sum(S2==1 & tau < 0)
uts0 = sum(S2==0 & tau < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal2 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{\tau}>0$)'),TeX(r'($\hat{\tau}<0$)')), main = TeX(r'(OTR w.r.t.$S_2$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S_2=1$)'), TeX(r'($S_2=0$)')),
        col=pal2, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)
par(mfrow=c(1,1))