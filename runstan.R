library(rstan)
rstan_options(auto_write = TRUE)
source('~/Dropbox/monkeybris/rscript/setup.R')
fv <- cullSNP(fv$X,fv$SNP,fv$MON,mthresh=0.25)
X <- fv$X[,-SNPdelink(fv$SNP)]
inds <- which(!is.na(X),arr.ind = T)

#flat
standat <- list(N=nrow(inds),n=nrow(X),m=ncol(X),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=X[inds])

fit <- stan(file="flipperbaby.stan",data=standat,iter=200,chains = 3,cores = 3)



#beta-heirarchical
standat <- list(N=nrow(inds),n=nrow(X),m=ncol(X),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=X[inds],v0=0.075)

fit3 <- stan(file="flipperking.stan",data=standat,iter=30,chains = 1,cores = 4)


#logit-link heirarchical
standat <- list(N=nrow(inds),n=nrow(X),m=ncol(X),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=X[inds],
                mu0=10,s0=10)
getpars <- c("maf","f_mu","f_s","f")
fit5 <- stan(file="flipperking_logit.stan",data=standat,pars = getpars,
             iter=200,chains = 3,cores = 3)

plt <- stan_plot(fit5,pars="f",point_est="mean",ci_level=0.5)
plt$data$y[order(plt$data$mean)] <- 1:387
plt <- plt + theme_minimal() + scale_y_discrete(labels="",breaks=c(),name="monkey") +
  scale_x_continuous(limits=c(0,0.25),name="F (SNP)")


#truncated normal heirarchical
standat <- list(N=nrow(inds),n=nrow(X),m=ncol(X),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=X[inds],s0=0.25)
fit6 <- stan(file="flipperking_normtrunc.stan",data=standat,
             iter=500,chains = 3,cores = 3)
