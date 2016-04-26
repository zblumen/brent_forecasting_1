rm(list=ls());
#install.packages("rstan")
#install.packages("rstudioapi")
#install.packages("gtools")

library(rstan);
library(rstudioapi);
library(gtools);

setwd("D:/oil_pred_3/Brent_Futures_Monthly_6");
source("adapt_filter.R");
samp.info=list(sims=NA,burn.in=NA,chains=NA,cores=NA,adapt_delta=NA);
samp.info$sims=1e4*6; #number of samples for HMC
samp.info$burn.in=2e3; 
samp.info$chains=1; #number of Markov chains, NOTE: using more than 1 chain will likely cause label switching
samp.info$cores=1; #number of cores (actually logical processors?) to utiize in parellel
samp.info$adapt_delta=.8
horizons=5:6;
G=(samp.info$sims-samp.info$burn.in)*samp.info$chains;
n.floor=G/10;
model.names=c("normal_ls_mh.stan","normal_ls_ar_mh.stan",
              'mixture_ls_mh.stan',"mixture_ls_ar_mh.stan",
              "markov_mixture_ls_new_mh.stan",
              "markov_mixture_ls_new_ar_mh.stan")
N=6
h=4

stan.mod = stan_model(model.names[N]);

  ###sample from WTI PRIOR (run STAN model)
  d=read.csv(paste("h_",h,".csv",sep=''),as.is=TRUE);
  TT = sum(d$brent.ind);
  M=TT-h+1; #number of 
  T = nrow(d)-TT;
  fe=log(d$exp.spot)-log(d$fut);
  fit=sampling(stan.mod,data=list(H=h,T=T,TT=TT,y=fe),
               iter=samp.info$sims,warmup=samp.info$burn.in,
               chains = samp.info$chains,cores=samp.info$cores, 
               control=list(adapt_delta = samp.info$adapt_delta));
  ext=extract(fit);
  ext$out_prob=exp(ext$out_log_prob);
  ##run adaptive filter (as a function)
  s1=adapt.filter(ext$out_prob,n.floor)
  #collect predictive likelihood, PIT and Samples
  p=numeric(M);
  u=numeric(M);
  Samps=matrix(NA,G,M);
  
  for(m in 1:M){
    p[m]=crossprod(s1$w[,m],ext$pred_like[s1$ind[,m],m])
    u[m]=crossprod(s1$w[,m],ext$unif_like[s1$ind[,m],m])
    Samps[,m]=log(d$fut[T+m+h-1])+
      sample(ext$y_sims[s1$ind[,m],m],G,replace=TRUE,prob=s1$w[,m])
  }
pp=colMeans(ext$pred_like);
sum(log(pp))
sum(log(p))
cc=log(c(p,pp))
plot(log(pp),type='l',ylim=c(min(cc),max(cc)))
points(log(p),type='l',col=2,lty=2)

cc=c(cumsum(log(p)),cumsum(log(pp)))
plot(cumsum(log(pp)),type='l',ylim=c(min(cc),max(cc)))
points(cumsum(log(p)),type='l',col=2,lty=2)
  ##go forward and compare
  T2=T+170;
  TT2=TT-170;
  
  M2=TT2-h+1;

  fit2=sampling(stan.mod,data=list(H=h,T=T2,TT=TT2,y=fe),
               iter=samp.info$sims,warmup=samp.info$burn.in,
               chains = samp.info$chains,cores=samp.info$cores, 
               control=list(adapt_delta = samp.info$adapt_delta));
  ext2=extract(fit2,c('out_log_prob','pred_like','unif_like'));
  ext2$out_prob=exp(ext2$out_log_prob);
  ##run adaptive filter (as a function)
  s2=adapt.filter(ext2$out_prob,n.floor)
  #collect predictive likelihood, PIT and Samples
  p2=numeric(M2);
  u2=numeric(M2)

  for(m in 1:M2){
    p2[m]=crossprod(s2$w[,m],ext2$pred_like[s2$ind[,m],m])
    u2[m]=crossprod(s2$w[,m],ext2$unif_like[s2$ind[,m],m])
  }
