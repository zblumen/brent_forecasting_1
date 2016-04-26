rm(list=ls());
#install.packages("rstan")
#install.packages("rstudioapi")
#install.packages("gtools")

library(rstan);
library(rstudioapi);
library(gtools);

setwd("D:/oil_pred_3/Brent_Futures_Monthly_Git");
source("adapt_filter.R");
samp.info=list(sims=NA,burn.in=NA,chains=NA,cores=NA,adapt_delta=NA);
samp.info$sims=1e4; #number of samples for HMC
samp.info$burn.in=2e3; 
samp.info$chains=1; #number of Markov chains, NOTE: using more than 1 chain will likely cause label switching
samp.info$cores=1; #number of cores (actually logical processors?) to utiize in parellel
samp.info$adapt_delta=.9
G=(samp.info$sims-samp.info$burn.in)*samp.info$chains;
n.floor=G/5;
model.names=c("normal_ls_mh.stan","normal_ls_ar_mh.stan",
              'mixture_ls_mh.stan',"mixture_ls_ar_mh.stan",
              "markov_mixture_ls_new_mh.stan",
              "markov_mixture_ls_new_ar_mh.stan")
N=1
stan.mod = stan_model(model.names[N]);
MN=sub("\\.","_",model.names[N]);
for(h in 1:12){
  message("==============================================")
  message(paste("Forecast Horizon",h,"of",12))
  message("==============================================")
  samp.info$h=h;
  ########################################################
  d=read.csv(paste("h_",samp.info$h,".csv",sep=''),as.is=TRUE);
  samp.info$TT = sum(d$brent.ind);
  samp.info$T = nrow(d)-samp.info$TT;
  samp.info$y=log(d$exp.spot)-log(d$fut);
  dates=tail(d$exp.date, samp.info$TT-h+1)
  ########################################################
  s1=adapt.filter2(stan.mod,samp.info,n.floor);
  ########################################################
  message("writting w.out")
  w.out=as.data.frame(s1$w);  names(w.out)=dates;
  write.csv(w.out,paste("nw_",MN,"_",h,".csv",sep=''),row.names=FALSE);

  message("writting y.out")  
  y.out=as.data.frame(s1$y);  names(y.out)=dates;
  write.csv(y.out,paste("y_",MN,"_",h,".csv",sep=''),row.names=FALSE);
  
  message("writting likleihoods and PIT values")  
  exp.out=data.frame(p=s1$p,u=s1$u,po=s1$non_updated_p,
                     uo=s1$non_updated_u,n.eff=s1$n.eff)
  row.names(exp.out)=dates;
  write.csv(exp.out,paste("pu_",MN,"_",h,".csv",sep=''));
  ########################################################
}
