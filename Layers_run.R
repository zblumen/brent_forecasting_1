rm(list=ls())
setwd('D:/oil_pred_3/Brent_Futures_Monthly_6');
source("Layers_functions.R")

th=10;
tnl=3;
d=read.csv(paste('p_',th,'.csv',sep=''),header=TRUE,as.is=TRUE)
tp=as.matrix(d[,2:13])
talpha_def=list();
for(i in 1:tnl){
  talpha_def[[i]]=list();
  talpha_def[[i]]$alpha=2;
  talpha_def[[i]]$def.w=1;  
  talpha_def[[i]]$control=list(reltol=1e-32);
}


mm=Layered_scores(tp,tnl,th,talpha_def)

colSums(tail(mm$lp[-c(1:(th*tnl)),],10))

aa=colSums(tail(mm$lp[-c(1:(th*tnl)),],10))
rank(-aa)
