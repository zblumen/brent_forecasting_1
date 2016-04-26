rm(list=ls());
setwd("D:/oil_pred_3/Brent_Futures_Monthly_6");

model.names=c("normal_ls_mh_stan","normal_ls_ar_mh_stan",
  "mixture_ls_mh_stan","mixture_ls_ar_mh_stan",
  "markov_mixture_ls_new_mh_stan",
  "markov_mixture_ls_new_ar_mh_stan")
horizons = 1:12

for(h in horizons){
  p.mat=NULL;
  u.mat=NULL;
  windows()
  par(mfrow=c(1,1))
  for(i in model.names){
    in.file = paste("pu_",i,"_",h,".csv",sep='')  
    temp.frame=read.csv(in.file,header=TRUE,row.names=1)
    p.mat=cbind(p.mat,temp.frame$p,temp.frame$po);
    u.mat=cbind(u.mat,temp.frame$u,temp.frame$uo);

    #plot(density(qnorm(temp.frame$u)),main=i)
    #plot(density(qnorm(temp.frame$uo)),main=paste(i,"Not_up"))
  }
  uu=rowMeans(u.mat)
 # plot(density(qnorm(uu)),main=h)
  #points(sort(qnorm(uu)),dnorm(sort(qnorm(uu))),
     #    type='l',col=2,lty=2)
  plot(sort(qnorm(uu)),(1:length(uu))/length(uu),type='l',main=h)
  points(sort(qnorm(uu)),pnorm(sort(qnorm(uu))),type='l',col=2)
  print(shapiro.test(qnorm(uu)))
  print(mean(qnorm(uu)))
  p.frame = data.frame(p.mat);
 # names(p.frame)=model.names;
  row.names(p.frame)=row.names(temp.frame);
  write.csv(p.frame,paste("p_",h,".csv",sep=''));
  
  u.frame = data.frame(u.mat);
 # names(u.frame)=model.names;
  row.names(u.frame)=row.names(temp.frame);
  write.csv(u.frame,paste("u_",h,".csv",sep=''));
}