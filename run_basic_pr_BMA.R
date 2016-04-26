rm(list=ls());
z.ddir.unorm=function(x,a){ #unormalized log dir density
  sum((a-1)*log(x))
}
m.dir=function(a){ #dericlet mode
  (a-1)/(sum(a)-length(a))
}
setwd('D:/oil_pred_3/Brent_Futures_Monthly_6');
w.fun=function(r){ #invertabel function from R^n to unit (n+1) weight space
  b=exp(c(0,r));
  b/sum(b)
}
winv.fun=function(w){
  #w[1]=exp(0)/sum(b)=>
  #sum.b=1/w[1];
  #w=exp(r)/sum.b=> r=
  #log(w[2:length(w)]*sum.b)=>
  #log(w[2:length(w)])+log(sum.b)=>
  #log(w[2:length(w)])+log(1)-log(w[1])+>
  log(w[2:length(w)])-log(w[1])
}
BMA_wheights=function(p,prior_prob){
  if(is.matrix(p)){
    like=apply(p,2,prod);    
  }else{
    like=p;
  }
  target=like*prior_prob;
  total_like=sum(target)
  target/total_like;
}
n=12;
alpha=rep(2,n)
SHRINK=1; ##1=> Shrinkage is used, 0 => no Shrinkage 
p.index=1:n
run_code='';

t_log_score=matrix(NA,nrow=n+3,ncol=12); #total log scores for ind models + equal weight + optimal pool

reltol=1e-32 ;
equal_weight=rep(1/n,n);
for(h in 1:12){ #h is the forecast horizon
  message("==================================")
  message(paste("horizon",h));
  message("----------------------------------")
  infile=paste('p_',h,'.csv',sep=''); #file containing pred.likes
  p_data=read.csv(infile,row.names=1,header=TRUE,as.is=TRUE)
  p=as.matrix(p_data[,p.index])
  T=nrow(p);
  
  log_score=matrix(NA,T,n+3); #log score for ind models + equal weight+ alt. weights + optimal pool
  weights_pool=matrix(NA,T,n);
  weights_BMA=matrix(NA,T,n)
  par.init=rnorm(n-1); #random initialization of opt. pool weights
  for(i in 1:T){
    #print(paste(i,'of',T))
    w_pool = m.dir(alpha);
    w_BMA = equal_weight;
    wmat=cbind(diag(n),equal_weight);
    
    if(i>h){
      F = 1:(i-h);
      pp=p[F,];
      nll=function(par){
        ww=w.fun(par);
        -sum(log(pp%*%ww))-z.ddir.unorm(ww,alpha)*SHRINK
      }
      o1=optim(rep(0,n-1),nll,method="BFGS",hessian=TRUE,
               control=list(reltol=reltol));
      
      w_pool=w.fun(o1$par);
      w_BMA=BMA_wheights(pp,w_BMA)
      
    }
    wmat=cbind(wmat,w_pool,w_BMA);
    log_score[i,] = log( p[i,]%*%wmat );
    weights_pool[i,]=w_pool;
    weights_BMA[i,]=w_BMA;
  }
  t_log_score[,h]=colSums(log_score);
  woutfile=paste(run_code,'w_pool',h,'.csv',sep='');
  write.csv(weights_pool,woutfile,row.names=FALSE);
  woutfile=paste(run_code,'w_BMA',h,'.csv',sep='');
  write.csv(weights_BMA,woutfile,row.names=FALSE);
    
  lsoutfile=paste(run_code,'ls_',h,'.csv',sep='');
  write.csv(log_score,lsoutfile,row.names=FALSE);
}
tlsoutfile=paste(run_code,'tls.csv',sep='');
write.csv(t_log_score,tlsoutfile,row.names=FALSE)
apply(t_log_score,2,function(x){which.max(x)})


