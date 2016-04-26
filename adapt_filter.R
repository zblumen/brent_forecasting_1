
adapt.filter2=function(sm,mi,n.floor=2000){
  ##sm: stan model
  ##mi: model information (data, forecast horizon, test/training breakpoint,..etc)
  ##n.floor: the minum effective sample size before re-sampling in STAN
  #calculate G and M
  M= mi$TT-mi$h+1;
  G=(mi$sims-mi$burn.in)*mi$chains;
  #define storage objects
  n.eff=numeric(M);
  p=numeric(M);
  u=numeric(M);
  aw.mat=matrix(NA,G,M);
  y.sims=matrix(NA,G,M);
  joint.like.mat=matrix(NA,G,M);
  joint.pl=numeric(M);
  #Estimate STAN
  fit=sampling(sm,data=list(H=mi$h,T=mi$T,TT=mi$TT,y=mi$y),
                 iter=mi$sims,warmup=mi$burn.in,
                 chains = mi$chains,cores=mi$cores, 
                 control=list(adapt_delta = mi$adapt_delta));
  ext=extract(fit,c('out_log_prob','pred_like','unif_like','y_sims'));
  ext$out_prob=exp(ext$out_log_prob);
  non_updated_p = colMeans(ext$pred_like);
  non_updated_u = colMeans(ext$unif_like);
  #estimate first iteration pred_like and other
  #3 calculate pred_like and unif_like
  p[1]=mean(ext$pred_like[,1]);
  u[1]=mean(ext$unif_like[,1]);
  #4 sample y_sims
  y.sims[,1]=ext$y_sims[,1]
  #5 reset other values
  aw.mat[,1]=rep(1/G,G);
  n.eff[1]=G;
  joint.like.mat[,1]=ext$out_prob[,1];
  joint.pl[1]=mean(joint.like.mat[,1]);
  
  i=1;
  benchmark=10;
  step.size=10;
  for(m in 2:M){
    i=i+1; #this iterator will reset with resampling...hack carefully!!!
    #take a particle filter step 
    #1: calculate importance weights
    aw.mat[,m] = joint.like.mat[,m-1]/(G*joint.pl[m-1])
    #2: collect the n.eff (effective sample size)
    n.eff[m]=1/sum(aw.mat[,m]^2);
    #3: is n.eff higher than n.floor?
    if(n.eff[m]>n.floor){
      #1 calculate pred_like and unif_like
      p[m]=crossprod(aw.mat[,m],ext$pred_like[,i]);
      u[m]=crossprod(aw.mat[,m],ext$unif_like[,i]);
      #2 sample y_sims
      y.sims[,m]=ext$y_sims[sample(1:G,G,replace=TRUE,prob=aw.mat[,m]),i]
      #3 calculate next joint mat
      joint.like.mat[,m]=joint.like.mat[,m-1]*ext$out_prob[,i]
      joint.pl[m]=mean(joint.like.mat[,m])
    }else{
      #1 Reset indicies
      i=1;
      ##################################
      new.T = mi$T+m-1;
      new.TT = mi$TT-m+1;
      ####################################
      #2 run new stan model
      fit=sampling(sm,data=list(H=mi$h,T=new.T,TT=new.TT,y=mi$y),
                   iter=mi$sims,warmup=mi$burn.in,
                   chains = mi$chains,cores=mi$cores, 
                   control=list(adapt_delta = mi$adapt_delta));
      ext=extract(fit,c('out_log_prob','pred_like','unif_like','y_sims'));
      ext$out_prob=exp(ext$out_log_prob);
      #3 calculate pred_like and unif_like
      p[m]=mean(ext$pred_like[,i]);
      u[m]=mean(ext$unif_like[,i]);
      #4 sample y_sims
      y.sims[,m]=ext$y_sims[,i]
      #5 reset other values
      aw.mat[,m]=rep(1/G,G);
      n.eff[m]=G;
      joint.like.mat[,m]=ext$out_prob[,i];
      joint.pl[m]=mean(joint.like.mat[,m]);
    }
    perc.done=100*m/M;
    if(perc.done>=benchmark){
      message(paste("adaptive filter ",round(perc.done),"% complete",sep=''));
      benchmark=benchmark+step.size;
    }
  }
  return(list(p=p,u=u,n.eff=n.eff,y=y.sims,w=aw.mat,
              non_updated_p=non_updated_p,
              non_updated_u=non_updated_u ))
}
#################################################################
#################################################################
