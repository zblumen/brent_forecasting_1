run_layer_pre_dw=function(lp,h,def.w,alpha,control){
  n=nrow(lp);
  m=ncol(lp);
  log_scores=matrix(NA,n,2); 
  w=list();
  log_scores=as.data.drame(log_scores);
  names(log_scores)=c("bma","opt");
  
  BMA = bma(lp,h,def.w);
  log_scores$bma=BMA$lp;
  w[["bma"]]=BMA$w;
  
  opt = opt_pool(lp,h,def.w,alpha,control)
  log_scores$opt=opt$lp;
  w[["opt"]]=opt$w;
  return(lp=log_scores,w=w)
}
###############################################
###############################################
run_layer_dw=function(lp,def.w=1){
  m=ncol(lp);
  if(length(def.w)==1){
    def.w=rep(def.w/m,m)
  }else if(length(def.w)!=m){
    stop("def.w must be equal to ncol(lp)")
  }else{
    def.w=def.w/sum(def.w)    
  }
 return(list(lp=log(exp(lp)%*%def.w),w=def.w))
}
###############################################
###############################################
run_layer_post_dw=function(lp,h,def.w,alpha=NULL,control=NULL){
  n=nrow(lp);
  m=ncol(lp);
  log_scores=matrix(NA,n,2); 
  w=list();
  log_scores=as.data.drame(log_scores);
  names(log_scores)=c("ms");
  
  
  ms = model_selector(lp,h,def.w);
  log_scores$ms=ms$lp;
  w[["ms"]]=ms$ind;
  return(lp=log_scores,w=w)
}

