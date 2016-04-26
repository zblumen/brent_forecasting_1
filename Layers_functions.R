Layered_scores=function(p,nl=1,h=1,alpha_def=NULL){
  if(is.null(alpha_def)){
    alpha_def=list();
    for(i in 1:nl){
      alpha_def[[i]]=list();
      alpha_def[[i]]$alpha=NULL;
      alpha_def[[i]]$def.w=1;  
      alpha_def[[i]]$control=list(reltol=1e-32);
    }
  }
  n=ncol(p);
  TT=nrow(p);
  lp=log(p);
  
  dw_scores=matrix(NA,TT,nl);
  W=list()
  for(i in 1:nl){   ##loop throgh layers
    W[[i]]=list();
    row.st=(i-1)*h+1;  
    col.end=n+(i-1)*3;
    pre_dw_lp = lp[row.st:TT,1:col.end];
    print("pre m1")
    m1=run_layer_pre_dw( pre_dw_lp,h,
                         alpha_def[[i]]$def.w,alpha_def[[i]]$alpha,
                         alpha_def[[i]]$control)
    W[[i]]$pre=m1$w;
    print("m1")
    m2=run_layer_dw(pre_dw_lp);
    W[[i]]$dw=m2$w;
    print("m2")
    post_dw_lp=lp[row.st:TT,1:col.end]
    if(i>1){
      post_dw_lp=cbind(post_dw_lp,dw_scores[row.st:TT,1:(i-1)]);
    }
    
    m3=run_layer_post_dw( post_dw_lp,h,
                          alpha_def[[i]]$def.w,alpha_def[[i]]$alpha);
    W[[i]]$post=m3$w;
    print("m3")
    new_lp=as.matrix(cbind(m1$lp,m3$lp));
    new_dw=m2$lp;
    if(i>1){
      print("begin NA merge")
      na.mat = matrix(NA,(i-1)*h,3);
      new_lp=rbind(na.mat,new_lp);
      print("second NA merge")
      new_dw=c(na.mat[,1],new_dw)
    }
    lp=cbind(lp,new_lp)
    dw_scores[,i]=new_dw;
    print("post final merge")
  }
  names.vec=paste("m",1:n,sep='_');
  for(i in 1:nl){
    names.vec=c(names.vec,paste(c("bma","opt","ms"),i,sep='_'));
  }
  names.vec=c(names.vec,paste("dw",1:nl,sep='_'));
  print("name.vec created")
  log_scores=as.data.frame(cbind(lp,dw_scores)); 
  print(names(log_scores))
  names(log_scores)=names.vec;
  return(list(lp=log_scores,W=W))
}
###############################################
###############################################
run_layer_pre_dw=function(lp,h,def.w,alpha,control){
  n=nrow(lp);
  m=ncol(lp);
  log_scores=matrix(NA,n,2); 
  w=list();
  log_scores=as.data.frame(log_scores);
  names(log_scores)=c("bma","opt");
  BMA = bma(lp,h,def.w);
  log_scores$bma=BMA$lp;
  w[["bma"]]=BMA$w;
  
  opt = opt_pool(lp,h,def.w,alpha,control)
  log_scores$opt=opt$lp;
  w[["opt"]]=opt$w;
  return(list(lp=log_scores,w=w))
}
###############################################
###############################################
run_layer_dw=function(lp,def.w=1){
  lp=as.matrix(lp);
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
  log_scores=matrix(NA,n,1); 
  w=list();
  log_scores=as.data.frame(log_scores);
  names(log_scores)=c("ms");
  
  
  ms = model_selector(lp,h,def.w);
  log_scores$ms=ms$lp;
  w[["ms"]]=ms$ind;
  return(list(lp=log_scores,w=w))
}
###############################################
###############################################
z.ddir.unorm=function(x,a){ #unormalized log dir density
  sum((a-1)*log(x))
}
m.dir=function(a){ #dericlet mode
  (a-1)/(sum(a)-length(a))
}
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
#########################################################
model_selector=function(lp,h,def.w){
  lp=as.matrix(lp);
  n=nrow(lp);
  m=ncol(lp);
  if(length(def.w==1)){
    def.w=rep(def.w,m)/(m*def.w)
  }else if(length(def.w)!=m){
    stop("def.w must be the same length as ncol(lp)")
  }
  cum_lp=apply(lp,2,cumsum);
  
  out_lp=numeric(n);
  ind=numeric(n)
  for(i in 1:n){
    if(i>h){
      ind[i]=which.max(cum_lp[i-h,])
      out_lp[i]=lp[i,ind[i]]
    }else{
      ind[i]=0;
      out_lp[i] = crossprod(lp[i,],def.w)
    }
  }
  return(list(lp=out_lp,ind=ind))
}
###############################################
###############################################
bma=function(lp,h,def.w=1,alpha=NULL){
  lp=as.matrix(lp);
  n=nrow(lp);
  m=ncol(lp);
  if(is.null(alpha)){
    SHRINK=FALSE;
    if(length(def.w)==1){
      def.w=rep(def.w,m)/(m*def.w)
    }else if(length(def.w)!=m){
      stop("def.w must be the same length as ncol(lp)")
    }
  }else{
    SHRINK=TRUE;
    if(any(alpha<=1)){
      stop("all alpha must be greatrer than 1 for mode to exist")
    }else if(length(alpha)==1){
      alpha=rep(alpha,m)
    }else if(length(alpha)!=m){
      stop("alpha must be the same length as ncol(lp)")
    }else{
      stop("unrecognized error")
    }
    def.w=m.dir(alpha);
  }
  #################################################
  cum_lp=apply(lp,2,cumsum);
  if(SHRINK){
    log_marginal=matrix(NA,n,m);
    for(i in 1:m){
      log_marginal[,i]=cum_lp[,i]+log(def.w[i])
    }    
  }else{
    log_marginal=cum_lp;
  }
  out_lp=numeric(n);
  w=matrix(NA,n,m);
  post_prob=matrix(NA,n,m);
  for(i in 1:n){
    for(j in 1:m){
      #to avoid underflow/overflow calculate in logs
      post_prob[i,j]=1/(sum(exp(log_marginal[i,]-log_marginal[i,j])))      
    }
    if(i>h){
      w[i,]=post_prob[i-h,]
    }else{w[i,]=def.w}
    out_lp[i]=log(crossprod(w[i,],exp(lp[i,])))
  }
  return(list(lp=out_lp,w=w,post_prob=post_prob) )
}

############################################
############################################

opt_pool=function(lp,h,def.w=1,alpha=NULL,control=list(reltol=1e-32)){
  lp=as.matrix(lp);
  n=nrow(lp);
  m=ncol(lp);
  if(is.null(alpha)){
    SHRINK=FALSE;
    if(length(def.w)==1){
      def.w=rep(def.w,m)/(m*def.w)
    }else if(length(def.w)!=m){
      stop("def.w must be the same length as ncol(lp)")
    }
  }else{
    SHRINK=TRUE;
    if(any(alpha<=1)){
      stop("all alpha must be greatrer than 1 for mode to exist")
    }else if(length(alpha)==1){
      alpha=rep(alpha,m)
    }else if(length(alpha)!=m){
      stop("alpha must be the same length as ncol(lp)")
    }else{
      stop("Logic Fail:  unrecognized error")
    }
    def.w=m.dir(alpha);
  }
  #################################################
  p=exp(lp);
  log_score=numeric(n);
  w=matrix(NA,n,m)
  for(i in 1:n){
    if(i>h){
      pp=p[1:(i-h),];
      if(SHRINK){
        nll=function(par){
          ww=w.fun(par);
          -sum(log(pp%*%ww))-z.ddir.unorm(ww,alpha)
        }        
      }else{
        nll=function(par){
          ww=w.fun(par);
          -sum(log(pp%*%ww))
        }
      }
      o1=optim(rep(0,m-1),nll,method="BFGS",hessian=TRUE,
               control=control);
      w[i,]=w.fun(o1$par);
    }else{
      w[i,]=def.w;
    }
    log_score[i] = log( crossprod(p[i,],w[i,]));
  }
  return(list(lp=log_score,w=w))
}