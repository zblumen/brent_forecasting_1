/*functions {
  real state_to_fut_prob(int z, real t_p){
    real OUT[2];
    if(z==1){
      OUT[1] <- t_p[1];
      OUT[2] <- 1.0-t_p[1];
    }else if(z==2){
      OUT[1] <- t_p[2];
      OUT[2] <- 1.0-t_p[2];      
    }else{
      reject("z can only be equal to 1 or 2")
    }
    return OUT;
  }
}*/
data {
  int<lower=1> H; //forecast horizon
  int<lower=H> T; // number of observations in-sample
  int<lower=H> TT; //number of observations ouyt-of-sample
  real y[T+TT]; // observations
}
parameters {
  real mu[2]; // locations of mixture components
  real<lower=0> sigma[2]; // scales of mixture components
  real x; //latent component
  real<lower=0,upper=1> trans_p[2]; //transitional probabilties for switching (a,(1-b))
  real<lower=0,upper=1> ar_trans; 
}
transformed parameters {
  real ar;
  real<lower=0,upper=1> prior_theta[T-1]; //@ time t
  real<lower=0,upper=1> post_theta[T]; //@ time t-1
  real pm[T-1]; //temp marginal components
  matrix[T-1,2] means;
  real ps[2]; // temp for log component densities
  real lpr_sig;
  real lpr_mu;
  real lpr_x;
  real lpr_trans_p;
  real lpr_ar;
  real lpr_fraction;
  real y_log_prob[T-1];
  
  ar <- 2.0*ar_trans - 1.0;
  //initialize bernoulli markov chain
  lpr_sig <- cauchy_log(sigma,0.0,5.0)+log(4.0);//2 half cuachys => log(2)+log(2)=log(4)
  lpr_mu <- normal_log(mu,0.0,30.0);
  lpr_x <- normal_log(x,0.0,30.0);
  lpr_trans_p <- beta_log(trans_p,1.0,1.0);
  lpr_ar <- beta_log(ar_trans,1.0,1.0);
  lpr_fraction <- (lpr_sig+lpr_mu+lpr_x+lpr_trans_p+lpr_ar)/T;
  post_theta[1] <- normal_cdf(x,0.0,1.0); //initilize at time t=0
  for(t in 2:T){
    //prior: The forecasting step-the forecasted prob of state
    prior_theta[t-1] <- trans_p[1]*post_theta[t-1]+
      trans_p[2]*(1.0-post_theta[t-1]);
    means[t-1,1] <- mu[1]+ar*y[t-1];
    means[t-1,2] <- mu[2]+ar*y[t-1];    
    ps[1] <- log(prior_theta[t-1])
        + normal_log(y[t],means[t-1,1],sigma[1]);
    ps[2] <- log(1.0-prior_theta[t-1])
        + normal_log(y[t],means[t-1,2],sigma[2]);
        
      pm[t-1] <- log_sum_exp(ps);
      y_log_prob[t-1] <- pm[t-1]+lpr_fraction;
      //posterior: The updating step-the updated prob of state
      post_theta[t] <- exp(ps[1]-pm[t-1]);
  }
}
model {
  increment_log_prob(sum(y_log_prob));
}
generated quantities {
  //declarations
  int<lower=1,upper=2> z; //latent state
  real temp_mean;
  real temp_sim;
  real<lower=0,upper=1> prior_theta_out[TT-H+1]; //@ time t
  real<lower=0,upper=1> post_theta_out[TT-H+2]; //@ time t-1
  real pm_out[TT-H+1]; //temp marginal components
  real ps_out[2]; // temp for log component densities
 # real y_log_prob_out[T];
  simplex[2] theta;
  simplex[2] theta_lag;
  real y_sims[TT-H+1];
  real out_log_prob[TT-H+1]; //1-step ahead predictive likelihoods
  real pred_like[TT-H+1]; //simulated predicitve likelihoods
  real unif_like[TT-H+1]; //simulated cummulative like
  real pll[2];    //vector of log component pred like
  real pll_1[2];    //vector of 1-step ahead log component pred like 
  real<upper=0> ull[2];//vector of log component cumm. pred like
  //initialize bernoulli markov chain
  post_theta_out[1] <- post_theta[T]; //initialize with last obs
  for(tt in H:TT){
    prior_theta_out[tt-H+1] <- trans_p[1]*post_theta_out[tt-H+1]+
      trans_p[2]*(1.0-post_theta_out[tt-H+1]);
    theta[1] <- prior_theta_out[tt-H+1];
    theta[2] <- 1.0 - prior_theta_out[tt-H+1];
    y_sims[tt-H+1] <- y[T+tt-H]; //start at contemporanious time point
    for(h in 1:H){//update temp prior state and simulate y
      z <- categorical_rng(theta);
      theta_lag <- theta; 
      temp_mean <- mu[z]+ar*y_sims[tt-H+1];
      y_sims[tt-H+1] <-  normal_rng(temp_mean,sigma[z]);
      if(z==1){
        theta[1] <- trans_p[1];
        theta[2] <- 1.0-trans_p[1];
      }else if(z==2){
        theta[1] <- trans_p[2];
        theta[2] <- 1.0-trans_p[2];  
      }
    }
    temp_sim <- temp_mean - mu[z]; //move back one
    //predictive likelihood + unif like
    for(k in 1:2){
      pll_1[k] <- log(theta[k])
        + normal_log(y[T+tt-H+1],mu[k]+ar*y[T+tt-H],sigma[k]);
      pll[k] <- log(theta_lag[k])
        + normal_log(y[T+tt],mu[k]+temp_sim,sigma[k]);
      ull[k] <- log(theta_lag[k])
        + normal_cdf_log(y[T+tt],mu[k]+temp_sim,sigma[k]);
    }
    out_log_prob[tt-H+1] <- log_sum_exp(pll_1);    
    pred_like[tt-H+1] <- exp(log_sum_exp(pll));
    unif_like[tt-H+1] <- exp(log_sum_exp(ull));
    
    //update posterior
    ps_out[1] <- log(prior_theta_out[tt-H+1])
        + normal_log(y[T+tt-H+1],mu[1],sigma[1]);
    ps_out[2] <- log(1.0-prior_theta_out[tt-H+1])
        + normal_log(y[T+tt-H+1],mu[2],sigma[2]);
        
    pm_out[tt-H+1] <- log_sum_exp(ps_out);
    post_theta_out[tt-H+2] <- exp(ps_out[1]-pm_out[tt-H+1]);
  }
}