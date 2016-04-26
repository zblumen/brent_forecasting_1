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
  real<lower=0,upper=1> trans_p[2]; //transitional probabilties for switching
}
transformed parameters { 
  real<lower=0,upper=1> prior_theta[T]; //@ time t
  real<lower=0,upper=1> post_theta[T+1]; //@ time t-1
  real pm[T]; //temp marginal components
  real ps[2]; // temp for log component densities
  real lpr_sig;
  real lpr_mu;
  real lpr_x;
  real lpr_trans_p;
  real lpr_fraction;
  real y_log_prob[T];
  //initialize bernoulli markov chain
  lpr_sig <- cauchy_log(sigma,0.0,5.0)+log(4.0);//2 half cuachys => log(2)+log(2)=log(4)
  lpr_mu <- normal_log(mu,0.0,30.0);
  lpr_x <- normal_log(x,0.0,30.0);
  lpr_trans_p <- beta_log(trans_p,1.0,1.0);
  lpr_fraction <- (lpr_sig+lpr_mu+lpr_x+lpr_trans_p)/T;
  post_theta[1] <- normal_cdf(x,0.0,1.0); //initilize at time t=0
  for(t in 1:T){
    //prior: The forecasting step-the forecasted prob of state
    prior_theta[t] <- trans_p[1]*post_theta[t]+
      trans_p[2]*(1.0-post_theta[t]);
   
    ps[1] <- log(prior_theta[t])
        + normal_log(y[t],mu[1],sigma[1]);
    ps[2] <- log(1.0-prior_theta[t])
        + normal_log(y[t],mu[2],sigma[2]);
        
      pm[t] <- log_sum_exp(ps);
      y_log_prob[t] <- pm[t]+lpr_fraction;
      //posterior: The updating step-the updated prob of state
      post_theta[t+1] <- exp(ps[1]-pm[t]);
  }
}
model {
  increment_log_prob(sum(y_log_prob));
}
generated quantities {
  //declarations
  int<lower=1,upper=2>z; //latent state
  real<lower=0,upper=1> temp_prior;
  real<lower=0,upper=1> prior_theta_out[TT-H+1]; //@ time t
  real<lower=0,upper=1> post_theta_out[TT-H+2]; //@ time t-1
  real pm_out[TT-H+1]; //temp marginal components
  real ps_out[2]; // temp for log component densities
 # real y_log_prob_out[T];
  simplex[2] theta;
  real y_sims[TT-H+1];
  real out_log_prob[TT-H+1]; //1-step ahead predictive likelihoods
  real pred_like[TT-H+1]; //simulated predicitve likelihoods
  real unif_like[TT-H+1]; //simulated cummulative like
  real pll[2];    //vector of log component pred like
  real pll_1[2];    //vector of 1-step ahead log component pred like
  real<upper=0> ull[2];//vector of log component cumm. pred like
  //initialize bernoulli markov chain
  post_theta_out[1] <- post_theta[T+1]; //initialize with last obs
  for(tt in H:TT){
    temp_prior <- post_theta_out[tt-H+1];
    for(h in 1:H){//update temp prior state
      temp_prior <- trans_p[1]*temp_prior+
      trans_p[2]*(1.0-temp_prior);
    }
    //simulate y
    theta[1] <- temp_prior;
    theta[2] <- 1.0 - temp_prior;
    z <- categorical_rng(theta);
    y_sims[tt-H + 1] <- normal_rng(mu[z],sigma[z]);
    //predictive likelihood + unif like
    for(k in 1:2){
      pll_1[k] <- log(theta[k])
        + normal_log(y[T+tt-H+1],mu[k],sigma[k]);
      pll[k] <- log(theta[k])
        + normal_log(y[T+tt],mu[k],sigma[k]);
      ull[k] <- log(theta[k])
        + normal_cdf_log(y[T+tt],mu[k],sigma[k]);
    }
    out_log_prob[tt-H+1] <- log_sum_exp(pll_1);
    pred_like[tt-H+1] <- exp(log_sum_exp(pll));
    unif_like[tt-H+1] <- exp(log_sum_exp(ull));
    //update posterior
    prior_theta_out[tt-H+1] <- trans_p[1]*post_theta_out[tt-H+1]+
      trans_p[2]*(1.0-post_theta_out[tt-H+1]);
      
    ps_out[1] <- log(prior_theta_out[tt-H+1])
        + normal_log(y[T+tt-H+1],mu[1],sigma[1]);
    ps_out[2] <- log(1.0-prior_theta_out[tt-H+1])
        + normal_log(y[T+tt-H+1],mu[2],sigma[2]);
        
    pm_out[tt-H+1] <- log_sum_exp(ps_out);
    post_theta_out[tt-H+2] <- exp(ps_out[1]-pm_out[tt-H+1]);
  }
}