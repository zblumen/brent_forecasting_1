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
  real<lower=0,upper=1> ar_trans;
}
transformed parameters { //to keep the means from shifting between components
  simplex[2] theta; // mixing proportions plus origninal
  real ar;
  real<lower=0,upper=1> post_theta[T-1];
  matrix[T-1,2] means;
  real pm[T-1]; //log-likrlihood components
  real ps[2]; // temp for log component densities
  real lpr_sig;
  real lpr_mu;
  real lpr_x;
  real lpr_ar;
  real lpr_fraction; //log-prior components
  real y_log_prob[T-1]; //log-target component
  theta[1] <- normal_cdf(x,0.0,1.0); #/3.0 +2.0/3.0; //first model must have greatest weight
  theta[2] <- 1.0-theta[1];
  ar <- 2.0*ar_trans- 1.0;
  lpr_sig <- cauchy_log(sigma,0.0,5.0)+log(4.0);
  lpr_mu <- normal_log(mu,0.0,5.0);
  lpr_x <- normal_log(x,0.0,30.0);
  lpr_ar <- beta_log(ar_trans,1.0,1.0);
  lpr_fraction <- (lpr_sig+lpr_mu+lpr_x+lpr_ar)/(T-1);
  for(t in 2:T){
    for(k in 1:2){
      means[t-1,k] <- ar*y[t-1] + mu[k];
      ps[k] <- log(theta[k])
        + normal_log(y[t],means[t-1,k],sigma[k]);
    }
      pm[t-1] <- log_sum_exp(ps);
      post_theta[t-1] <- exp(ps[1]-pm[t-1]);
      y_log_prob[t-1] <- pm[t-1]+lpr_fraction;
  }

}
model {
  increment_log_prob(sum(y_log_prob));
}
generated quantities {
  real out_log_prob[TT-H+1]; //1-step ahead predictive likelihoods
  real temp_mean;
  real temp_sim;
  real y_sims[TT-H+1]; //simulated predicictions
  real pred_like[TT-H+1]; //simulated predicitve likelihoods
  real unif_like[TT-H+1]; //simulated cummulative like
  int<lower=1,upper=2> z; //latent model selector
  real pll[2];    //vector of log component pred like
  real pll_1[2];    //vector of log component pred like
  real<upper=0> ull[2];//vector of log component cumm. pred like
  
  for(tt in H:TT){
  //simulate future values
  y_sims[tt-H + 1] <- y[T+tt-H];//start with contemporanius time period
    for(h in 1:H){ //forecast out H-steps
      //forecqst state
      z <- categorical_rng(theta);
      //forecast mean
      temp_mean <- mu[z] + ar*y_sims[tt-H + 1];
      //forecast y
      y_sims[tt-H + 1] <- normal_rng(temp_mean,sigma[z]);
    }//this loop is inneficient...but I can live with it
    temp_sim <- temp_mean  - mu[z]; //move backward one sim
  //pred like
    for(k in 1:2){
      temp_mean <- mu[k] + temp_sim;
      pll_1[k] <- log(theta[k])
        + normal_log(y[T+tt-H+1],mu[z] + ar*y[T+tt-H],sigma[k]);
      pll[k] <- log(theta[k])
        + normal_log(y[T+tt],temp_mean,sigma[k]);
      ull[k] <- log(theta[k])
        + normal_cdf_log(y[T+tt],temp_mean,sigma[k]);
    }
    out_log_prob[tt-H+1] <- log_sum_exp(pll_1);   
    pred_like[tt-H+1] <- exp(log_sum_exp(pll));
    unif_like[tt-H+1] <- exp(log_sum_exp(ull));
  }
}