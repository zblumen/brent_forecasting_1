data {
  int<lower=1> H; //forecast horizon
  int<lower=H> T; // number of observations in-sample
  int<lower=H> TT; //number of observations out-of-sample
  real y[T + TT]; // observations
}
parameters {
  real mu;
  real<lower=0> sigma; // scales of mixture components
  real<lower=0,upper=1> ar_trans;
}
transformed parameters { 
  real<lower=-1,upper=1> ar;
  real means[T - 1];
  real pm[T - 1]; //log likelihood component
  real lpr_sig;
  real lpr_mu;
  real lpr_ar;
  
  lpr_sig <- cauchy_log(sigma,0.0,5.0)+log(2.0);
  lpr_mu <- normal_log(mu,0.0,5.0);
  lpr_ar <- beta_log(ar_trans,1.0,1.0);
  ar <- 2.0*ar_trans - 1.0; //ensures stationarity; ar \in (-1,1)
  
  for(t in 2:T){ //I use the unconditional likelihood
    means[t - 1] <- mu + ar*y[t - 1];
    pm[t - 1] <- normal_log(y[t],means[t - 1],sigma);
  }
}
model {
  increment_log_prob(sum(pm)+lpr_sig+lpr_mu+lpr_ar);
}
generated quantities {
  real temp_mean;
  real out_log_prob[TT-H+1]; //1-step ahead predictive likelihoods
  real y_sims[TT-H+1]; //simulated predicictions
  real pred_like[TT-H+1]; //H-step ahead predicitve likelihoods
  real unif_like[TT-H+1]; //H-step ahead PIT's
  //pred like
  for(tt in H:TT){
    y_sims[tt-H + 1] <- y[T+tt-H];//start with contemporanius time period
    for(h in 1:H){ //forecast out H-steps
      //forecast mean
      temp_mean <- mu + ar*y_sims[tt-H + 1];
      //forecast y
      y_sims[tt-H + 1] <- normal_rng(temp_mean,sigma);
    }//this loop is inneficient...but I can live with it
    
    //predictive likelihood and PIT 
    out_log_prob[tt-H+1] <- normal_log(y[T+tt-H+1],mu +ar*y[T+tt-H],sigma);
    pred_like[tt-H+1] <- exp(normal_log(y[T+tt],temp_mean,sigma));
    unif_like[tt-H+1] <- normal_cdf(y[T+tt],temp_mean,sigma);   
  }
}