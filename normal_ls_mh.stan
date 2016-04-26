data {
  int<lower=1> H; //forecast horizon
  int<lower=H> T; // number of observations in-sample
  int<lower=H> TT; //number of observations ouyt-of-sample
  real y[T+TT]; // observations
}
parameters {
  real<lower=0> mu; 
  real<lower=0> sigma; 
}
transformed parameters { 
  real pm[T]; //log-likelihood components
  real lpr_mu;
  real lpr_sig;
  real lpr_fraction;
  real y_log_prob[T]; //log target components
  
  lpr_mu <- normal_log(mu,0.0,5.0);
  lpr_sig <- cauchy_log(sigma,0.0,5.0)+log(2.0);
  lpr_fraction <- (lpr_sig+lpr_mu)/T;
  for(t in 1:T){
    pm[t] <- normal_log(y[t],mu,sigma);
    y_log_prob[t] <- pm[t]+lpr_fraction;
  }
}
model {
  increment_log_prob(sum(y_log_prob));
}
generated quantities {
  //declarations
  real out_log_prob[TT-H+1]; //1-step ahead predictive likelihoods
  real y_sims[TT-H+1]; //simulated predicictions
  real pred_like[TT-H+1]; //simulated predicitve likelihoods
  real unif_like[TT-H+1]; //simulated cummulative like
  //pred like
  for(tt in H:TT){
  y_sims[tt-H+1] <- normal_rng(mu,sigma);
  out_log_prob[tt-H+1] <- normal_log(y[T+tt-H+1],mu,sigma);
  pred_like[tt-H+1] <- exp(normal_log(y[T+tt],mu,sigma));
  unif_like[tt-H+1] <- normal_cdf(y[T+tt],mu,sigma);   
  }
}