// non-spatial route-level trends model
// with structures to support cross-validation of next year's observations
// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only



data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
 
   real<lower=0> sd_alpha_prior;


  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
// values for predicting next years data
  int<lower=1> ncounts_pred;
  array[ncounts_pred] int<lower=0> count_pred;              // count observations
  array[ncounts_pred] int<lower=1> route_pred; // route index
  array[ncounts_pred] int<lower=0> firstyr_pred; // first year index
  array[ncounts_pred] int<lower=0> observer_pred; 


}

parameters {

  vector[nroutes] beta_raw_rand;
  real BETA; 

  vector[nroutes] alpha_raw;
  real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
 real<lower=0> sdbeta_rand;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts

  
}


model {


  vector[ncounts] E;           // log_scale additive likelihood
 vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  real phi;

// covariate effect on intercepts and slopes
 
   beta = (sdbeta_rand*beta_raw_rand) + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i];
  }
  
  
  beta_raw_rand ~ normal(0,1);//random slope effects
  sum(beta_raw_rand) ~ normal(0,0.001*nroutes);
  alpha_raw ~ normal(0,1);//random slope effects
  sum(alpha_raw) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ student_t(3,0,1); //prior on scale of inverse squared dispersion of NBinomial distribution phi = 1/sqrt(sdnoise)

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ normal(0,0.3); //prior on sd of observer effects
 
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
 sdalpha ~ normal(0,sd_alpha_prior); //prior on sd of intercept variation

  sdbeta_rand ~ gamma(3,30);//zero-avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.3



}

 generated quantities {

  vector[nroutes] beta;
  vector[nroutes] alpha;
  real phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  vector[ncounts_pred] log_lik;
  vector[ncounts_pred] E_pred;           // log_scale additive likelihood

// intercepts and slopes
 
   
   beta = (sdbeta_rand*beta_raw_rand) + BETA;
    alpha = (sdalpha*alpha_raw) + ALPHA;
  

  
    for(i in 1:ncounts_pred){
    real obs_tmp;
    if(observer_pred[i] == 0) 
    obs_tmp = normal_rng(0,sdobs);
    else
    obs_tmp = sdobs*obs_raw[observer_pred[i]];

    
    E_pred[i] =  beta[route_pred[i]] * ((nyears+1)-fixedyear) + alpha[route_pred[i]] + obs_tmp + eta*firstyr_pred[i];
   log_lik[i] = neg_binomial_2_log_lpmf(count_pred[i] | E_pred[i],phi);
  }
  
 }

