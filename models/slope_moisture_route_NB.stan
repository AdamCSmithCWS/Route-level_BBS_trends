// Model used in the Horned Grebe covariate example
// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only


//iCAR function
 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }

data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=0> nobservers;
  
  real<lower=0> sd_alpha_prior;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators

   array [ncounts] real pond;              // count-level pond predictors - centered and scaled (number of ponds near route each year)

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]


}


parameters {

  vector[nroutes] beta_raw_space;
  real BETA; 

  vector[nroutes] alpha_raw;
  real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta_space;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts


//climate stuff

  vector[nroutes] rho_raw; // route-level effect of ponds 
  real<lower=0> sdrho;    // sd of pond effects
  real RHO; // mean of pond effect

  
}


model {


  vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  real phi;
  vector[nroutes] rho;

// covariate effect on intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);

   beta = beta_space + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
   obs = sdobs*obs_raw;
   
   
   rho = sdrho*rho_raw + RHO;


   
  for(i in 1:ncounts){

    E[i] =  beta[route[i]] * (year[i]-fixedyear) +
    alpha[route[i]] +
    obs[observer[i]] +
    rho[route[i]]*pond[i] +
    eta*firstyr[i];
  }
  
  

  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ normal(0,0.3); //prior on sd of observer effects
 
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers); //soft sum-zero constraint

  sdrho ~ gamma(3,30);//prior on sd of covariate effect
  rho_raw ~ icar_normal(nroutes, node1, node2);// ~ std_normal();//covariate effects
  RHO ~ std_normal();// prior on fixed effect of ponds

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
 
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect
 
  
  //spatial iCAR intercepts and slopes by strata
 // sdalpha ~ gamma(2,2); // alternate zero-avoiding prior on sd of intercept variation with similar limits to normal
 sdalpha ~ normal(0,sd_alpha_prior); //prior on sd of intercept variation

  sdbeta_space ~ gamma(3,30);//zero-avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.3

  beta_raw_space ~ icar_normal(nroutes, node1, node2);
  alpha_raw ~ icar_normal(nroutes, node1, node2);


}

 generated quantities {

  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nroutes] rho;

   // intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);

    beta = beta_space + BETA;
    alpha = (sdalpha*alpha_raw) + ALPHA;
    rho = sdrho*rho_raw + RHO;




 }

