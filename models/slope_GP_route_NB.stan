// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only

// GP function - courtesy of https://rpubs.com/NickClark47/stan_geostatistical 
functions {
  // Function to compute the squared exponential covariance matrix and
  // take the cholesky decomposition for more efficient computation
  // of multivariate normal realisations
  matrix L_cov_exp_quad_dis(
          int N,
          matrix distances,
                real alpha,
                real rho,
                real delta) {
    matrix[N, N] K;
    K = square(alpha) * exp(-0.5 * (square(distances / rho))) +
            diag_matrix(rep_vector(delta, N));
    return cholesky_decompose(K); //this explodes the computational requirements for large N (n^3)
  }
}
data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
 // spatial distance matrix information
 
  matrix[nroutes, nroutes] distances;   // distance matrix (in km/1000)


}

transformed data {
  // Small offset to ensure the covariance matrix is positive definite
  real delta = 1e-9;
}

parameters {

  real BETA; 

  real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  //real<lower=0> sdbeta_space;    // sd of slopes 
  real<lower=0> gp_rho_beta; // slope spatial  GP
  real<lower=0> gp_alpha_beta;
  real<lower=0> gp_rho_alpha; // intercepts spatial GP
  real<lower=0> gp_alpha_alpha;

  vector[nroutes] gp_eta_beta;
  vector[nroutes] gp_eta_alpha;
  
}

transformed parameters {
  //this is the aspect of the GP that explodes the memory requirements nroutes^2 *2
  // Calculate the latent Gaussian process function
  matrix[nroutes, nroutes] beta_LK = L_cov_exp_quad_dis(nroutes, distances, gp_alpha_beta, gp_rho_beta, delta);
  vector[nroutes] beta_space = beta_LK * gp_eta_beta;

  matrix[nroutes, nroutes] alpha_LK = L_cov_exp_quad_dis(nroutes, distances, gp_alpha_alpha, gp_rho_alpha, delta);
  vector[nroutes] alpha_space = alpha_LK * gp_eta_alpha;

   vector[nroutes] beta = beta_space + BETA;
   vector[nroutes] alpha = alpha_space + ALPHA;
   real phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations


}

model {


  vector[ncounts] E;           // log_scale additive likelihood

  // Prior for the GP length scale (i.e. spatial decay)
  // very small values will be inidentifiable, so an informative prior
  // is a must
  // gp_rho_beta ~ normal(2, 2.5);
  // gp_rho_alpha ~ normal(2, 2.5);
  gp_rho_beta ~ inv_gamma(5,5);
  gp_rho_alpha ~ inv_gamma(5,5);
  // Prior for the GP covariance magnitude
  gp_alpha_beta ~ student_t(5,0,3);
  gp_alpha_alpha ~ student_t(5,0,3);
  // Multiplier for non-centred GP parameterisation
  gp_eta_beta ~ normal(0,0.1);
  gp_eta_alpha ~ student_t(5,0,3);
  
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance

  sdobs ~ std_normal(); //prior on sd of observer effects
 
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

 
  BETA ~ normal(0,0.2);// prior on fixed effect mean slope
  ALPHA ~ student_t(10,0,3);;// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect
  
  
   for(i in 1:ncounts){
      real obs = sdobs*obs_raw[observer[i]];

    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs + eta*firstyr[i];
  }

   count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
}



