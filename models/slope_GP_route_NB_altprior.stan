// Gaussian Process, route-level trend model with alternative prior on the distance of trend and abundance covariance
// with structures to support cross-validation of next year's observations
// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only

// GP function - courtesy of rethinking  
functions {
  // Function to compute the squared exponential covariance matrix and
  // cholesky deomposition
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        matrix[N,N] L_K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        
        L_K = cholesky_decompose(K);
        return L_K;
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
  real<lower=0> gp_sq_rho_beta; // slope spatial  GP
  real<lower=0> gp_sq_alpha_beta;
  real<lower=0> gp_sq_rho_alpha; // intercepts spatial GP
  real<lower=0> gp_sq_alpha_alpha;

  vector[nroutes] gp_eta_beta;
  vector[nroutes] gp_eta_alpha;
  
}

transformed parameters {
  //this is the aspect of the GP that explodes the memory requirements nroutes^2 *2
  // Calculate the latent Gaussian process function
  matrix[nroutes, nroutes] beta_LK = cov_GPL2(distances, gp_sq_alpha_beta, gp_sq_rho_beta, delta);
  vector[nroutes] beta_space = beta_LK * gp_eta_beta;

  matrix[nroutes, nroutes] alpha_LK = cov_GPL2(distances, gp_sq_alpha_alpha, gp_sq_rho_alpha, delta);
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
  // gp_sq_rho_beta ~ normal(2, 2);
  // gp_sq_rho_alpha ~ normal(2, 2);
  gp_sq_rho_beta ~ gamma(2,2);
  gp_sq_rho_alpha ~ gamma(2,2);
  // Prior for the GP covariance magnitude
  gp_sq_alpha_beta ~ student_t(5,0,1);
  gp_sq_alpha_alpha ~ student_t(5,0,1);
  // Multiplier for non-centred GP parameterisation
  gp_eta_beta ~ normal(0,0.1);
  gp_eta_alpha ~ std_normal();
  
  sdnoise ~ student_t(3,0,1); //prior on scale of inverse squared dispersion of NBinomial distribution phi = 1/sqrt(sdnoise)

  sdobs ~ normal(0,0.3); //prior on sd of observer effects
 
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

 
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect
  
  
   for(i in 1:ncounts){
      real obs = sdobs*obs_raw[observer[i]];

    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs + eta*firstyr[i];
  }

   count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
}



