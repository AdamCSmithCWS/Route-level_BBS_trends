// Model used in the Rufous Hummingbird habitat covariate example
// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, and covariates on the intercepts and trends
// 

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
  int<lower=1> nobservers;
 
  array [ncounts] int<lower=0> count;   // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index =1 if observer's first year on route, 0 otherwise
  array [ncounts] int<lower=1> observer;   // observer indicators
  
  // mean annual habitat suitability on route (centered and scaled)
  array [nroutes] real route_habitat;              
  // mean rate of change in annual habitat suitability on route (centered)
  array [nroutes] real route_habitat_slope;              

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]

  int<lower=0, upper=1> fit_spatial; 
  // conditional: 
  // if 1 then use spatial component to model intercept residual
  // if 0 then use simple random effect

}

parameters {

  vector[nroutes] beta_raw;// non-centered residual trend
  real BETA; // hyperparameter mean residual trend
  vector[nroutes] rho_beta_raw_hab;// non-centered habitat-based trend
  real rho_BETA_hab; // hyperparameter mean habitat-based trend

  vector[nroutes] alpha_raw;// non-centered residual intercept
  real ALPHA; // hyperparameter mean residual intercept
  vector[nroutes] rho_alpha_raw_hab;// non-centered habitat-based intercept
  real rho_ALPHA_hab; // hyperparameter mean habitat-based intercept

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta;    // sd of residual slopes 
  real<lower=0> sdrho_beta_hab;    // sd of habitat-change effect on slopes 
  real<lower=0> sdalpha;    // sd of residual intercepts
  real<lower=0> sdrho_alpha_hab;    // sd of habitat effect on intercepts

  
}

transformed parameters{
  
   vector[nroutes] beta; // full slope
   vector[nroutes] beta_resid; //residual component of slope
   vector[nroutes] beta_hab; //habitat component of slope
   
  vector[nroutes] alpha; // full intercept
  vector[nroutes] alpha_hab; //habitat component intercepts
  vector[nroutes] alpha_resid; // residual component intercepts
  vector[nobservers] obs; // observer effects
  real phi;//dipsersion of negative binomial
  vector[ncounts] E;           // predicted log-scale counts (lambda)


// covariate effect on intercepts and slopes
 
   beta_resid = (sdbeta*beta_raw) + BETA;
   alpha_resid = (sdalpha*alpha_raw) + ALPHA;
   
   for(s in 1:nroutes){
   beta_hab[s] = ((sdrho_beta_hab*rho_beta_raw_hab[s]) + rho_BETA_hab) * (route_habitat_slope[s]);
   alpha_hab[s] = ((sdrho_alpha_hab*rho_alpha_raw_hab[s]) + rho_ALPHA_hab) * (route_habitat[s]);
   }
   
   alpha =  alpha_resid + alpha_hab;
   beta = beta_resid  + beta_hab;
  
   
   obs = sdobs*obs_raw;

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  for(i in 1:ncounts){
    E[i] =  alpha[route[i]] + beta[route[i]] * (year[i]-fixedyear)  + obs[observer[i]] + eta*firstyr[i];
  }
  
  
}

model {




  // habitat-change effects on slope
  rho_beta_raw_hab ~  normal(0,1); 
  sum(rho_beta_raw_hab) ~ normal(0,0.001*nroutes);//soft zero-sum constraint
  // habitat-change effects on slope
  rho_alpha_raw_hab ~ normal(0,1);
  sum(rho_alpha_raw_hab) ~ normal(0,0.001*nroutes);//soft zero-sum constraint

  
    // spatially varying residual trend and intercepts
  beta_raw ~ icar_normal(nroutes, node1, node2); 
  if(fit_spatial){
   alpha_raw ~ icar_normal(nroutes, node1, node2);
  }else{
  alpha_raw ~ normal(0,1);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  }

  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance

  sdobs ~ normal(0,0.3); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers); //soft-zero-sum constraint

  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  rho_BETA_hab ~ normal(0,0.1);// prior on haibtat-change effect on slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  rho_ALPHA_hab ~ normal(0,1);// prior on habitat effect on intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  

  sdalpha ~ normal(0,1); //prior on sd of intercept 
  sdrho_alpha_hab ~ normal(0,1); //prior on sd of habitat effect on intercept 

  sdbeta ~ normal(0,0.1);// prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  sdrho_beta_hab ~ normal(0,0.1);// prior on sd of habitat-change effect on slope 


  // full count likelihood
  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
 


}

 generated quantities {


   array[nroutes,nyears] real<lower=0> nsmooth;
   array[nroutes,nyears] real<lower=0> nsmooth_no_habitat;
    array[nyears] real<lower=0> NSmooth;
    array[nyears] real<lower=0> NSmooth_no_habitat;
    real T; // full end-point trend of slope-based trajectory
    real T_no_habitat; // end-point trend of slope-based trajectory without habitat-component of slope
    real CH; // total change transformation of trend 
    real CH_no_habitat; //total change transformation of trend without habitat-component of slope
    real CH_dif; // difference between CH and CH_no_habitat
    real T_dif; // difference between T and T_no
    real retrans_obs = 0;// 0.5*(sdobs^2);

// intercepts and slopes


 for(y in 1:nyears){

      for(s in 1:nroutes){



      nsmooth[s,y] = exp(alpha[s] + beta[s] * (y-fixedyear)  + retrans_obs);//
      nsmooth_no_habitat[s,y] = exp(alpha[s] + beta_resid[s] * (y-fixedyear)  + retrans_obs);//
        }
  NSmooth[y] = mean(nsmooth[,y]);
  NSmooth_no_habitat[y] = mean(nsmooth_no_habitat[,y]);
  
    }
    T = 100*(((NSmooth[nyears]/NSmooth[1])^(1.0/nyears))-1);
    
    
    T_no_habitat = 100*(((NSmooth_no_habitat[nyears]/NSmooth_no_habitat[1])^(1.0/nyears))-1);
    
        CH = 100*((NSmooth[nyears]/NSmooth[1])-1);
    
    
    CH_no_habitat = 100*((NSmooth_no_habitat[nyears]/NSmooth_no_habitat[1])-1);
    
    CH_dif = CH-CH_no_habitat;
    T_dif = T-T_no_habitat;
    
  }



 

