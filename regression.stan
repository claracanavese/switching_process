functions{
  real zmin_ode(real t, 
                vector z0,
                real lambda_minus,
                real lambda_plus,
                real omega_minus,
                real omega_plus) { 
    real delta = (lambda_minus - lambda_plus + omega_minus - omega_plus)^2 + 4*omega_minus*omega_plus;
    real c1 = (1/2 - (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] + omega_minus/sqrt(delta)*z0[2];
    real c2 = (1/2 + (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] - omega_minus/sqrt(delta)*z0[2];
    real zmin = exp((lambda_minus + lambda_plus - omega_minus - omega_plus)/2*t)*(c1*exp(sqrt(delta)/2*t) + c2*exp(-sqrt(delta)/2*t));
    
    return zmin;
  }
  
  real zplus_ode(real t, 
                vector z0,
                real lambda_minus,
                real lambda_plus,
                real omega_minus,
                real omega_plus) { 
    real delta = (lambda_minus - lambda_plus + omega_minus - omega_plus)^2 + 4*omega_minus*omega_plus;
    real c1 = (1/2 - (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] + omega_minus/sqrt(delta)*z0[2];
    real c2 = (1/2 + (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] - omega_minus/sqrt(delta)*z0[2];
    real zplus = ;
    return zplus;
  }
}

data {
  int<lower=1> n_times; 
  int zminus[n_times]; 
  int zplus[n_times]; 
  real t[n_times]; 
  vector[2] z0;
}

parameters {
  real<lower=0> lambda_minus;
  real<lower=0> lambda_plus;
  real<lower=0> effomega_minus;
  real<lower=0> effomega_plus;
}

transformed parameters {
  real<lower=0> omega_minus = effomega_minus*lambda_plus;
  real<lower=0> omega_plus = effomega_plus*lambda_minus;
}

model {
  // Priors
  lambda_minus ~ 
  lambda_plus ~
  effomega_minus ~ 
  effomega_plus ~
  
  // Likelihood
  for(t in 1:n_times) {
    zminus[t] ~ poisson(zmin_ode(t,z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
    zplus[t] ~ poisson();
  }
  
}

