functions{
  real zmin_ode(real t, 
                vector z0,
                real lambda_minus,
                real lambda_plus,
                real omega_minus,
                real omega_plus) { 
    real delta = (lambda_minus - lambda_plus + omega_minus - omega_plus)^2 + 4*omega_minus*omega_plus;
    real c1 = (0.5 - (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] + omega_minus/sqrt(delta)*z0[2];
    real c2 = (0.5 + (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] - omega_minus/sqrt(delta)*z0[2];
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
    real c1 = (0.5 - (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] + omega_minus/sqrt(delta)*z0[2];
    real c2 = (0.5 + (lambda_plus - lambda_minus + omega_plus - omega_minus)/(2*sqrt(delta)))*z0[1] - omega_minus/sqrt(delta)*z0[2];
    real zplus = exp((lambda_minus + lambda_plus - omega_minus - omega_plus)/2*t)/(2*omega_minus)*((lambda_plus-lambda_minus+omega_plus-omega_minus)*(c1*exp(sqrt(delta)/2*t) + c2*exp(-sqrt(delta)/2*t))+sqrt(delta)*(c1*exp(sqrt(delta)/2*t) - c2*exp(-sqrt(delta)/2*t)));
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
  real<lower=0, upper=0.01> effomega_minus;
  real<lower=0, upper=0.01> effomega_plus;
}

transformed parameters {
  real<lower=0> omega_minus = effomega_minus*lambda_plus;
  real<lower=0> omega_plus = effomega_plus*lambda_minus;
}

model {
  // Priors
  lambda_minus ~ gamma(6,1.0/2.5);
  lambda_plus ~ gamma(6,1.0/2.5);
  effomega_minus ~ lognormal(-5,0.5);
  effomega_plus ~ lognormal(-5,0.5);
  
  // Likelihood
  for (i in 1:n_times){
    zminus[i] ~ poisson(zmin_ode(i,z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
    zplus[i] ~ poisson(zplus_ode(i,z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
  }
}

generated quantities {
  real pred_minus[n_times];
  real pred_plus[n_times];
  
  for (i in 1:n_times){
    pred_minus[i] = zmin_ode(i,z0,lambda_minus,lambda_plus,omega_minus,omega_plus);
    pred_plus[i] = zplus_ode(i,z0,lambda_minus,lambda_plus,omega_minus,omega_plus);
  }
}