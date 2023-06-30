functions{
  real zmin_ode(real t, 
                vector z0,
                real lambda_minus,
                real lambda_plus,
                real omega_minus,
                real omega_plus) { 
    real delta = (lambda_minus - lambda_plus)^2 + 4*omega_minus*omega_plus;
    real c1 = ((lambda_minus - lambda_plus + sqrt(delta))*z0[1] + 2*omega_minus*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus);
    real c2 = (2*omega_plus*z0[1] + (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus);
    real zmin = exp((lambda_minus + lambda_plus)*t/2.)*(c1*(lambda_minus - lambda_plus + sqrt(delta))*exp(sqrt(delta)*t/2.) + c2*2*omega_minus*exp(-sqrt(delta)*t/2.));
    
    return zmin;
  }
  
  real zplus_ode(real t, 
                vector z0,
                real lambda_minus,
                real lambda_plus,
                real omega_minus,
                real omega_plus) { 
    real delta = (lambda_minus - lambda_plus)^2 + 4*omega_minus*omega_plus;
    real c1 = ((lambda_minus - lambda_plus + sqrt(delta))*z0[1] + 2*omega_minus*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus);
    real c2 = (2*omega_plus*z0[1] + (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus);
    real zplus = exp((lambda_minus + lambda_plus)*t/2.)*(c1*2*omega_plus*exp(sqrt(delta)*t/2.) + c2*(lambda_minus - lambda_plus + sqrt(delta))*exp(-sqrt(delta)*t/2.));
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
  // real<lower=0> lambda_minus;
  // real<lower=0> lambda_plus;
  real lambda_minus;
  real lambda_plus;
  real<lower=0,upper=0.1> omega_minus;
  real<lower=0,upper=0.1> omega_plus;
}

// transformed parameters {
//   real<lower=0> omega_minus = effomega_minus*lambda_plus;
//   real<lower=0> omega_plus = effomega_plus*lambda_minus;
// }

model {
  // Priors
  lambda_minus ~ gamma(4.,1./3.);
  lambda_plus ~ gamma(4.,1./3.);
  
  // lambda_minus ~ normal(0, 1);
  // lambda_minus ~ normal(0, 1);
  // omega_minus ~ lognormal(-3.,1.);
  // omega_plus ~ lognormal(-3.,1.);
  // omega_minus ~ student_t(1,0.05,0.05);
  // omega_plus ~ student_t(1,0.05,0.05);
  omega_minus ~ cauchy(0.05,0.1);
  omega_plus ~ cauchy(0.05,0.1);
  
  // Likelihood
  for (i in 1:n_times){
    zminus[i] ~ normal(zmin_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus),zmin_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
    zplus[i] ~ normal(zplus_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus),zplus_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
    // zminus[i] ~ poisson(zmin_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
    // zplus[i] ~ poisson(zplus_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus));
  
  }
}

generated quantities {
  real pred_minus[n_times];
  real pred_plus[n_times];
  real y_prior[n_times];
  
  for (i in 1:n_times){
    pred_minus[i] = zmin_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus);
    pred_plus[i] = zplus_ode(t[i],z0,lambda_minus,lambda_plus,omega_minus,omega_plus);
    y_prior[i] = cauchy_rng(0.05,0.1);
  }
}
