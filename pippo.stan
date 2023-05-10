
functions {
  real erlangm_log (real z, real t, real z0, real t0, real lambda_plus, real lambda_minus, real alpha_minus, real omega_plus, real omega_minus) {
    real prob;
    real rate;
    real lprob;
    rate = lambda_minus/(alpha_minus*cosh(sqrt(omega_minus*omega_plus)*t))*exp(-lambda_minus*t);
    prob = rate^z0*z^(z0-1)*exp(-rate*z)/tgamma(z0);
    lprob = log(lprob);
    return lprob;
  }
  
    real erlangp_log (real z, real t, real z0, real t0, real lambda_plus, real lambda_minus, real alpha_minus, real omega_plus, real omega_minus) {
    real prob;
    real rate;
    real lprob;
    rate = lambda_minus/(alpha_minus*sinh(sqrt(omega_minus*omega_plus)*(t-t0)))*sqrt(omega_minus/omega_plus)*exp(-lambda_minus*(t-t0));
    prob = rate^z0*z^(z0-1)*exp(-rate*z)/tgamma(z0);
    lprob = log(lprob);
    return lprob;
  }
}

data {
  int N;
  int Zplus[N];
  int Zminus[N];
  real times[N];
}

parameters {
  // real lambda_plus;
  real <lower=0> lambda_minus;
  real <lower=0> alpha_minus;
  real <lower=0> omega_plus;
  real <lower=0> omega_minus;
  // real t0;
}

model {
  real lambda_plus;
  lambda_minus ~ gamma(20,1);
  alpha_minus ~ gamma(20,1);
  omega_plus ~ beta(1,10);
  omega_minus ~ beta(1,10);
 
  for (i in 1:N) {
      Zminus[i] ~ exponential(lambda_minus/(alpha_minus*cosh(sqrt(omega_minus*omega_plus)*times[i]))*exp(-lambda_minus*times[i]));
      Zplus[i] ~ exponential(lambda_minus/(alpha_minus*sinh(sqrt(omega_minus*omega_plus)*times[i]))*sqrt(omega_minus/omega_plus)*exp(-lambda_minus*times[i]));  
    }

}
