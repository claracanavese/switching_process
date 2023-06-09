functions {
  real[] switching_process(real t,
                          real[] z,
                          real[] theta, //lambda_minus,lambda_plus, omega_minus, omega_plus
                          real[] x_r,
                          int[] x_i
                          ) {
    real dzdt[5];
    real lambda_minus = theta[1];
    real lambda_plus = theta[2];
    real omega_minus = theta[3];
    real omega_plus = theta[4];
    
    dzdt[1] = lambda_minus*z[1] + omega_minus*z[2];
    dzdt[2] = lambda_plus*z[2] + omega_plus*z[1];
    dzdt[3] = 2*lambda_minus*z[3] + lambda_minus*z[1] + 2*omega_minus*z[5] + omega_minus*z[2];
    dzdt[4] = 2*lambda_plus*z[4] + lambda_plus*z[2] + 2*omega_plus*z[5] + omega_plus*z[1];
    dzdt[5] = (lambda_minus + lambda_plus)*z[5] + omega_plus*z[3] + omega_minus*z[4];
    
    return dzdt;
  } 
}

data {
  int<lower=1> n_times; 
  real z0[5];
  real t0;
  int zminus[n_times]; 
  int zplus[n_times]; 
  real t[n_times]; 
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> lambda_minus;
  real<lower=0> lambda_plus;
  real<lower=0,upper=0.02> omega_minus;
  real<lower=0,upper=0.02> omega_plus;
}

transformed parameters {
  real theta[4];
  theta[1] = lambda_minus;
  theta[2] = lambda_plus;
  theta[3] = omega_minus;
  theta[4] = omega_plus;
}

model {
  real z_hat[n_times, 5];
  
  lambda_minus ~ gamma(2.,1.);
  lambda_plus ~ gamma(2.,1.);
  omega_minus ~  cauchy(0.01,0.01);
  omega_plus ~ cauchy(0.01,0.01);
  
  z_hat = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);
  
  for (i in 1:n_times) {
    zminus[i] ~ normal(z_hat[i,1], sqrt(z_hat[i,3]));
    zplus[i] ~ normal(z_hat[i,2], sqrt(z_hat[i,4]));
  }
  
}

generated quantities {
  real pred_minus[n_times];
  real pred_plus[n_times];
  real pred[n_times, 5] = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);
  
  for (i in 1:n_times){
    pred_minus[i] = pred[i,1];
    pred_plus[i] = pred[i,2];
  }
  
  // for (i in 1:n_times){
  //   pred_minus[i] = zmin_ode(t[i],z0,theta[1],theta[2],theta[3],theta[4]);
  //   pred_plus[i] = zplus_ode(t[i],z0,theta[1],theta[2],theta[3],theta[4]);
  // }
}
