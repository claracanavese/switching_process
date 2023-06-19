

// The input data is a vector 'z' of length 'N'.
data {
  int<lower=1> n_times; // number of time steps 
  int z[n_times]; // observed data
  real t[n_times]; // time steps
}

// The parameters accepted by the model
parameters {
  real lambda;
  real<lower=0> alpha;
}


model {
  lambda ~ normal(5,1);
  alpha ~ normal(15,5);
  for (i in 1:n_times) {
    z[i] ~ exponential(lambda/alpha*exp(-lambda*t[i]));
  }
}


