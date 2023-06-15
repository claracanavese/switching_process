library(rstan)
library(deSolve)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

N <- 10
simulation <- readRDS("./simulations_time.rds")
samples <- simulation[sample(nrow(simulation),N),]

dlist <- list(
  n_times = nrow(samples), # number of time steps
  z = samples$Z, # observed data
  t = samples$t # timesteps
)

model <- rstan::stan_model("exp_growth.stan")
fit <- rstan::sampling(model, dlist, chains=4, warmup=3000, iter = 6000)
print(fit)


rstan::extract(fit, pars=c("lambda"))
bayesplot::mcmc_areas(fit, pars = c("lambda"))
bayesplot::mcmc_areas(fit, pars = c("alpha"))
bayesplot::mcmc_trace(fit)
