library(rstan)
library(rstanarm)
N <- 10

time_simulation <- readRDS("./simulations_time/output_[20_20_001][0.8].rds")

time_simulation <- time_simulation %>% filter("Z+" > 1000)

samples = time_simulation[seq(400000,800000, by = 40000),]


# samples <- time_simulation[sample(nrow(time_simulation),N),]

data_list <- list(
  n_times = nrow(samples),
  t = samples$t,
  zminus = samples$`Z-`,
  zplus = samples$`Z+`,
  z0 = c(1,0)
)

model <- rstan::stan_model("./GitHub/switching_process/regression.stan")
fit <- rstan::sampling(model, data_list, init = list(list(lambda_minus = 20, lambda_plus=20, omega_plus = 0.01, omega_minus = 0.1),list(lambda_minus = 20, lambda_plus=20, omega_plus = 0.01, omega_minus = 0.1),list(lambda_minus = 20, lambda_plus=20, omega_plus = 0.01, omega_minus = 0.1),list(lambda_minus = 20, lambda_plus=20, omega_plus = 0.01, omega_minus = 0.1)), chains=4, warmup=3000, iter = 6000)
fit <- rstan::sampling(model, data_list, chains=4, warmup=3000, iter = 6000)
print(fit)



bayesplot::ppc_intervals(
  y = samples$`Z-`,
  yrep = rstan::extract(fit, pars = c("pred_minus"))$pred_minus %>% as.matrix(),
  x = samples$t,
  prob = 0.5
)

bayesplot::ppc_intervals(
  y = samples$`Z+`,
  yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix(),
  x = samples$t,
  prob = 0.5
)

posterior = as.data.frame(fit)
print(posterior)
x = posterior %>% dplyr::select(starts_with("pred"))


rstan::extract(fit, pars=c("lambda_minus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_plus"))
bayesplot::mcmc_areas(fit, pars = c("omega_minus"))

bayesplot::mcmc_trace(fit)

Y <- rexp(6, 1/7)
typeof(Y)
Y
