library(rstan)
library(rstanarm)
N <- 10

time_simulation <- readRDS("./simulations_time/output_[20_20_001][0.8].rds")
time_simulation <- time_simulation %>% filter("Z+" > 1000)
samples = time_simulation[seq(400000,800000, by = 40000),]

simulation_py <- read.csv("./GitHub/switching_process/Gillespy2/switching_results_1.csv") %>%
  tibble::as_tibble()
simulation_py <- simulation_py[,-1]
t_samples = seq(0.7,1.0, by = 0.03) %>% round(., 2)
samples = simulation_py %>% filter(simulation_py$time %in% t_samples)

data_list <- list(
  n_times = nrow(samples),
  t = samples$time,
  zminus = samples$z_minus,
  zplus = samples$z_plus,
  z0 = c(1,0)
)

model <- rstan::stan_model("./GitHub/switching_process/regression.stan")
fit <- rstan::sampling(model, data_list, init = list(list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01)), chains=4, warmup=3000, iter = 6000)
# fit <- rstan::sampling(model, data_list, chains=4, warmup=3000, iter = 6000)
print(fit)

rstan::extract(fit, pars=c("lambda_minus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_plus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_minus"))
bayesplot::mcmc_areas(fit, pars = c("omega_minus"))
bayesplot::mcmc_trace(fit, pars = c("lambda_minus"))

bayesplot::ppc_intervals(
  y = samples$z_minus,
  yrep = rstan::extract(fit, pars = c("pred_minus"))$pred_minus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
)

bayesplot::ppc_intervals(
  y = samples$z_plus,
  yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
)

posterior = as.data.frame(fit)
print(posterior)
x = posterior %>% dplyr::select(starts_with("pred"))




Y <- rexp(6, 1/7)
typeof(Y)
Y
