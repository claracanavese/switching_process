library(rstan)

N <- 10

time_simulation <- readRDS("./simulations_time/output_[20_20_01][0.694].rds")
time_simulation <- filter(time_simulation, `Z+` > 1000)
samples <- time_simulation[sample(nrow(time_simulation),N),]

data_list <- list(
  N = N,
  times = samples$t,
  Zminus = samples$`Z-`,
  Zplus = samples$`Z+`
)

model <- rstan::stan_model("pippo.stan")
fit <- rstan::sampling(model, data_list, init = list(list(lambda_minus = 20, alpha_minus=20, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 20, alpha_minus=20, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 20, alpha_minus=20, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 20, alpha_minus=20, omega_plus = 0.1, omega_minus = 0.01)), chains=4, warmup=3000, iter = 6000)
print(fit)

rstan::extract(fit, pars=c("lambda_minus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_minus"))

bayesplot::mcmc_trace(fit)

Y <- rexp(6, 1/7)
typeof(Y)
Y
