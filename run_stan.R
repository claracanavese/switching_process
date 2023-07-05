library(rstan)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(bayesplot)

time_simulation <- readRDS("./simulations_time/output_[20_20_001][0.8].rds")
time_simulation <- time_simulation %>% filter("Z+" > 1000)
samples = time_simulation[seq(400000,800000, by = 40000),]

simulation_py <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.0_005/switching_results_1.5_1.0_avg.csv") %>% tibble::as_tibble()
# simulation_py <- read.csv("./GitHub/switching_process/Gillespy2/no_switching.csv") %>% tibble::as_tibble()
simulation_py <- simulation_py[,-1]
t_samples = seq(1.0 ,4.0, by = 0.2) %>% round(., 2)
samples = simulation_py %>% filter(simulation_py$time %in% t_samples)
samples = simulation_py[seq(71,121,by=5),]

data_list <- list(
  n_times = nrow(samples),
  t = samples$time,
  zminus = as.integer(samples$z_minus),
  zplus = as.integer(samples$z_plus),
  z0 = c(1000,100)
)

model <- rstan::stan_model("./GitHub/switching_process/regression.stan")
# fit <- rstan::sampling(model, data_list, init = list(list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01),list(lambda_minus = 15, lambda_plus=10, omega_plus = 0.1, omega_minus = 0.01)), chains=4, warmup=3000, iter = 6000)
fit <- rstan::sampling(model, data_list, chains=4, warmup=5000, iter=10000, cores=4)

bayesplot::mcmc_trace(fit, pars = c("lambda_minus","lambda_plus","omega_minus","omega_plus"))
print(fit, pars = c("lambda_minus","lambda_plus","omega_minus","omega_plus"), digits_summary = 3)
print(fit, digits_summary = 1e-8)


rstan::extract(fit, pars=c("omega_plus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_minus"))
bayesplot::mcmc_areas(fit, pars = c("lambda_plus"))
bayesplot::mcmc_areas(fit, pars = c("omega_minus"))
bayesplot::mcmc_areas(fit, pars = c("omega_plus"))

prior_lambda = ggplot() +
  stat_function(fun=dgamma, args = list(shape = 2., rate = 1.)) +
  xlim(0,5) + ggtitle("Prior")
prior_omega = ggplot() +
  stat_function(fun=dgamma, args = list(shape = 2., rate = 1/0.01)) +
  xlim(0.,0.2) + ggtitle("Prior")

ggplot() +
  stat_function(fun=dgamma, args = list(shape = 2., rate = 100.)) +
  xlim(0.,0.2) + ggtitle("Prior")

bayesplot::mcmc_trace(fit, pars = c("lambda_minus"))
bayesplot::mcmc_trace(fit, pars = c("lambda_plus"))
bayesplot::mcmc_trace(fit, pars = c("omega_plus"))
bayesplot::mcmc_trace(fit, pars = c("omega_minus"))

bayesplot::ppc_intervals(
  y = samples$z_minus,
  yrep = rstan::extract(fit, pars = c("pred_minus"))$pred_minus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
) +
  bayesplot::ppc_intervals(
  y = samples$z_plus,
  yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
)

posterior = as.data.frame(fit)

# plot samples from PRIOR
y_prior = posterior %>% dplyr::select(starts_with("y_prior"))
y_prior <- reshape2::melt(y_prior)
ggplot(y_prior) + geom_density(aes(x=value,y=after_stat(density))) + xlim(-0.2,0.2)#+ stat_function(fun=dgamma, args = list(shape = 8.5, rate = 1./1.8))

# POSTERIOR vs PRIOR
posterior_lambda_min = posterior %>% ggplot() + geom_density(aes(x = lambda_minus, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5)
posterior_lambda_plus = posterior %>% ggplot() + geom_density(aes(x = lambda_plus, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5)
posterior_omega_min = posterior %>% ggplot() + geom_density(aes(x = omega_minus, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.2)
posterior_omega_plus = posterior %>% ggplot() + geom_density(aes(x = omega_plus, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.2)

posterior_lambda_min / prior_lambda
posterior_lambda_plus / prior_lambda
posterior_omega_min / prior_omega
posterior_omega_plus / prior_omega

ggplot() +
  geom_density(data = posterior, aes(x = lambda_minus, y = after_stat(density))) +
  stat_function(fun=dgamma, args = list(shape = 8.5, rate = 1.8)) 

posterior <- posterior[,1:4]
mcmc_intervals(posterior,)

x = posterior %>% dplyr::select(starts_with("pred_minus"))
x <- x[nrow(x),]
x <- reshape2::melt(x)
# x[,-1]
min_pred_df <- data.frame(time = t_samples, z_min_pred = x[,-1])
ggplot(min_pred_df, aes(x = time, y = z_min_pred)) + geom_point()

y = posterior %>% dplyr::select(starts_with("pred_plus"))
y <- y[nrow(y),]
y <- reshape2::melt(y)
plus_pred_df <- data.frame(time = t_samples, z_plus_pred = y[,-1])
ggplot(plus_pred_df, aes(x = time, y = z_plus_pred)) + geom_point()

ggplot() +
  geom_line(data = simulation_py, aes(x = time, y = z_minus),  color = "forestgreen") +
  geom_point(data = min_pred_df, aes(x = time, y = z_min_pred)) +
  xlim(0.5,1)

ggplot() +
  geom_line(data = simulation_py, aes(x = time, y = z_plus),  color = "forestgreen") +
  geom_point(data = plus_pred_df, aes(x = time, y = z_plus_pred)) +
  xlim(0.5,1)

# SECOND MODEL
data_list <- list(
  n_times = nrow(samples),
  z0 = c(1000,100,0,0,0),
  t0 = simulation_py$time[1],
  zminus = as.integer(samples$z_minus),
  zplus = as.integer(samples$z_plus),
  t = samples$time
)
model <- rstan::stan_model("./GitHub/switching_process/regressionODE.stan")
fit <- rstan::sampling(model, data_list, chains=4, warmup=5000, iter=10000, cores=4)

bayesplot::mcmc_trace(fit, pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
print(fit, pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"), digits_summary = 3)
print(fit) 

bayesplot::ppc_intervals(
  y = samples$z_minus,
  yrep = rstan::extract(fit, pars = c("pred_minus"))$pred_minus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
) +
  bayesplot::ppc_intervals(
    y = samples$z_plus,
    yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix(),
    x = samples$time,
    prob = 0.5
  )      
  
posterior = as.data.frame(fit)
posterior_lambda_min = posterior %>% ggplot() + geom_density(aes(x = `theta[1]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5)
posterior_lambda_plus = posterior %>% ggplot() + geom_density(aes(x = `theta[2]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5)
posterior_omega_min = posterior %>% ggplot() + geom_density(aes(x = `theta[3]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.2)
posterior_omega_plus = posterior %>% ggplot() + geom_density(aes(x = `theta[4]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.2)

posterior_lambda_min / prior_lambda
posterior_lambda_plus / prior_lambda
posterior_omega_min / prior_omega
posterior_omega_plus / prior_omega

yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix()
View(yrep)
