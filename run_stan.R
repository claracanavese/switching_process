library(rstan)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(bayesplot)

simulation_py <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.0_001_005/switching_results_avg.csv") %>% tibble::as_tibble()
simulation_py <- simulation_py[,-1]
t_samples = seq(2.0, 5.0, by = 0.25) %>% round(., 3)
samples = simulation_py %>% filter(simulation_py$time %in% t_samples)

prior_lambda = ggplot() +
  stat_function(fun=dgamma, args = list(shape = 2., rate = 1.)) +
  xlim(0,5) + ggtitle("Prior") + xlab("lambda") + theme(plot.title = element_text(hjust = 0.5)) + ylab("density")
prior_omega = ggplot() +
  stat_function(fun=dgamma, args = list(shape = 2., rate = 1/0.02)) +
  xlim(0.,0.1) + ggtitle("Prior") + xlab("omega") + theme(plot.title = element_text(hjust = 0.5)) + ylab("density")

# NEW MODEL
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

print(fit, pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"), digits_summary = 3)
print(fit) 
saveRDS(fit,"./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/fit.rds")

fit = readRDS("./GitHub/switching_process/Gillespy2/1.5_1.0_001/fit.rds")
print(fit, pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"), digits_summary = 3)

traceplot <- bayesplot::mcmc_trace(fit, pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/mcmc_trace.png", width = 10, height = 7, dpi = 600)

options(scipen = 1)
color_scheme_set("pink")
minuspred = bayesplot::ppc_intervals(
  y = samples$z_minus,
  yrep = rstan::extract(fit, pars = c("pred_minus"))$pred_minus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
) + xlab("t") + ylab("z-") + scale_x_continuous(breaks = pretty) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/zminus_pred.png", width = 10, height = 7, dpi = 600)

pluspred = bayesplot::ppc_intervals(
  y = samples$z_plus,
  yrep = rstan::extract(fit, pars = c("pred_plus"))$pred_plus %>% as.matrix(),
  x = samples$time,
  prob = 0.5
) + xlab("t") + ylab("z+") + scale_x_continuous(breaks = pretty) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))     
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/zplus_pred.png", width = 10, height = 7, dpi = 600)


posterior = as.data.frame(fit)
posterior_lambda_min = posterior %>% ggplot() + geom_density(aes(x = `theta[1]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5) + xlab("lambda_minus") + theme(plot.title = element_text(hjust = 0.5))
posterior_lambda_plus = posterior %>% ggplot() + geom_density(aes(x = `theta[2]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,5) + xlab("lambda_plus") + theme(plot.title = element_text(hjust = 0.5))
posterior_omega_min = posterior %>% ggplot() + geom_density(aes(x = `theta[3]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.1) + xlab("omega_minus") + theme(plot.title = element_text(hjust = 0.5))
posterior_omega_plus = posterior %>% ggplot() + geom_density(aes(x = `theta[4]`, y = after_stat(density))) + ggtitle("Posterior") + xlim(0,0.1) + xlab("omega_plus") + theme(plot.title = element_text(hjust = 0.5))

posterior_lambda_min / prior_lambda
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/lambda_minus_posterior.png", width = 12, height = 7, dpi = 600)
posterior_lambda_plus / prior_lambda
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/lambda_plus_posterior.png", width = 12, height = 7, dpi = 600)
posterior_omega_min / prior_omega
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/omega_minus_posterior.png", width = 12, height = 7, dpi = 600)
posterior_omega_plus / prior_omega
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/gamma_2_100/omega_plus_posterior.png", width = 12, height = 7, dpi = 600)

(minuspred + pluspred) / (posterior_lambda_min + posterior_lambda_plus) / (prior_lambda + prior_lambda) / (posterior_omega_min + posterior_omega_plus) / (prior_omega + prior_omega)
ggsave("./GitHub/switching_process/Gillespy2/1.5_1.0_001_005/gamma_2_50/panel.png", width = 12, height = 14, dpi = 600)



# OLD MODEL
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
  stat_function(fun=dgamma, args = list(shape = 2., rate = 1/0.005)) +
  xlim(0.,0.1) + ggtitle("Prior")

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
