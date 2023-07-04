library(deSolve)
library(scatterplot3d)

# EXAMPLE
Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a*X + Y*Z
    dY <- b*(Y-Z)
    dZ <- -X*Y + c*Y - Z
    list(c(dX, dY, dZ))
  })
}

parameters <- c(a = -8/3, b = -10, c = 28)
state <- c(X = 1, Y = 1, Z = 1)
times <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

plot(out)

scatterplot3d(out[,-1], type = "l")

# MY ODE
# rho, sigma
switching_process1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    drho <- rho*rho*(lambda_plus + omega_min - lambda_min - omega_plus) + rho*(lambda_min - lambda_plus - 2*omega_min) + omega_min
    dsigma <- sigma*((lambda_min + omega_plus - lambda_plus - omega_min)*rho + lambda_plus + omega_min)
    list(c(drho, dsigma))
  })
}

parameters <- c(lambda_min = 15, lambda_plus = 10, omega_min = 0.01, omega_plus = 0.1)

state <- c(rho = 1, sigma = 1)
times <- seq(0, 2, by = 0.001)
out <- ode(y = state, times = times, func = switching_process1, parms = parameters)
plot(out)

switching_process1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    drho <- (lambda_min - omega_plus)*rho + omega_min*(1 - rho) - rho*(lambda_min*rho + lambda_plus*(1 - rho))
    dsigma <- sigma*(lambda_min*rho + lambda_plus*(1 - rho))
    list(c(drho, dsigma))
  })
}

# x,y
switching_process2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- lambda_min*X + omega_min*Y
    dY <- lambda_plus*Y + omega_plus*X
    list(c(dX, dY))
  })
}
state <- c(X = 1000, Y = 100)
times <- seq(0, 4., by = 0.01)
parameters1 <- c(lambda_min = 1.5, lambda_plus = 1.0, omega_min = 0.01, omega_plus = 0.05)
parameters2 <- c(lambda_min = 1.5, lambda_plus = 0.768, omega_min = 0.02, omega_plus = 0.072)
out2 <- ode(y = state, times = times, func = switching_process2, parms = parameters1)
out2b <- ode(y = state, times = times, func = switching_process2, parms = parameters2)
plot(out2)

out2_df <- data.frame(t = out2[,1], ZM = out2[,2], ZP = out2[,3])
out2b_df <- data.frame(t = out2b[,1], ZM = out2b[,2], ZP = out2b[,3])
ggplot() + 
  #geom_line(data = out2_df, aes(x=t,y=ZM),color="red") + 
  geom_line(data = out2_df, aes(x=t,y=ZP),color="red") +
  #geom_line(data = out2b_df, aes(x=t,y=ZM),color="blue") + 
  geom_line(data = out2b_df, aes(x=t,y=ZP),color="blue") +
  ylab("Z")

n_obs <- rpois(n = length(times),
               lambda = out[,2])
plot(n_obs ~ times, xlab = "Time", ylab = "Z-")
points(times, out[,2], type = "l", lwd=2)

saveRDS(out2_df, file = paste0("./simulations_time/ode_[1.5_1.0_001]_4.0_1000.100.rds"))

# (co)variances
covariances <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    dM1 <- lambda_min*M1 + omega_min*M2
    dM2 <- lambda_plus*M2 + omega_plus*M1
    dV1 <- 2*lambda_min*V1 + (alpha_min+beta_min)*M1 + 2*omega_min*C + omega_min*M2
    dV2 <- 2*lambda_plus*V2 + (alpha_plus+beta_plus)*M2 + 2*omega_plus*C + omega_plus*M1
    dC <- (lambda_min + lambda_plus)*C + omega_plus*V1 + omega_min*V2
    list(c(dM1, dM2, dV1, dV2, dC))
  })
}
parameters_cov <- c(lambda_min = 1.5, lambda_plus = 1.0, omega_min = 0.01, omega_plus = 0.05, alpha_min = 1.5, beta_min = 0, alpha_plus = 1.0, beta_plus = 0)
state <- c(M1 = 1000, M2 = 100, V1 = 0, V2 = 0, C = 0)
out <- ode(y = state, times = seq(0, 4, by = 0.01), func = covariances, parms = parameters_cov)
plot(out)

out_df <- data.frame(t = out[,1], M1 = out[,2], M2 = out[,3], V1 = out[,4], V2 = out[,5], C = out[,6])
out_df <- out_df %>% mutate(D1 = sqrt(V1), D2 = sqrt(V2))

ggplot() + 
  geom_line(data = out_df, aes(x=t,y=M1), color = "red") +
  geom_line(data = out2_df, aes(x=t,y=ZM),color = "blue")
ggsave("./imgs/ODEvsCOV.png", dpi=800)
ggplot() + 
  geom_line(data = out_df, aes(x=t,y=M2), color = "red") +
  geom_line(data = out2_df, aes(x=t,y=ZP),color = "blue")

ggplot(out_df, aes(x = t, y = sqrt(C))) + geom_line()

plotr1 <- out_df %>% mutate(R1 = D1/M1) %>%
  ggplot() +
  geom_line(aes(x = t,y = R1), color = "red")
plotd1 <- out_df %>% ggplot(aes(x = t, y = D1)) + geom_line(color = "red")
plotm1 <- out_df %>% ggplot(aes(x = t, y = M1)) + geom_line(color = "red")
plotd1

plotr2 <- out_df %>% mutate(R2 = D2/M2) %>%
  ggplot() +
  geom_line(aes(x = t,y = R2), color = "blue") 
plotd2 <- out_df %>% ggplot(aes(x = t, y = D2)) + geom_line(color = "blue") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
plotm2 <- out_df %>% ggplot(aes(x = t, y = M2)) + geom_line(color = "blue")
plotd2

plotr1 / plotr2
(plotm1 + plotd1) / (plotm2 + plotd2)

ggplot(out_df) +
  geom_line(aes(x = t, y = M1, color = "mean1"), linewidth = 1) +
  geom_ribbon(aes(x = t, ymin = M1 - D1, ymax = M1 + D1), fill = "red", alpha = 0.2)

ggplot(out_df) +
  geom_line(aes(x = t, y = M2, color = "mean2"), linewidth = 1, color = "blue") +
  geom_ribbon(aes(x = t, ymin = M2 - D2, ymax = M2 + D2), fill = "blue", alpha = 0.2)
  
zmin_ode <- function(t, z0, lambda_minus, lambda_plus, omega_minus, omega_plus){
  delta = (lambda_minus - lambda_plus)^2 + 4*omega_minus*omega_plus
  c1 = ((lambda_minus - lambda_plus + sqrt(delta))*z0[1] + 2*omega_minus*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  c2 = (2*omega_plus*z0[1] + (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  zmin = exp((lambda_minus + lambda_plus)*t/2.)*(c1*(lambda_minus - lambda_plus + sqrt(delta))*exp(sqrt(delta)*t/2.) + c2*2*omega_minus*exp(-sqrt(delta)*t/2.))
  return(zmin)
}

zplus_ode <- function(t, z0, lambda_minus, lambda_plus, omega_minus, omega_plus){
  delta = (lambda_minus - lambda_plus)^2 + 4*omega_minus*omega_plus
  c1 = ((lambda_minus - lambda_plus + sqrt(delta))*z0[1] + 2*omega_minus*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  c2 = (2*omega_plus*z0[1] + (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  zplus = exp((lambda_minus + lambda_plus)*t/2.)*(c1*2*omega_plus*exp(sqrt(delta)*t/2.) + c2*(lambda_plus - lambda_minus - sqrt(delta))*exp(-sqrt(delta)*t/2.))
  return(zplus)
}

t <- seq(0,1., by = 0.01)
z0 <- as.array(c(1000,100))

zmin_ode(0, z0, 1.5, 1.0, .01, .01)

y = lapply(t, zplus_ode, z0=z0, lambda_minus=1.5, lambda_plus=1.0, omega_minus=0.01, omega_plus=0.01) %>% unlist()

plot(t, y)

py_simulation <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.0_005/switching_results_1.5_1.0_avg.csv") %>%
  tibble::as_tibble()
colnames(py_simulation) <- c("step","t","Z-","Z+")


ggplot() +
  geom_point(data = out2b_df, aes(x = t, y = ZP), size = 0.5) +
  #geom_point(data = py_simulation, aes(x = t, y = `Z-`), size = 0.8, color = "red") +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 0.768, omega_minus = 0.02, omega_plus = 0.072))

ggplot() +
  geom_point(data = out2_df, aes(x = t, y = ZP), size = 0.5) +
  #geom_point(data = py_simulation, aes(x = t, y = `Z-`), size = 0.8, color = "red") +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = 0.05))


omega_plus = c(0.0, 0.01, 0.05, 0.5)

ggplot() +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[1]), color = 'red') +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[2]), color = "blue") +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[3]), color = "green") +
  #stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[4])) +
  xlim(0,4)

ggplot() +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[1]), color = 'red') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[1]), color = 'red')
  
ggplot() +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[2]), color = 'blue') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = omega_plus[2]), color = 'blue')

ggplot() +
  #stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = 0.05), color = 'blue') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = 0.05), color = 'blue') +
  #stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 0.768, omega_minus = 0.02, omega_plus = 0.072), color = 'red') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 0.768, omega_minus = 0.02, omega_plus = 0.072), color = 'red') +
  xlim(1,4)

