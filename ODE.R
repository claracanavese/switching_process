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
times <- seq(0, 1, by = 0.001)
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
state <- c(X = 1, Y = 0)
times <- seq(0, 1, by = 0.01)
out <- ode(y = state, times = times, func = switching_process2, parms = parameters)
plot(out)

n_obs <- rpois(n = length(times),
               lambda = out[,2])
plot(n_obs ~ times, xlab = "Time", ylab = "Z-")
points(times, out[,2], type = "l", lwd=2)

saveRDS(out, file = paste0("./simulations_time/ode_[15_10_01].rds"))

# (co)variances
covariances <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    dM1 <- (lambda_min - omega_plus)*M1 + omega_min*M2
    dM2 <- (lambda_plus - omega_min)*M2 + omega_plus*M1
    dV1 <- 2*((lambda_min - omega_plus)*V1 + omega_min*C)
    dV2 <- 2*((lambda_plus - omega_min)*V2 + omega_plus*C)
    dC <- (lambda_min + lambda_plus - omega_min - omega_plus)*C + omega_plus*(V1 - M1) + omega_min*(V2 - M2)
    list(c(dM1, dM2, dV1, dV2, dC))
  })
}
state <- c(M1 = 1, M2 = 0, V1 = 0, V2 = 0, C = 0)
out <- ode(y = state, times = times, func = covariances, parms = parameters)
plot(out)