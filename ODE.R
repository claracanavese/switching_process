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
    drho <- (lambda_min - omega_plus)*rho + omega_min*(1 - rho) - rho*(lambda_min*rho + lambda_plus*(1 - rho))
    dsigma <- sigma*(lambda_min*rho + lambda_plus*(1 - rho))
    list(c(drho, dsigma))
  })
}

parameters <- c(lambda_min = 15, lambda_plus = 10, omega_min = 0.01, omega_plus = 0.1)
state <- c(rho = 1, sigma = 1)
times <- seq(0, 3, by = 0.01)
out <- ode(y = state, times = times, func = switching_process, parms = parameters)
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
    dX <- (lambda_min - omega_plus)*X + omega_min*Y
    dY <- (lambda_plus - omega_min)*Y + omega_plus*X
    list(c(dX, dY))
  })
}
state <- c(X = 1, Y = 0)
times <- seq(0, 1, by = 0.001)
out <- ode(y = state, times = times, func = switching_process2, parms = parameters)
plot(out)

saveRDS(out, file = paste0("./simulations_time/ode_[15_10_01].rds"))
