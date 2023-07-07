library(deSolve)
library(RColorBrewer)

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
parameters1 <- c(lambda_min = 1.5, lambda_plus = 1.0, omega_min = 0.01, omega_plus = 0.01)
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

state <- c(M1 = 1000, M2 = 100, V1 = 0, V2 = 0, C = 0)
lambda_min = list(0.5,1.0,1.5)
lambda_plus = list(0.5,1.0,1.5)
omega_min = list(0.01,0.05,0.1)
omega_plus = list(0.01,0.05,0.1)

for (lp in lambda_plus) {
  for (op in omega_plus) {
    for (lm in lambda_min) {
      for (om in omega_min) {
        param <- c(lambda_min = lm, lambda_plus = lp, omega_min = om, omega_plus = op, alpha_min = lm, beta_min = 0, alpha_plus = lp, beta_plus = 0)
        out <- paste("out", lm, lp, om, op, sep = "_")
        assign(out,ode(y = state, times = seq(0, 5, by = 0.01), func = covariances, parms = param))
        out_df <- paste("out_df", lm, lp, om, op, sep = "_")
        assign(out_df, data.frame(t = get(out)[,1], M1 = get(out)[,2], M2 = get(out)[,3])) 
      }
    }
  }
}

rm(list = ls())

ggplot() + 
  geom_line(data = out_df_0.5_0.01, aes(x=t,y=M1, color = "1")) +
  geom_line(data = out_df_0.5_0.05, aes(x=t,y=M1, color = "2")) +
  geom_line(data = out_df_1_0.01, aes(x=t,y=M1, color = "3")) +
  geom_line(data = out_df_1_0.05, aes(x=t,y=M1, color = "4")) +
  geom_line(data = out_df_1.5_0.01, aes(x=t,y=M1, color = "5")) +
  geom_line(data = out_df_1.5_0.05, aes(x=t,y=M1, color = "6")) +
  geom_line(data = out_df_0.5_0.1, aes(x=t,y=M1, color = "7")) +
  geom_line(data = out_df_1_0.1, aes(x=t,y=M1, color = "8")) +
  geom_line(data = out_df_1.5_0.1, aes(x=t,y=M1, color = "9")) +
  geom_line(data = out_df_0.5_0.001, aes(x=t,y=M1, color = "10")) +
  geom_line(data = out_df_1_0.001, aes(x=t,y=M1, color = "11")) +
  geom_line(data = out_df_1.5_0.001, aes(x=t,y=M1, color = "12")) +
  xlim(4,5)
  #scale_color_brewer(palette = "Dark2")
ggplot() + 
  geom_line(data = out_df_0.5_0.01, aes(x=t,y=M2, color = "1")) +
  #geom_line(data = out_df_0.5_0.05, aes(x=t,y=M2, color = "2")) +
  #geom_line(data = out_df_1_0.01, aes(x=t,y=M2, color = "3")) +
  geom_line(data = out_df_1_0.05, aes(x=t,y=M2, color = "4")) +
  #geom_line(data = out_df_1.5_0.01, aes(x=t,y=M2, color = "5")) +
  #geom_line(data = out_df_1.5_0.05, aes(x=t,y=M2, color = "6")) +
  geom_line(data = out_df_0.5_0.1, aes(x=t,y=M2, color = "7")) +
  #geom_line(data = out_df_1_0.1, aes(x=t,y=M2, color = "8")) +
  #geom_line(data = out_df_1.5_0.1, aes(x=t,y=M2, color = "9")) +
  #geom_line(data = out_df_0.5_0.001, aes(x=t,y=M2, color = "10")) +
  geom_line(data = out_df_1_0.001, aes(x=t,y=M2, color = "11")) +
  geom_line(data = out_df_1.5_0.001, aes(x=t,y=M2, color = "12")) +
  xlim(4,5)

ggplot() + 
  geom_line(data = out_df_1_0.05, aes(x=t,y=M2, color = "4")) +
  geom_line(data = out_df_0.5_0.1, aes(x=t,y=M2, color = "7")) +
  xlim(4,5)

zplus_ode(5, c(1000,100), 1.5, 1, 0.01,0.05)
zplus_ode(5, c(1000,100), 1.5, 0.5, 0.01,0.1)

mean(out_df_1_0.05[,2] - out_df_0.5_0.1[,2])

mean(out_df_0.5_0.05[,2] - out_df_1_0.05[,2])


out_df <- data.frame()
for (lp in lambda_plus) {
  for (op in omega_plus) {
    for (lm in lambda_min) {
      for (om in omega_min) {
        param <- c(lambda_min = lm, lambda_plus = lp, omega_min = om, omega_plus = op, alpha_min = lm, beta_min = 0, alpha_plus = lp, beta_plus = 0)
        out <- paste("out", lm, lp, om, op, sep = "_")
        assign(out,ode(y = state, times = seq(0, 5, by = 0.01), func = covariances, parms = param))
        if (nrow(out_df) == 0) {
          M1 <- paste("M1")
          out_df <- data.frame(t = get(out)[,1], M1 = get(out)[,2], M2 = get(out)[,3])
          newM1 = paste("M1",lm,lp,om,op,sep="_")
          newM2 = paste("M2",lm,lp,om,op,sep="_")
          names(out_df)[names(out_df) == 'M1'] <- newM1
          names(out_df)[names(out_df) == 'M2'] <- newM2
        }
        else {
          out_df <- cbind(out_df, M1 = get(out)[,2], M2 = get(out)[,3])
          newM1 = paste("M1",lm,lp,om,op,sep="_")
          newM2 = paste("M2",lm,lp,om,op,sep="_")
          names(out_df)[names(out_df) == 'M1'] <- newM1
          names(out_df)[names(out_df) == 'M2'] <- newM2
        }
      }
    }
  }
}

options(scipen = 999)
evenindex <- seq(2, ncol(out_df)-2, by = 2)
oddindex <- seq(3, ncol(out_df)-2, by = 2)

diff_M1 <- data.frame(matrix(NA,
                          nrow = (ncol(out_df)-1)/2 - 1, 
                          ncol = (ncol(out_df)-1)/2 - 1))
diff_M2 <- data.frame(matrix(NA,
                             nrow = (ncol(out_df)-1)/2 - 1, 
                             ncol = (ncol(out_df)-1)/2 - 1))

for (i in evenindex) {
  second_index <- seq(i+2, ncol(out_df), by = 2)
  for (j in second_index) {
    diff_M1[(j-2)/2,i/2] = mean(out_df[,i] - out_df[,j])
    colnames(diff_M1)[i/2] = colnames(out_df)[i]
    rownames(diff_M1)[(j-2)/2] = colnames(out_df)[j]
  }
}
View(diff_M1)

which(diff_M1 < 20 & diff_M1 > -20, arr.ind = TRUE)

colnames(diff_M1)[21]

for (i in oddindex) {
  second_index <- seq(i+2, ncol(out_df), by = 2)
  for (j in second_index) {
    diff_M2[(j-3)/2,(i-1)/2] = mean(out_df[,i] - out_df[,j])
    colnames(diff_M2)[(i-1)/2] = colnames(out_df)[i]
    rownames(diff_M2)[(j-3)/2] = colnames(out_df)[j]
  }
}

View(diff_M2)

which(diff_M2 < 20 & diff_M2 > -20, arr.ind = TRUE)
colnames(diff_M1)[2]

parameters_cov1 <- c(lambda_min = 1.5, lambda_plus = 1.0, omega_min = 0.01, omega_plus = 0.05, alpha_min = 1.5, beta_min = 0, alpha_plus = 1.0, beta_plus = 0)
parameters_cov2 <- c(lambda_min = 1.501, lambda_plus = 0.835, omega_min = 0.025, omega_plus = 0.067, alpha_min = 1.501, beta_min = 0, alpha_plus = 0.835, beta_plus = 0)
parameters_cov3 <- c(lambda_min = 1.495, lambda_plus = 0.504, omega_min = 0.077, omega_plus = 0.100, alpha_min = 1.495, beta_min = 0, alpha_plus = 0.504, beta_plus = 0)


out1 <- ode(y = state, times = seq(0, 5, by = 0.01), func = covariances, parms = parameters_cov1)
out2 <- ode(y = state, times = seq(0, 5, by = 0.01), func = covariances, parms = parameters_cov2)
out3 <- ode(y = state, times = seq(0, 5, by = 0.01), func = covariances, parms = parameters_cov3)

plot(out)

out1_df <- data.frame(t = out1[,1], M1 = out1[,2], M2 = out1[,3], V1 = out1[,4], V2 = out1[,5], C = out1[,6])
out1_df <- out1_df %>% mutate(D1 = sqrt(V1), D2 = sqrt(V2))



out2_df <- data.frame(t = out2[,1], M1 = out2[,2], M2 = out2[,3], V1 = out2[,4], V2 = out2[,5], C = out2[,6])
out2_df <- out2_df %>% mutate(D1 = sqrt(V1), D2 = sqrt(V2))

out3_df <- data.frame(t = out3[,1], M1 = out3[,2], M2 = out3[,3], V1 = out3[,4], V2 = out3[,5], C = out3[,6])
out3_df <- out3_df %>% mutate(D1 = sqrt(V1), D2 = sqrt(V2))

ggplot() + 
  geom_line(data = out1_df, aes(x=t,y=M1), color = "red") +
  geom_line(data = out2_df, aes(x=t,y=M1),color = "blue") +
  geom_line(data = out3_df, aes(x=t,y=M1),color = "green") +
  geom_point(data = simulation_py, aes(x = time, y = z_minus)) +
  xlim(4,5)

simulation_py <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/switching_results_avg1.csv") %>% tibble::as_tibble()

ggplot() + 
  geom_line(data = out1_df, aes(x=t,y=M2), color = "red") +
  geom_line(data = out2_df, aes(x=t,y=M2),color = "blue") +
  geom_line(data = out3_df, aes(x=t,y=M2),color = "green") +
  geom_point(data = simulation_py, aes(x = time, y = z_plus)) +
  xlim(4,5)

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
  c2 = (2*omega_plus*z0[1] - (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  zmin = exp((lambda_minus + lambda_plus)*t/2.)*(c1*(lambda_minus - lambda_plus + sqrt(delta))*exp(sqrt(delta)*t/2.) + c2*2*omega_minus*exp(-sqrt(delta)*t/2.))
  return(zmin)
}

zplus_ode <- function(t, z0, lambda_minus, lambda_plus, omega_minus, omega_plus){
  delta = (lambda_minus - lambda_plus)^2 + 4*omega_minus*omega_plus
  c1 = ((lambda_minus - lambda_plus + sqrt(delta))*z0[1] + 2*omega_minus*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  c2 = (2*omega_plus*z0[1] - (lambda_minus - lambda_plus + sqrt(delta))*z0[2])/((lambda_minus - lambda_plus + sqrt(delta))^2 + 4*omega_minus*omega_plus)
  zplus = exp((lambda_minus + lambda_plus)*t/2.)*(c1*2*omega_plus*exp(sqrt(delta)*t/2.) + c2*(lambda_plus - lambda_minus - sqrt(delta))*exp(-sqrt(delta)*t/2.))
  return(zplus)
}



t <- seq(0,1., by = 0.01)
z0 <- as.array(c(1000,100))

zmin_ode(0, z0, 1.5, 1.0, .01, .01)

y = lapply(t, zplus_ode, z0=z0, lambda_minus=1.5, lambda_plus=1.0, omega_minus=0.01, omega_plus=0.01) %>% unlist()

plot(t, y)

py_simulation <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.0_005_001/switching_results_avg.csv") %>%
  tibble::as_tibble()
colnames(py_simulation) <- c("step","t","Z-","Z+")


ggplot() +
  geom_point(data = out_df, aes(x = t, y = M1), size = 0.5) +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = 0.01)) +
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 0.938, omega_minus = 0.02, omega_plus = 0.013), color = "red")

ggplot() +
  geom_point(data = out_df, aes(x = t, y = M2), size = 0.5) +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 1.0, omega_minus = 0.01, omega_plus = 0.01)) +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.5, lambda_plus = 0.938, omega_minus = 0.02, omega_plus = 0.013), color = "red")


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
  stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 1.0, lambda_plus = 1.5, omega_minus = 0.001, omega_plus = 0.01), color = 'blue') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 1.0, lambda_plus = 1.5, omega_minus = 0.001, omega_plus = 0.01), color = 'blue') +
  geom_point(data = py_simulation, aes(x = t, y = `Z-`)) +
  #geom_point(data = py_simulation, aes(x = t, y = `Z+`)) +
  #stat_function(fun = zmin_ode, args = list(z0, lambda_minus = 0.995, lambda_plus = 1.474, omega_minus = 0.016, omega_plus = 0.015), color = 'red') +
  stat_function(fun = zplus_ode, args = list(z0, lambda_minus = 0.995, lambda_plus = 1.474, omega_minus = 0.016, omega_plus = 0.015), color = 'red') +
  xlim(1,4)

