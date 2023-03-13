library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)

# define parameters
alpha_min = 20;beta_min = 0;alpha_plus = 20;beta_plus = 0;omega_p = 0.1;omega_m = 0.01

# population starting with 1 cell in state -
Z_minus = 1; Z_plus = 0; t = 0
Z <- c(Z_minus,Z_plus)

# define stoichiometric vectors
o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
o <- rbind(o1,o2,o3,o4,o5,o6)

# create tibble to store Z values for each t
output <- tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2])

while (t < 0.65) { 
  a1 = alpha_min*Z[1]
  a2 = beta_min*Z[1]
  a3 = omega_p*Z[1]
  a4 = alpha_plus*Z[2]
  a5 = beta_plus*Z[2]
  a6 = omega_m*Z[2]
  a <- c(a1,a2,a3,a4,a5,a6)
  a0 <- sum(a)
  anorm <- a/a0
  tau <- rexp(n = 1, rate = a0)
  t <- t + tau
  i <- rcat(1,anorm)
  Z <- Z+o[i,]
  print(t)
  output <- bind_rows(output,tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]))
}