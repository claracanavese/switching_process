library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)
library(ggmuller)
library(patchwork)
library(tibble)
library(devtools)
library(easypar)
library(RColorBrewer)

# define parameters
alpha_min = 15;beta_min = 0;alpha_plus = 10;beta_plus = 0;omega_p = 0.1;omega_m = 0.01

# population starting with 1 cell in state -
Z_minus = 1; Z_plus = 0; time = 0
Z <- c(Z_minus,Z_plus)

# define stoichiometric vectors
o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(0,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,0)
o <- rbind(o1,o2,o3,o4,o5,o6)

# create tibble to store Z values for each t
output <- tibble("t" = time,"Z-" = Z[1],"Z+" = Z[2])

while (time < 1.0) { 
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
  time <- time + tau
  print(time)
  i <- rcat(1,anorm)
  Z <- Z+o[i,]
  output <- bind_rows(output,tibble("t" = time,"Z-" = Z[1],"Z+" = Z[2]))
}

times_simulation <- function(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, index) {
  # population starting with 1 cell in state -
  Z_minus = 1; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)
  
  # create tibbles to store switching times
  min_plus <- tibble(t1 = numeric())
  plus_min <- tibble(t2 = numeric())
  
  while (t < 1.0) { 
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
    if (i == 3) {
      min_plus <- bind_rows(min_plus,tibble("t1"=t))
    } else if (i == 6) {
      plus_min <- bind_rows(plus_min,tibble("t2"=t))
    }
    Z <- Z+o[i,]
  }
  switching_times <- bind_rows(min_plus[1,],plus_min[1,])
  return(switching_times)
}

# failed attempt to use easypar :(
times_simulation2 <- function(alpha_min) {
  # population starting with 1 cell in state -
  Z_minus = 1; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)
  
  # create tibbles to store switching times
  min_plus <- tibble(t1 = numeric())
  plus_min <- tibble(t2 = numeric())
  
  while (t < 0.5) { 
    a1 = alpha_min*Z[1]
    a2 = 0*Z[1]
    a3 = 0.01*Z[1]
    a4 = 10*Z[2]
    a5 = 0*Z[2]
    a6 = 0.1*Z[2]
    a <- c(a1,a2,a3,a4,a5,a6)
    a0 <- sum(a)
    anorm <- a/a0
    tau <- rexp(n = 1, rate = a0)
    t <- t + tau
    i <- rcat(1,anorm)
    if (i == 3) {
      min_plus <- bind_rows(min_plus,tibble("t1"=t))
    } else if (i == 6) {
      plus_min <- bind_rows(plus_min,tibble("t2"=t))
    }
    Z <- Z+o[i,]
  }
  switching_times <- bind_rows(min_plus[1,],plus_min[1,])
  return(min_plus[1,])
}

el <- list(alpha_min)
inputs <- list(el)[rep(1,10)]

easypar::run(FUN = times_simulation2, PARAMS = inputs, parallel = TRUE)

# standard procedure
times = lapply(1:100, function(i){
  x <- times_simulation(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, i)
  x$step <- i
  print(i)
  x
})
times = times %>% do.call(rbind, .)
sapply(times, mean, na.rm = TRUE)
sapply(times, sd, na.rm = TRUE)
mean(times$t1, na.rm= TRUE)
mean(times$t2, na.rm = TRUE)

saveRDS(times, file = paste0("./simulations_time/times",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,".rds"))

# save data in rds file
saveRDS(output, file = paste0("./simulations_time/output_[15_10_01][0.96].rds"))
output <- melt(output[,1:3], value.name = "Z", id = "t")
plot <- output %>% ggplot(aes(time,Z, col=variable)) + ylab("Z") + geom_point() +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139"))
plot

# ODE comparison

for (i in 0:14) {
  final_time <- paste("final_time", i, sep = "_")
  assign(final_time,read.csv(paste("./GitHub/switching_process/Gillespy2/1.2_1.5_0.015_0.005_8t_81p/switching_results_",i,".csv", sep="")) %>%
    tibble::as_tibble())
}

for (i in 0:14) {
  final_time <- paste("final_time1", i, sep = "_")
  assign(final_time,read.csv(paste("./GitHub/switching_process/Gillespy2/1.2_1.5_0.005_0.015_8t_81p/switching_results_",i,".csv", sep="")) %>%
           tibble::as_tibble())
}

final_time_15 <- read.csv("./GitHub/switching_process/Gillespy2/1.5_1.2_0.015_0.005_8t_81p/switching_results_avg.csv") %>%
  tibble::as_tibble()

ode_sol <- readRDS("./simulations_time/ode_[10_15_00]_1.2.rds")
ode_sol <- data.frame(t = ode_sol[,1], X = ode_sol[,2], Y = ode_sol[,3])
options(scipen = 0)

# plotting together
plotmin <- ggplot() +
  geom_line(data = final_time_0[,-1], aes(time,z_minus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_1[,-1], aes(time,z_minus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time_2[,-1], aes(time,z_pminus, color = "2"), linewidth=0.8) +
  geom_line(data = final_time_3[,-1], aes(time, z_minus, color = "3"), linewidth=0.8) +
  geom_line(data = final_time_4[,-1], aes(time,z_minus, color = "4"), linewidth=0.8) +
  geom_line(data = final_time_5[,-1], aes(time,z_minus, color = "5"), linewidth=0.8) +
  geom_line(data = final_time_6[,-1], aes(time,z_minus, color = "6"), linewidth=0.8) +
  geom_line(data = final_time_7[,-1], aes(time,z_minus, color = "7"), linewidth=0.8) +
  geom_line(data = final_time_8[,-1], aes(time,z_minus, color = "8"), linewidth=0.8) +
  geom_line(data = final_time_9[,-1], aes(time,z_minus, color = "9"), linewidth=0.8) +
  geom_line(data = final_time_10[,-1], aes(time,z_minus, color = "10"), linewidth=0.8) +
  geom_line(data = final_time_11[,-1], aes(time,z_minus, color = "11"), linewidth=0.8) +
  geom_line(data = final_time_12[,-1], aes(time,z_minus, color = "12"), linewidth=0.8) +
  geom_line(data = final_time_13[,-1], aes(time,z_minus, color = "13"), linewidth=0.8) +
  geom_line(data = final_time_14[,-1], aes(time,z_minus, color = "14"), linewidth=0.8) +
  geom_line(data = final_time_15[,-1], aes(time,z_minus, color = "15"), linewidth=0.8) #+
  # geom_line(data = ode_sol, aes(x = t, y = X), color = "black", linewidth=1.5) +
  #scale_color_brewer(palette = "Paired")
plotmin

plotmin <- ggplot() +
  geom_line(data = final_time_0[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_1[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_2[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_3[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_4[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_5[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_6[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_7[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_8[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_9[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_10[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time1_0[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_1[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_2[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_3[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_4[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_5[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_6[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_7[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_8[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_9[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time1_10[,-1], aes(time,z_plus, color = "1"), linewidth=0.8)

plotmin

plotplus <- ggplot() +
  geom_line(data = final_time_0[,-1], aes(time,z_plus, color = "0"), linewidth=0.8) +
  geom_line(data = final_time_1[,-1], aes(time,z_plus, color = "1"), linewidth=0.8) +
  geom_line(data = final_time_2[,-1], aes(time,z_plus, color = "2"), linewidth=0.8) +
  geom_line(data = final_time_3[,-1], aes(time,z_plus, color = "3"), linewidth=0.8) +
  geom_line(data = final_time_4[,-1], aes(time,z_plus, color = "4"), linewidth=0.8) +
  geom_line(data = final_time_5[,-1], aes(time,z_plus, color = "5"), linewidth=0.8) +
  geom_line(data = final_time_6[,-1], aes(time,z_plus, color = "6"), linewidth=0.8) +
  geom_line(data = final_time_7[,-1], aes(time,z_plus, color = "7"), linewidth=0.8) +
  geom_line(data = final_time_8[,-1], aes(time,z_plus, color = "8"), linewidth=0.8) +
  geom_line(data = final_time_9[,-1], aes(time,z_plus, color = "9"), linewidth=0.8) +
  geom_line(data = final_time_10[,-1], aes(time,z_plus, color = "10"), linewidth=0.8) +
  geom_line(data = final_time_11[,-1], aes(time,z_plus, color = "11"), linewidth=0.8) +
  geom_line(data = final_time_12[,-1], aes(time,z_plus, color = "12"), linewidth=0.8) +
  geom_line(data = final_time_13[,-1], aes(time,z_plus, color = "13"), linewidth=0.8) +
  geom_line(data = final_time_14[,-1], aes(time,z_plus, color = "14"), linewidth=0.8) +
  geom_line(data = final_time_15[,-1], aes(time,z_plus, color = "15"), linewidth=0.8)
  #geom_line(data = ode_sol, aes(x = time, y = Y), color = "black", linewidth=1) +
  #scale_color_brewer(palette = "Paired")
plotplus

plotmin / plotplus


plot1 <- final_time[,-1] %>% 
  reshape2::melt(id=c("t"), variable.name="type", value.name="Z") %>%
  ggplot() +
  geom_point(aes(x=time, y=Z, color = type)) +
  scale_fill_manual(values=my_palette) +
  xlim(0,xmax) +
  facet_grid(~type)
plot1
plot2 <- ode_sol %>% 
  reshape2::melt(id=c("t"), variable.name="type", value.name="Z") %>%
  ggplot() +
  geom_point(aes(x=time, y=Z, color = type)) +
  scale_fill_manual(values=my_palette) +
  xlim(0,xmax) +
  facet_grid(~type)  

plot1 / plot2



Zt_plots <- function(dat) {
  #add columns with sum and proportions
  dat <- mutate(dat, sum = z_plus+z_plus, "Z- rate" = z_plus/sum, "Z+ rate" = z_plus/sum)
  
  # plot the evolution of both populations on same graph 
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(time,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) #+ 
    #scale_y_continuous(trans = 'log10')

  plot2 <- dat1 %>% ggplot(aes(time,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10')
  
  # plot relative proportions
  plot3 <- dat %>% 
    dplyr::select(time, dplyr::contains("rate")) %>% 
    melt(id="t", variable.name="type", value.name="rate") %>% 
    ggplot() +
    geom_line(aes(x=time, y=rate, color=type)) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139"))
  output <- list(plot1,plot2,plot3)
  return(output)
}
final <- Zt_plots(final_time)
final[1]
final[2]
final[3]


# plots to verify asymptotic limit

# NO
final_time <- readRDS("./simulations_time/output_[10_20_005][1.05].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot3 <- final_time1 %>% 
  mutate(eexp1=omega_p/alpha_plus*pi/sin(pi*alpha_min/alpha_plus)*omega_m/(alpha_plus-alpha_min)*exp((alpha_plus+omega_m*omega_p/(alpha_plus-alpha_min))*time)+exp(alpha_min*time)) %>% 
  mutate(eexp2=omega_p/alpha_plus*pi/sin(pi*alpha_min/alpha_plus)*exp((alpha_plus+omega_m*omega_p/(alpha_plus-alpha_min))*(time-0.1))) %>% 
  ggplot(aes(time, col=variable)) + 
  ylab("Z") + 
  geom_point(aes(y=Z),size=0.5) +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10') +
  geom_line(aes(y=eexp1), color='black') +
  geom_line(aes(y=eexp2), color='black')
plot3  

# YES
alpha_min = 20;alpha_plus = 20;omega_p = 0.1;omega_m = 0.01
final_time <- readRDS("./simulations_time/output_[20_20_01][0.673].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot2 <- final_time1 %>% 
  filter(Z<1e6) %>% 
  mutate(eexp1=cosh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time)) %>% 
  mutate(eexp2=sinh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time)*sqrt(omega_p/omega_m)) %>% 
  ggplot(aes(time, col=variable)) + 
  ylab("Z") + 
  geom_point(aes(y=Z),size=1.5) +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(color = NULL) +
  theme(legend.text = element_text(size = 14)) +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  geom_line(aes(y=eexp1), color='black', linewidth = 1) +
  geom_line(aes(y=eexp2), color='black', linewidth = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(family = "Arial")) +
  ggtitle(bquote(~ lambda['-']==20 ~ lambda['+']==20)) +
  ylim(0,1e5) #+
  #scale_y_continuous(trans = 'log10')
plot2

alpha_min = 15;alpha_plus = 10;omega_p = 0.15;omega_m = 0.01
final_time <- readRDS("./simulations_time/output_[15_10_001][0.983].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot1 <- final_time1 %>% 
  filter(Z<1e6) %>% 
  mutate(eexp1=exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time)) %>% 
  mutate(eexp2=omega_p/(alpha_min-alpha_plus)*exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time)) %>% 
  ggplot(aes(time, col=variable)) + 
  ylab("Z") + 
  #geom_point(aes(y=Z),size=1.5) +
  #scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  geom_line(aes(y=eexp1,color="Z-"), color=rgb(102,204,102,maxColorValue = 255), linewidth = 1) +
  geom_line(aes(y=eexp2,color="Z+"), color='#D5D139', linewidth = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #labs(color = NULL) +
  theme(legend.text = element_text(size = 14)) +
  #ggtitle(bquote(~ lambda['-']==15 ~ lambda['+']==10)) +
  ylim(0,1e5) #+
  #theme(plot.title = element_text(family = "Arial")) #+
  #scale_y_continuous(trans = 'log10')
plot1

ggplot() +
  geom_function(fun = function(time) exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time),color = rgb(102,204,102,maxColorValue = 255), linewidth = 1.2) +
  geom_function(fun = function(time) omega_p/(alpha_min-alpha_plus)*exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time),color='#D5D139', linewidth = 1.2) +
  geom_function(fun = function(time) exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time)+omega_p/(alpha_min-alpha_plus)*exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time),color='black', linewidth = 1.2) +
  xlim(0,0.5)# +
  #scale_y_continuous(trans = 'log10')

ggplot() +
  geom_function(fun = function(time) cosh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time),color = rgb(102,204,102,maxColorValue = 255), linewidth = 1.2) +
  geom_function(fun = function(time) sinh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time)*sqrt(omega_p/omega_m),color='#D5D139', linewidth = 1.2) +
  geom_function(fun = function(time) cosh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time)+sinh(sqrt(omega_m*omega_p)*time)*exp(alpha_min*time)*sqrt(omega_p/omega_m),color='black', linewidth = 1.2) +
  xlim(0,0.5)

plot1 + plot2
ggsave("./imgs/asymptotic.png",dpi=600)

# MULLER PLOT
row_odd <- seq_len(nrow(output)) %% 2; output <- output[row_odd == 1, ]
row_odd <- seq_len(nrow(output)) %% 2; output <- output[row_odd == 1, ]
row_odd <- seq_len(nrow(output)) %% 2; output <- output[row_odd == 1, ]
row_odd <- seq_len(nrow(output)) %% 2; output <- output[row_odd == 1, ]
row_odd <- seq_len(nrow(output)) %% 2; output <- output[row_odd == 1, ]
# row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
output <- output %>% reshape2::melt(id=c("t"), variable.name="type", value.name="Z") 

final_time1 <- readRDS("./simulations_time/output_[10_20_01][1.0].rds")
row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
# row_odd <- seq_len(nrow(final_time1)) %% 2; final_time1 <- final_time1[row_odd == 1, ]
final_time1 <- final_time1 %>% reshape2::melt(id=c("t"), variable.name="type", value.name="Z") 

final_time2 <- readRDS("./simulations_time/output_[10_20_005][1.0].rds")
row_odd <- seq_len(nrow(final_time2)) %% 2; final_time2 <- final_time2[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time2)) %% 2; final_time2 <- final_time2[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time2)) %% 2; final_time2 <- final_time2[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time2)) %% 2; final_time2 <- final_time2[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time2)) %% 2; final_time2 <- final_time2[row_odd == 1, ]
final_time2 <- final_time2 %>% reshape2::melt(id=c("t"), variable.name="type", value.name="Z") 

final_time3 <- readRDS("./simulations_time/output_[10_20_001][1.0].rds")
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
row_odd <- seq_len(nrow(final_time3)) %% 2; final_time3 <- final_time3[row_odd == 1, ]
final_time3 <- final_time3 %>% reshape2::melt(id=c("t"), variable.name="type", value.name="Z") 

# MULLER
# create edges dataframe
edges <- data.frame(Parent = c("Z.0","Z-"), Identity = c("Z-","Z+"))

# rename columns and add row with ancestral population
colnames(output) <- c("Time","Identity","Population")
dat1 = rbind(output, data.frame("Time"=-1,"Identity"="Z.0","Population"=1))

colnames(final_time1) <- c("Time","Identity","Population")
dat1 = rbind(final_time1, data.frame("Time"=-1,"Identity"="Z.0","Population"=1))
colnames(final_time2) <- c("Time","Identity","Population")
dat2 = rbind(final_time2, data.frame("Time"=-1,"Identity"="Z.0","Population"=1))
colnames(final_time3) <- c("Time","Identity","Population")
dat3 = rbind(final_time3, data.frame("Time"=-1,"Identity"="Z.0","Population"=1))

# create and plot muller plots
Muller_df1 <- get_Muller_df(edges,dat1)
mp1 <- Muller_plot(Muller_df1, add_legend = TRUE, xlab = "Time") + 
  scale_fill_manual(values = c(rgb(102,204,102,maxColorValue = 255),"lightblue","#D5D139")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(bquote(~ omega['-']==1.1 ~ omega['+']==1.0))
mp1

Muller_df2 <- get_Muller_df(edges,dat2)
mp2 <- Muller_plot(Muller_df2, add_legend = TRUE, xlab = "Time") + 
  scale_fill_manual(values = c(rgb(102,204,102,maxColorValue = 255),"lightblue","#D5D139")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(bquote(~ omega['-']== ~ omega['+']==0.05))
mp2

Muller_df3 <- get_Muller_df(edges,dat3)
mp3 <- Muller_plot(Muller_df3, add_legend = TRUE, xlab = "Time") + 
  scale_fill_manual(values = c(rgb(102,204,102,maxColorValue = 255),"lightblue","#D5D139")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(bquote(~ omega['-']==0.1 ~ omega['+']==0.01))
mp3

mp1 + mp2 + mp3


mp2 <-ggplot(Muller_df, aes_string(x = "Time", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
  geom_area() +
  theme(legend.position = "right") +
  guides(linetype = FALSE, color = FALSE) + 
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage") +
  scale_fill_manual(name = "Identity", values = c(rgb(102,204,102,maxColorValue = 255),"lightblue","#D5D139")) +
  scale_color_manual(values = c(rgb(102,204,102,maxColorValue = 255),"lightblue","#D5D139"))

ggplot(Muller_df, aes_string(x = "Time", y = "Frequency", fill = "Identity", colour = "Identity")) +
  geom_area() + theme(legend.position = "right") + 
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = 25 * (0:4), name = "Percentage")

rm(list=ls())
