library(LaplacesDemon)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(reshape2)
library(patchwork)
library(plotly)

simulation <- function(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, index) {
  # initialize Z and t
  Z_minus = 1; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)
  # tibble to store evolution of populations with time
  # Zt_output <- tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]) 
  while (t < 1.2) {
    a1 = alpha_min*Z[1]
    a2 = beta_min*Z[1]
    a3 = omega_p*Z[1]
    a4 = alpha_plus*Z[2]
    a5 = beta_plus*Z[2]
    a6 = omega_m*Z[2]
    a <- c(a1,a2,a3,a4,a5,a6)
    a0 <- sum(a)
    anorm <- a/a0 # vector with normalized probabilities
    tau <- rexp(n = 1, rate = a0) # sample time step
    t <- t + tau
    i <- rcat(1,anorm) # sample event
    Z <- Z+o[i,]
    # Zt_output <- bind_rows(Zt_output,tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]))
  }
  # return final populations
  output = tibble("Z-" = Z[1],"Z+" = Z[2], "time" = t)
  return(output)
}

Zt_log_plots <- function(dat) {
  
  # plot the evolution of both populations on same graph (log scale)
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point(size=0.5) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
    scale_y_continuous(trans = 'log10') +
    theme(legend.position = "none")
  
  return(plot1)
}
Zt_norm_plots <- function(dat){
  #add columns with sum and proportions
  dat <- mutate(dat, sum = `Z+`+`Z-`, "Z- ratio" = `Z-`/sum, "Z+ ratio" = `Z+`/sum)
  # plot relative proportions
  plot2 <- dat %>% 
    dplyr::select(t, dplyr::contains("ratio")) %>% 
    melt(id="t", variable.name="type", value.name="ratio") %>% 
    ggplot() +
    geom_line(aes(x=t, y=ratio, color=type),linewidth=0.8) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) +
    theme(legend.position = "none") +
    ylab(bquote(Z/Z[tot]))
  
  return(plot2)
}
Pz_exp <- function(dat,alpha_min,alpha_plus,omega_m,omega_p) {
  time <- mean(dat$time)
  if (alpha_min > alpha_plus) {
    exprate_minus <- exp(-(alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time)
    exprate_plus <- (alpha_min - alpha_plus)/omega_p*exp(-(alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*time)
  } else if (alpha_min == alpha_plus) {
    exprate_minus <- 1/(cosh(sqrt(omega_m*omega_p)*time))*exp(-alpha_min*time)
    exprate_plus <- 1/(sinh(sqrt(omega_m*omega_p)*time))*sqrt(omega_m/omega_p)*exp(-alpha_min*time)
  }
  rates <- list(exprate_minus,exprate_plus)
  return(rates)
}
Pz_exp_plot <- function(dat,rate) {
  exprate_minus = as.numeric(rate[1])
  exprate_plus = as.numeric(rate[2])
  my_palette <- c(rgb(102,204,102,maxColorValue = 255),"#D5D139")
  # no line
  plot1 <- dat %>% 
    reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
    ggplot() +
    geom_histogram(aes(x=Z, y=after_stat(count), fill=type), color="black", bins=180) +
    scale_fill_manual(values=my_palette) +
    facet_grid(type~., scales="free_y") #+ xlim(0,0.75e5) + ylim(0,170) 
  
  # with line
  plot2 <- dat %>% 
    reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
    dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate=exprate_minus), dexp(Z,rate=exprate_plus))) %>%
    ggplot() +
    geom_histogram(aes(x=Z, y=after_stat(density), fill=type), color="black", bins=150) +
    scale_fill_manual(values=my_palette) +
    geom_line(aes(x=Z,y=rate)) +
    facet_grid(type~., scales="free_y") #+ ylim(0,5e-5) + xlim(0,5e5)
  
  plots <- list(plot1,plot2)
  return(plots)
}

my_palette <- c(rgb(102,204,102,maxColorValue = 255),"#D5D139")
alpha_min = 10.1;beta_min = 0;alpha_plus = 10.0;beta_plus = 0;omega_p = 0.1;omega_m = 0.5
lambda_min = alpha_min - beta_min
lambda_plus = alpha_plus - beta_plus

options(scipen = 0)

# run simulation 1000 times
output_list = lapply(1:1000, function(i){
  x <- simulation(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, i)
  x$step <- i
  print(i)
  x
})

final = output_list %>% do.call(rbind, .)

# save file in rds format
saveRDS(final, file = paste0("./simulations/final",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,".rds"))

# TIME EVOLUTION

# open rds file
final_time1 <- readRDS("./simulations_time/output_[20_20_01][0.673].rds")
final_time2 <- readRDS("./simulations_time/output_[15_10_005][0.989].rds")
final_time3 <- readRDS("./simulations_time/output_[15_10_001][0.983].rds")

Zt_norm_plots(final_time1) 

patch1 <- Zt_log_plots(final_time1) + Zt_log_plots(final_time2) + Zt_log_plots(final_time3)
patch1

patch2 <- Zt_norm_plots(final_time1) + Zt_norm_plots(final_time2) + Zt_norm_plots(final_time3)
patch2

patch3 <- patch1 / patch2
patch3
ggsave(path = "./imgs", width = 7, height = 3.5, device='tiff', dpi=500, filename = "Zt_15_10")

# MARGINAL DISTRIBUTIONS: HISTOGRAMS

final1 <- readRDS("./simulations/10_20_01_001 [0.8]/final.rds")
final2 <- readRDS("./simulations/final[10_20_005]1000[0.8].rds")
final3 <- readRDS("./simulations/10_20_001_01 [0.8]/final.rds")

# horizontal
plot1 <- final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 5.5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=120, color = 'black', na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") 
plot1

plot2 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 5.5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") 
plot2

plot3 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 5.5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=80, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(20000,40000))
plot3

patch1 <- plot1 + plot2 + plot3
patch1

# vertical
plot1 <-final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  #filter((type == "Z-" & Z < 8.25e5) | type == "Z+") %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),color = "black",bins=80, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(~type) +
  theme(legend.position = "none")
plot1

plot2 <-final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(Z < 5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),color = "black",bins=80, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(~type) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(20000,40000))
plot2

plot3 <-final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(Z < 5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),color = "black",bins=80, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(~type) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(20000,40000))
plot3

patch4 <- plot1 + plot2 + plot3
patch4

# FITS

# for exponential distributions
rate1 <- Pz_exp(final,alpha_min,alpha_plus,omega_m,omega_p)
rate2 <- Pz_exp(final2,alpha_min,alpha_plus,0.05,0.05)
rate3 <- Pz_exp(final3,alpha_min,alpha_plus,0.1,0.01)

final_pz_plots <- Pz_exp_plot(final,rate)
final_pz_plots[1]
final_pz_plots[2]

final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  #filter((type == "Z-" & Z < 3e5) | type == "Z+") %>% 
  dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate = as.numeric(rate1[1])), dexp(Z,rate = as.numeric(rate1[2])))) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(density), fill=type), bins=180) +
  scale_fill_manual(values=my_palette) +
  geom_line(aes(x=Z,y=rate)) +
  facet_grid(~type) 

# individual plots
# Z-
plot1a <- final %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), color = "black", fill = rgb(102,204,102,maxColorValue = 255), bins = 100, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[1])), linewidth=0.8)
plot1a

plot2a <- final2 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), color = "black", fill = rgb(102,204,102,maxColorValue = 255), bins = 100, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[1])), size=0.8) +
  xlim(0,60000) #+ ylim(0,0.00015)
plot2a

plot3a <- final3 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), color = "black", fill = rgb(102,204,102,maxColorValue = 255), bins = 100, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[1])), size=0.8) +
  xlim(0,60000) +
  scale_y_continuous(breaks = c(0.000025,0.000050,0.000075,0.000100,0.000125))
plot3a

patch5 <- plot1a + plot2a + plot3a
patch5

# Z+

plot1b <- ggplot(final,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[2])), size=0.8) 
  #scale_x_continuous(breaks=c(30000,60000,90000),limits = c(0,90000))
plot1b

plot2b <- ggplot(final2,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[2])), size=0.8) +
  xlim(0,1500) + ylim(0,0.025)
plot2b

plot3b <- ggplot(final3,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[2])), size=0.8) +
  xlim(0,1500) + ylim(0,0.025)
plot3b

patch6 <- plot1b + plot2b + plot3b
patch6

patch1 / patch2 / patch4 / patch5 / patch6

# when not exponential

# Z+ follows a POWER-LAW

time1 <- mean(final1$time)
time2 <- mean(final2$time)
time3 <- mean(final3$time)

powerlaw_coeff <- function(alpha_min,alpha_plus,lambda_min,lambda_plus,omega_m,omega_p,time) {
  return(alpha_min/(lambda_plus*alpha_plus)*omega_p*pi/sin(pi*lambda_min/lambda_plus)*(alpha_plus/lambda_plus)^(lambda_min/lambda_plus)*exp((lambda_min+omega_m*omega_p/(lambda_plus-lambda_min)*lambda_min/lambda_plus)*time1)/gamma(1-lambda_min/lambda_plus))
}

powerlaw_coeff1 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus, omega_m = 0.01,omega_p = 0.1,time = time1)
powerlaw_coeff2 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus,omega_m = 0.05,omega_p = 0.05,time = time2)
powerlaw_coeff3 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus,omega_m = 0.1,omega_p = 0.01,time = time3)

powerlaw <- function(x,coeff,exponent) {
  return(coeff*x^exponent)
}

ggplot(final1,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 150, na.rm = TRUE) +
  stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff1, exponent = -(1+lambda_min/lambda_plus)), size=0.8) +
  xlim(0,5e6) + ylim(0,0.4e-5)

ggplot(final2,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 150, na.rm = TRUE) +
  stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff2, exponent = -(1+lambda_min/lambda_plus)), size=0.8) +
  xlim(0,8e5) + ylim(0,2.5e-5)

ggplot(final3,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff3, exponent = -(1+lambda_min/lambda_plus)), size=0.8) +
  xlim(0,8e4) + ylim(0,8e-5)

# Z-



analytic1 <- read.csv("./imgs/[10_20_01]/P(Z-)[10_20_01_001].csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/[10_20_01]/P(Z+)[10_20_01_001][0.8] mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))

analytic2 <- read.csv("./imgs/[10_20_005]/P(Z-)[10_20_005][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/[10_20_005]/P(Z+)[10_20_005][0.8]2 mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))

analytic3 <- read.csv("./imgs/[10_20_001]/P(Z-)[10_20_001_01][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/[10_20_001]/P(Z+)[10_20_001_01][0.8] mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))

final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter((type=="Z+" & Z < 50000) | (type=="Z-")) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,fill=type, y=after_stat(density)),
                 color = "black",bins=80, position = "identity") +
  geom_line(data=analytic %>% 
              dplyr::group_by(type) %>% 
              dplyr::mutate(nn=dplyr::n()), aes(x=Z,y=P)) +
  facet_grid(~type) +
  scale_fill_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139"))

# individual fit
# Z-
plot4 <- final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255), na.rm = TRUE) +
  geom_line(data=analytic1 %>% filter(type=="Z-"), aes(x=Z,y=P), linewidth=0.5) +
  xlim(0,45000)
plot4

plot5 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-" & Z<40000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255)) +
  geom_line(data=analytic2 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) #+
#xlim(0,40000)
plot5

plot6 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255), na.rm = TRUE) +
  geom_line(data=analytic3 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) +
  xlim(0,45000)
plot6

patch2 <- plot4 + plot5 + plot6
patch2

# Z+
plot7 <- final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z+" & Z<60000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=100,fill="#D5D139", na.rm = TRUE) +
  geom_line(data=analytic1 %>% filter(type=="Z+" & Z<60000), aes(x=Z,y=P)) +
  xlim(0,60000) +
  ylim(0,2e-4)
plot7

analytic2 %>% filter(type == "Z+" & Z<80000) %>% ggplot(aes(x=Z,y=P)) + geom_line()


plot8 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z+" & Z < 60000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=100, fill="#D5D139", na.rm = TRUE) +
  geom_line(data=analytic2 %>% filter(type=="Z+" & Z < 60000), aes(x=Z,y=P), na.rm = TRUE) +
  xlim(0,60000) +
  ylim(0,2e-4)
plot8

plot9 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z+" & Z < 60000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=120, fill="#D5D139", na.rm = TRUE) +
  geom_line(data=analytic3 %>% filter(type=="Z+" & Z < 60000), aes(x=Z,y=P), na.rm = TRUE) +
  xlim(0,60000) +
  ylim(0,2e-4)
plot9

patch3 <- plot7 + plot8 + plot9
patch3

patch1 / patch2 / patch3

# plot histograms on same graph
p <- ggplot() +
  geom_histogram(aes(x = final$`Z-`, fill = "Z+"), alpha = 0.5) +
  geom_histogram(aes(x = final$`Z+`, fill = "Z-"), alpha = 0.5) +
  scale_fill_manual(values = c("Z-" = "red", "Z+" = "green")) +
  xlab("Z") +
  ylab("count") 
p

final <- final %>% 
  filter(`Z-` > 30)

jp <- matrix(0:0, nrow = 10000, ncol = 300)

for (i in 0:10000) {
  for (j in 0:300) {
    if (j == as.integer(0.02*i)) {
      jp[i,j] = rexp(1, rate = as.numeric(rate1[1]))
    } else {
      jp[i,j] = 0
    }
  }
}

x <- 1:10000
y <- 1:300
contour(x,y,jp)

# plot Z- VS Z+
final2 %>% 
  filter(`Z-` < 10000) %>% 
  mutate(sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum) %>% 
  ggplot(aes(x=`Z-`,y=`Z+`)) + 
  geom_point(color = 'cyan4') +
  geom_line(aes(x=`Z-`, y=`Z-`*omega_p/(alpha_min-alpha_plus)), color='black')
  #ylim(0,max(max(final$`Z-`),max(final$`Z+`)))
  
final2 %>% 
  mutate(sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum) %>% 
  ggplot(aes(x=`Z- rate`,y=`Z+ rate`)) + 
  geom_point(color = 'cyan4') +
  ylim(0,1)

rm(list = ls())

