library(LaplacesDemon)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(reshape2)
library(patchwork)
library(plotly)

simulation <- function(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, index) {
  # initialize Z and t
  Z_minus = 100; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)

  while (t < 0.5) {
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
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
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
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
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
alpha_min = 20;beta_min = 0;alpha_plus = 20;beta_plus = 0;omega_p = 0.1;omega_m = 0.01
lambda_min = alpha_min - beta_min
lambda_plus = alpha_plus - beta_plus

options(scipen = 0)

# run simulation 1000 times
output_list = lapply(1:500, function(i){
  x <- simulation(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, i)
  x$step <- i
  print(i)
  x
})

final = output_list %>% do.call(rbind, .)

# save file in rds format
saveRDS(final, file = paste0("./simulations/final",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,"start100.rds"))

t <- mean(final$time)
meanmin = alpha_min*cosh(sqrt(omega_m*omega_p)*t)/lambda_min*exp(lambda_min*t)
ratemin = lambda_min/(alpha_min*cosh(sqrt(omega_m*omega_p)*t))*exp(-lambda_min*t)

final %>% 
  ggplot(aes(x=`Z-`)) +
  geom_histogram(aes(y=after_stat(density)),bins=40,na.rm = TRUE, fill = "pink") +
  geom_function(fun=dnorm, args = list(mean = meanmin*100, sd = meanmin*10))

meanpl = sqrt(omega_p/omega_m)*alpha_min*sinh(sqrt(omega_m*omega_p)*t)/lambda_min*exp(lambda_min*t)

final %>% 
  ggplot(aes(x=`Z+`)) +
  geom_histogram(aes(y=after_stat(density)),bins=40,na.rm = TRUE, fill = "pink") +
  geom_function(fun=dnorm, args = list(mean = meanpl*100, sd = meanpl*10))

ggplot() +
  geom_function(fun=dnorm, args = list(mean = meanth, sd = meanth**2)) +
  xlim(0,50000)

erlang <- function(x,rate,k,t) {
  return(rate^k*x^(k-1)*exp(-rate*x)/factorial(k-1))
}

ggplot() +
  geom_function(fun = erlang, args = list(rate = rateth, k = 20, t = t)) +
  xlim(0,2000000)


# TIME EVOLUTION

# open rds file
final_time1 <- readRDS("./simulations_time/output_[10_20_01][1.0].rds")
final_time2 <- readRDS("./simulations_time/output_[10_20_005][1.0].rds")
final_time3 <- readRDS("./simulations_time/output_[10_20_001][1.0].rds")

patch1 <- Zt_log_plots(final_time1) + Zt_log_plots(final_time2) + Zt_log_plots(final_time3)
patch1

patch2 <- Zt_norm_plots(final_time1) + Zt_norm_plots(final_time2) + Zt_norm_plots(final_time3)
patch2

patch3 <- patch1 / patch2
patch3
ggsave("./imgs/15_10/15_10_timespatch2.png", dpi=600)

# MARGINAL DISTRIBUTIONS: HISTOGRAMS

final1 <- readRDS("./simulations/final[15_10_01]1000:0.6.rds")
final2 <- readRDS("./simulations/final[10_20_005]1000[0.8].rds")
final3 <- readRDS("./simulations/10_20_001_01 [0.8]/final.rds")

# horizontal
plot1 <- final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 6e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 11)) +
  scale_x_continuous(breaks=c(20000,40000)) +
  ylim(0,875)
plot1

plot2 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 5.5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  facet_grid(type~.) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.position = "none") +
  ylim(0,875)
plot2

plot3 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 5.5e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(20000,40000))
plot3

patch4 <- plot1 + plot2 + plot3
patch4

# vertical
plot1 <-final %>% 
  #filter(`Z-` < 8e5) %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=50, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  facet_grid(~type) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.position = "none") #+
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  #scale_x_continuous(breaks = c(0,20000,40000)) +
  #ylim(0,1000)
plot1



plot2 <-final2 %>% 
  #filter(`Z-` < 8e5) %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=50, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  facet_grid(~type) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.position = "none") +
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(breaks = c(0,20000,40000)) +
  ylim(0,1000)
  #xlim(0,8e5)
plot2

plot3 <-final3 %>% 
  #filter(`Z-` < 8e5) %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=50, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  facet_grid(~type) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.position = "none") +
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(breaks = c(0,20000,40000))
plot3

options(scipen = 0)
patch4 <- plot1 + plot2 + plot3
patch4

patch3/patch4
ggsave("./imgs/15_10/15_10_partial.png", dpi=600)

# FITS

# for exponential distributions

rate1 <- Pz_exp(final1,alpha_min,alpha_plus,0.01,0.1)
rate2 <- Pz_exp(final2,alpha_min,alpha_plus,0.05,0.05)
rate3 <- Pz_exp(final3,alpha_min,alpha_plus,0.1,0.01)

final_pz_plots <- Pz_exp_plot(final,rate)
final_pz_plots[1]
final_pz_plots[2]

final1 %>% 
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
plot1a <- final1 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 90, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[1])), linewidth=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = c(0,2e4,4e4))
plot1a

plot2a <- final2 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)),fill = rgb(102,204,102,maxColorValue = 255), bins = 90, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[1])), linewidth=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = c(0,2e4,4e4))
plot2a

plot3a <- final3 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 90, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[1])), linewidth=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = c(0,2e4,4e4))
plot3a

patch5 <- plot1a + plot2a + plot3a
patch5

# Z+

plot1b <- final1 %>% 
  filter(`Z+` < 1000) %>% 
  ggplot(aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[2])), size=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  ylim(0,0.015) + xlim(0,1000)
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-3)) +
  #scale_x_continuous(breaks = c(0,20000,40000),labels = function(x) format(x, scientific = TRUE),limits = c(0,4.5e4))
plot1b

plot2b <- final2 %>% 
  filter(`Z+` < 1000) %>% 
  ggplot(aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[2])), size=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  ylim(0,0.015) + xlim(0,1000)
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0,1e-3)) +
  #scale_x_continuous(breaks = c(0,20000,40000), labels = function(x) format(x, scientific = TRUE),limits = c(0,4.5e4))
plot2b

plot3b <- ggplot(final3,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 120, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[2])), size=0.6) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  ylim(0,0.015) +
  xlim(0,1000)
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  #scale_x_continuous(breaks = c(0,20000,40000),labels = function(x) format(x, scientific = TRUE),limits = c(0,4.5e4))
plot3b

patch6 <- plot1b + plot2b + plot3b
patch6

patch1 / patch2 / patch4 / patch5 / patch6
ggsave("./imgs/15_10/15_10_full.png", dpi=600)
# when not exponential

# Z+ follows a POWER-LAW

time1 <- mean(final1$time)
time2 <- mean(final2$time)
time3 <- mean(final3$time)

powerlaw_coeff <- function(alpha_min,alpha_plus,lambda_min,lambda_plus,omega_m,omega_p,time) {
  return(alpha_min/(lambda_plus*alpha_plus)*omega_p*pi/sin(pi*lambda_min/lambda_plus)*exp((lambda_min+omega_m*omega_p/(lambda_plus-lambda_min)*lambda_min/lambda_plus)*time)/gamma(1-lambda_min/lambda_plus))
}

powerlaw_coeff1 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus, omega_m = 0.01,omega_p = 0.1,time = time1)
powerlaw_coeff2 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus,omega_m = 0.05,omega_p = 0.05,time = time2)
powerlaw_coeff3 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus,omega_m = 0.1,omega_p = 0.01,time = time3)

powerlaw <- function(x,coeff,exponent) {
  return(coeff*x^exponent)
}

plot1p <- final1 %>% 
  filter(`Z+`<4e5) %>% 
  mutate(pw=powerlaw_coeff1*`Z+`**(-(1+lambda_min/lambda_plus))) %>%
  ggplot(aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 80, na.rm = TRUE) +
  geom_line(aes(y=pw), color='black') +
  ylim(0,2e-5) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,2.5e5),labels = function(x) format(x, scientific = TRUE),breaks = c(1e5,2e5))
  #stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff1, exponent = -(1+lambda_min/lambda_plus)), size=0.5) #+
  #ylim(0,8e-4)
plot1p 

plot2p <- final2 %>% 
  filter(`Z+`<4e5) %>% 
  mutate(pw=powerlaw_coeff2*`Z+`**(-(1+lambda_min/lambda_plus))) %>%
  ggplot(aes(x=`Z+`)) + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 100, na.rm = TRUE) +
  geom_line(aes(y=pw), color='black') + ylim(0,2e-5) +
  scale_x_continuous(limits = c(0,2.7e5),labels = function(x) format(x, scientific = TRUE),breaks = c(1e5,2e5))
  #stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff2, exponent = -(1+lambda_min/lambda_plus)), size=0.5) #+
  #xlim(0,8e5) + ylim(0,2.5e-5)
plot2p

plot3p <- final3 %>% 
  filter(`Z+`<4e5) %>% 
  mutate(pw=powerlaw_coeff3*`Z+`**(-(1+lambda_min/lambda_plus))) %>%
  ggplot(aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 120, na.rm = TRUE) +
  geom_line(aes(y=pw), color='black') + ylim(0,2e-5) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,2.7e5),labels = function(x) format(x, scientific = TRUE),breaks = c(1e5,2e5))
  #stat_function(fun = powerlaw, args = list(coeff = powerlaw_coeff3, exponent = -(1+lambda_min/lambda_plus)), size=0.5) +
  #scale_x_continuous(breaks = c(0,20000,40000,60000),limits = c(0,60000)) +
  #scale_y_continuous(limits = c(0,2e-4),labels = function(x) format(x, scientific = TRUE))
plot3p

patch6 <- plot1p + plot2p + plot3p
patch6


# Z-
# exponential regime
time1 = mean(final1$time)
time2 = mean(final2$time)
time3 = mean(final3$time)

ratem1 <- exp(-alpha_min*time1)
ratem2 <- exp(-alpha_min*time2)
ratem3 <- exp(-alpha_min*time3)

final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  dplyr::mutate(eexp = dexp(Z,rate=ratem1)) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(density), fill=type), bins=180) +
  scale_fill_manual(values=my_palette) +
  geom_line(aes(x=Z,y=eexp))

plot1 <- ggplot(final1,aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 150, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = ratem1), size=0.8)
plot1

plot2 <- ggplot(final2,aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 150, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = ratem2), size=0.8)
plot2

plot3 <- ggplot(final3,aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 150, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = ratem3), size=0.8)
plot3

plot1 / plot2 / plot3

# upload numerical solution
analytic1 <- read.csv("./imgs/10_20/[10_20_01]/P(Z-)[10_20_01_001].csv") %>%
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/10_20/[10_20_01]/P(Z+)[10_20_01_001][0.8] mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))

analytic2 <- read.csv("./imgs/10_20/[10_20_005]/P(Z-)[10_20_005][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/10_20/[10_20_005]/P(Z+)[10_20_005][0.8]2 mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))

analytic3 <- read.csv("./imgs/10_20/[10_20_001]/P(Z-)[10_20_001_01][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/10_20/[10_20_001]/P(Z+)[10_20_001_01][0.8] mathematica.csv") %>% 
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
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  xlim(0,45000)
plot4

plot5 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-" & Z<40000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255)) +
  geom_line(data=analytic2 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15))#+
#xlim(0,40000)
plot5

plot6 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255), na.rm = TRUE) +
  geom_line(data=analytic3 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15)) +
  xlim(0,45000)
plot6

patch5 <- plot4 + plot5 + plot6
patch5

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

patch6 <- plot7 + plot8 + plot9
patch6

patch1 / patch2 / patch4 / patch5 / patch6
ggsave("./imgs/10_20/10_20_full.png", dpi=600, width = 12.4, height = 11.7, units = c("in"))


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

# ANALYTICAL PLOTS

# Marginal individual

# first
alpha_min = 15; lambda_min = 15; alpha_plus = 13; lambda_plus = 13; omega_p = 0.1; omega_m = 0.01
t = 0.6

# minus
ratem9 = lambda_min/alpha_min*exp(-(lambda_min+omega_m*omega_p/(lambda_min - lambda_plus))*t)
# plus
ratep9 = ratem9*(lambda_min - lambda_plus)/omega_p

# 1 2 4 5
t11 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratem1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratem2),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratem4),aes(color="3")) +
  geom_function(fun = dexp, args = list(rate=ratem5),aes(color="4")) +
  scale_color_manual(values=c("1","2","3","4"),labels = c("15_10","20_10","18_10","16_10")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z-") + ylab("P(Z-)") +
  theme(legend.position = "none") +
  xlim(0,300000)

t12 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratep1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratep2),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratep4),aes(color="3")) +
  geom_function(fun = dexp, args = list(rate=ratep5),aes(color="4")) +
  scale_color_manual(values=c("1","2","3","4"),labels = c("15_10","20_10","18_10","16_10")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z+") + ylab("P(Z+)") +
  theme(legend.position = "none") +
  xlim(0,7000)

t11 + t12

# 1 6 7
t13 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratem1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratem6),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratem7),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("01_001","001_01","005_005")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z-") + ylab("P(Z-)") +
  xlim(0,20000)

t14 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratep1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratep6),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratep7),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("01_001","001_01","005_005")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z+") + ylab("P(Z+)") +
  xlim(0,2000)

t13 + t14

# 1 8 9
t15 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratem1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratem8),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratem9),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("15_10","15_8","15_13")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z-") + ylab("P(Z-)") +
  xlim(0,20000)

t16 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratep1),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratep8),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratep9),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("15_10","15_8","15_13")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Z+") + ylab("P(Z+)") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  theme(legend.position = "none") +
  xlim(0,3000)

t15 + t16

(t11 + t12)/(t13 + t14)/(t15 + t16)

ggsave("./imgs/theorem1/case1.png",dpi = 600)

# second
alpha = 20; lambda = 20; omega_p = 0.05; omega_m = 0.05
t = 0.6

# minus
ratem55 = lambda/(alpha*cosh(sqrt(omega_m*omega_p)*t))*exp(-lambda*t)
# plus
ratep55 = ratem55*sqrt(omega_m/omega_p)/tanh(sqrt(omega_m*omega_p)*t)

# 11 22 33
t21 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratem11),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratem22),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratem33),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("20_20","15_15","10_10")) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Z-") + ylab("P(Z-)") +
  theme(legend.position = "none") +
  xlim(0,10000)
ggsave("./imgs/theorem1/legend1.png",dpi=600)

t22 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratep11),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratep22),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratep33),aes(color="3")) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_color_manual(values=c("1","2","3"),labels = c("20_20","15_15","10_10")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Z+") + ylab("P(Z+)") +
  theme(legend.position = "none") +
  xlim(0,1000)

t21 + t22

# 11 44 55

t23 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratem11),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratem44),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratem55),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("01_001","001_01","005_005")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  xlab("Z-") + ylab("P(Z-)") +
  xlim(0,800000)

t24 <- ggplot() +
  geom_function(fun = dexp, args = list(rate=ratep11),aes(color="1")) +
  geom_function(fun = dexp, args = list(rate=ratep44),aes(color="2")) +
  geom_function(fun = dexp, args = list(rate=ratep55),aes(color="3")) +
  scale_color_manual(values=c("1","2","3"),labels = c("01_001","001_01","005_005")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  theme(legend.position = "none") +
  xlab("Z+") + ylab("P(Z+)") +
  xlim(0,20000)

t23 + t24

(t21+t22)/(t23+t24)
ggsave("./imgs/theorem1/case2.png",dpi=600)

# JOINT PLOT

lambda_min=10; lambda_plus=20; omega_m=0.01; omega_p=0.1; alpha_min=10; alpha_plus=20; beta_min=0; beta_min=0
final <- readRDS("./simulations/10_20_01_001 [0.8]/final.rds")
final <- final %>% filter(final$`Z+`>1)
final <- final %>% filter(final$`Z-`>1)
final <- final %>% filter(final$`Z+`<5e6)
#final <- final %>% filter(final$`Z-` != 265)

ggplot() +
  geom_point(final,mapping = aes(x=`Z+`, y=`Z-`),geom = "point",bins = 70) #+
  #stat_density_2d(joint_th, mapping = aes(x=x,y=y, fill = ..level..), geom = "polygon", alpha = 0.5)
  #scale_fill_continuous(type = "viridis") 
  #scale_fill_brewer()
  #theme_bw()

ggplot() +
  geom_point(final,mapping = aes(x=`Z+`, y=`Z-`)) +
  #stat_density_2d_filled(joint_th, mapping = aes(x=x,y=y), alpha = 0.4) +
  geom_density_2d(joint_th, mapping = aes(x=x,y=y),colour = "black") +
  #scale_x_continuous(trans = "log10", limits = c(1,1.2e6), breaks = c(1e2,1e4,1e6)) +
  #scale_y_continuous(trans = "log10", limits = c(10,1.2e6), breaks = c(1e2,1e4,1e6)) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
  #xlim(0,10e3)
ggsave("./imgs/pt/joint_15_10_005.png",dpi=600)

ggplot() +
  stat_density_2d_filled(joint_th, mapping = aes(x=x,y=y),alpha = 0.5) +
  scale_fill_brewer()

x <- rpareto(1000,lambda_min/lambda_plus, scale = 10000)
y <- sapply(x, function(x) x*omega_m/(lambda_plus - lambda_min))

joint_th <- data.frame(x,y)

ggplot() +
  stat_density_2d(aes(x=rexp(1000,rate = (lambda_min-lambda_plus)/omega_p*exp(-(lambda_min+omega_m*omega_p/(lambda_min-lambda_plus))*0.8)),y=rexp(1000,rate = exp(-(lambda_min+omega_m*omega_p/(lambda_min-lambda_plus))*0.8)),fill= ..level..))
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")

ggplot(joint_th,aes(x,y)) +
  geom_density_2d(alpha = 0.5)


xmin = 1000; alpha = 1.5

n = 1000
x = rplcon(n, xmin, alpha)
con_rns = sort(con_rns)
p = rep(1/n, n)
#Zipfs plot
plot(con_rns, rev(cumsum(p)), log="xy", type="l")

n = 1e5
x = rpldis(n, xmin, alpha)

# P(Z = Z- + Z+)

time <- mean(final$time)


powerlaw_coeff <- function(alpha_min,alpha_plus,lambda_min,lambda_plus,omega_m,omega_p,time) {
  return(alpha_min/(lambda_plus*alpha_plus)*omega_p*pi/sin(pi*lambda_min/lambda_plus)*exp((lambda_min+omega_m*omega_p/(lambda_plus-lambda_min)*lambda_min/lambda_plus)*time)/gamma(1-lambda_min/lambda_plus))
}

powerlaw_coeff1 <- powerlaw_coeff(alpha_min = alpha_min,alpha_plus = alpha_plus,lambda_min = lambda_min, lambda_plus = lambda_plus, omega_m = 0.01,omega_p = 0.1,time = time)

powerlaw <- function(x,coeff,exponent) {
  return(coeff*x^exponent)
}
ratesum <- (lambda_min - lambda_plus)/(lambda_min - lambda_plus + omega_p)*exp(-(lambda_min+omega_m*omega_p/(lambda_min-lambda_plus))*time)


final %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`+`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#9C5DA6", bins = 90, na.rm=TRUE) +
  stat_function(fun = dexp, args = list(rate = ratesum), linewidth=1) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Z") + ylab("P(Z)")

ggsave("./imgs/pt/sum.png",dpi=600)
