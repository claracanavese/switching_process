library(LaplacesDemon)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(reshape2)
library(patchwork)

simulation <- function(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, index) {
  # initialize Z and t
  Z_minus = 1; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)
  # tibble to store evolution of populations with time
  # Zt_output <- tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]) 
  while (t < 0.8) {
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

Zt_plots <- function(dat) {

  # plot the evolution of both populations on same graph (log scale)
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point(size=0.8) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
    scale_y_continuous(trans = 'log10') +
    theme(legend.position = "none") 

  return(plot1)
}
Zt_norm_plots <- function(dat) {
  #add columns with sum and proportions
  dat <- mutate(dat, sum = `Z+`+`Z-`, "Z- ratio" = `Z-`/sum, "Z+ ratio" = `Z+`/sum)
  
  # plot relative proportions
  plot2 <- dat %>% 
    dplyr::select(t, dplyr::contains("ratio")) %>% 
    melt(id="t", variable.name="type", value.name="ratio") %>% 
    ggplot() +
    geom_line(aes(x=t, y=ratio, color=type),size=0.8) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) +
    theme(legend.position = "none") +
    ylab(bquote(Z/Z[tot]))
  
  return(plot2)
}
Pz_exp <- function(dat,alpha_min,alpha_plus,omega_m,omega_p) {
  time <- mean(dat$time)
  if (alpha_min > alpha_plus) {
    exprate_minus <- exp(-(alpha_min+omega_m*omega_p/(alpha_min/alpha_plus))*time)
    exprate_plus <- (alpha_min - alpha_plus)/omega_p*exp(-(alpha_min+omega_m*omega_p/(alpha_min/alpha_plus))*time)
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

# open rds file

final_time1 <- readRDS("./simulations_time/output_[10_20_01][1.0].rds")
final_time2 <- readRDS("./simulations_time/output_[10_20_005][1.0].rds")
final_time3 <- readRDS("./simulations_time/output_[10_20_001][1.0].rds")

# TIME EVOLUTION

patch1a <- Zt_plots(final_time1) + Zt_plots(final_time2) + Zt_plots(final_time3)
patch1b <- Zt_norm_plots(final_time1) + Zt_norm_plots(final_time2) + Zt_norm_plots(final_time3)
patch1a / patch1b

# MARGINAL DISTRIBUTIONS: HISTOGRAMS

final1 <- readRDS("./simulations/10_20_01_001 [0.8]/final.rds")
final2 <- readRDS("./simulations/final[10_20_005]1000[0.8].rds")
final3 <- readRDS("./simulations/10_20_001_01 [0.8]/final.rds")
                                  

plot1 <- final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 6e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=120, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") #+
  #xlim(0,5e6)
  #scale_x_continuous(breaks=c(4e5,8e5), limits = c(0,8.2e5)) +
  #ylim(0,400)
plot1

plot2 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 6e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") #+
  #scale_x_continuous(breaks=c(4e5,8e5))
plot2

plot3 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter( Z < 6e4) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),bins=60, na.rm = TRUE) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(legend.position = "none") #+
  #scale_x_continuous(breaks=c(4e5,8e5))
plot3

patch1 <- plot1 + plot2 + plot3
patch1

# FITS

# for exponential distributions only
rate1 <- Pz_exp(final1,alpha_min,alpha_plus,omega_m,omega_p)
rate2 <- Pz_exp(final2,alpha_min,alpha_plus,0.05,0.05)
rate3 <- Pz_exp(final3,alpha_min,alpha_plus,0.1,0.01)

final_pz_plots <- Pz_exp_plot(final,rate)
final_pz_plots[1]
final_pz_plots[2]

final1 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  #filter((type == "Z-" & Z < 3e5) | type == "Z+") %>% 
  dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate=as.numeric(rate1[1])), dexp(Z,rate=as.numeric(rate1[2])))) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(density), fill=type), color="black", bins=180) +
  scale_fill_manual(values=my_palette) +
  geom_line(aes(x=Z,y=rate)) +
  facet_grid(~type) 

# individual plots
# Z-
plot4 <- final1 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[1])), linewidth=0.5) +
  scale_x_continuous(breaks = c(5e5,1e6), limits = c(0,1.4e6))
  #+ ylim(0,0.00015)
plot4

plot5 <- final2 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[1])), linewidth=0.5) #+
  #xlim(0,60500) + ylim(0,0.00015)
plot5

plot6 <- final3 %>% 
  #filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = rgb(102,204,102,maxColorValue = 255), bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[1])), linewidth=0.5) +
  scale_x_continuous(breaks = c(5e5,1e6), limits = c(0,1.4e6))
plot6

patch2 <- plot4 + plot5 + plot6
patch2

# Z+
plot7 <- final1 %>% 
  filter(`Z+`<5e4) %>% 
  ggplot(aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate1[2])), size=0.5, na.rm = TRUE) +
  xlim(0,5e4) + ylim(0,0.00065)
plot7

plot8 <- ggplot(final2,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate2[2])), size=0.5, na.rm = TRUE) +
  xlim(0,5e4) + ylim(0,0.00065)
plot8

plot9 <- ggplot(final3,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)), fill = "#D5D139", bins = 80, na.rm = TRUE) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate3[2])), size=0.5, na.rm = TRUE) +
  xlim(0,5e4) + ylim(0,0.00065)
plot9

patch3 <- plot7 + plot8 + plot9
patch3

patch1 / patch2 / patch3

# when not exponential, upload analytic distribution for Z-
analytic1 <- read.csv("./imgs/[10_20_01]/P(Z-)[10_20_01_001][0.8] mathematica.csv") %>% 
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
  geom_line(data=analytic1 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) +
  xlim(0,45000)
plot1

plot5 <- final2 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-" & Z<40000) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255)) +
  geom_line(data=analytic2 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) #+
  #xlim(0,40000)
plot2

plot6 <- final3 %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=80, fill=rgb(102,204,102,maxColorValue = 255), na.rm = TRUE) +
  geom_line(data=analytic3 %>% filter(type=="Z-"), aes(x=Z,y=P), size=0.5) +
  xlim(0,45000)
plot3

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

# plot analytic distribution over the histograms
Pz_noexp_plot <- function(dat, analytic) {
  plot1 <- dat %>% 
    reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
    ggplot() +
    geom_histogram(aes(x=Z, y=after_stat(count), fill=type), color="black", bins=180) +
    scale_fill_manual(values=my_palette) +
    facet_grid(type~., scales="free_y")
  
  plot2 <- dat %>% 
    reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
    filter((type=="Z+" & Z < 80000) | (type=="Z-")) %>% 
    ggplot() +
    geom_histogram(aes(x=Z,fill=type, y=after_stat(density)),
                   color = "black", bins=100, position = "identity", alpha=0.5) +
    # xlim(0,80000) +
    geom_line(data=analytic %>% 
                dplyr::group_by(type) %>% 
                dplyr::mutate(nn=dplyr::n()), aes(x=Z,y=P)) + #+
    facet_grid(~type, scales="free")
}


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
# plot Z- VS Z+
final %>% 
  mutate(sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum) %>% 
  ggplot(aes(x=`Z-`,y=`Z+`)) + 
  geom_point(color = 'cyan4') #+
  #ylim(0,max(max(final$`Z-`),max(final$`Z+`)))
  
final %>% 
  mutate(sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum) %>% 
  ggplot(aes(x=`Z- rate`,y=`Z+ rate`)) + 
  geom_point(color = 'cyan4') +
  ylim(0,1)

rm(list = ls())

