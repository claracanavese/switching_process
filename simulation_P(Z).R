library(LaplacesDemon)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(reshape2)
library(magrittr)

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
  #add columns with sum and proportions
  dat <- mutate(dat, sum = `Z+`+`Z-`, "Z- ratio" = `Z-`/sum, "Z+ ratio" = `Z+`/sum)
  
  # plot the evolution of both populations on same graph (log scale)
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
    scale_y_continuous(trans = 'log10') +
    theme(axis.title=element_text(size=20),
          axis.text = element_text(size=14),
          legend.position = "none")
  
  # plot relative proportions
  plot2 <- dat %>% 
    dplyr::select(t, dplyr::contains("ratio")) %>% 
    melt(id="t", variable.name="type", value.name="ratio") %>% 
    ggplot() +
    geom_line(aes(x=t, y=ratio, color=type),size=2) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) +
    theme(axis.title=element_text(size=20),
          axis.text = element_text(size=14),
          legend.position = "none") +
    ylab(bquote(Z/Z[tot]))
  output <- list(plot1,plot2)
  return(output)
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
alpha_min = 20;beta_min = 0;alpha_plus = 20;beta_plus = 0;omega_p = 0.1;omega_m = 0.1

options(scipen = 1)

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
final <- readRDS("./simulations/20_20_01_001/final.rds")
final_time <- readRDS("./simulations_time/output_[20_20_001][0.8].rds")

# TIME EVOLUTION

final_time_plots <- Zt_plots(final_time)
final_time_plots[1]
final_time_plots[2]

# MARGINAL DISTRIBUTIONS: HISTOGRAMS

# horizontal
final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter((type == "Z-" & Z < 3.8e4) | type == "Z+") %>% 
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type), color = 'black',bins=200) +
  scale_fill_manual(values=my_palette) +
  facet_grid(type~.) +
  theme(axis.title=element_text(size=20),
      axis.text = element_text(size=14),
      legend.position = "none")

# vertical
final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter((type == "Z-" & Z < 8.25e5) | type == "Z+") %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(count), fill=type),color = "black",bins=100) +
  scale_fill_manual(values=my_palette) +
  facet_grid(~type) +
  theme(axis.title=element_text(size=20),
        axis.text = element_text(size=14),
        legend.position = "none") +
  scale_x_continuous(breaks=c(200000,400000,600000,800000),limits = c(0,820000))

# FITS

# for exponential distributions only
rate <- Pz_exp(final,alpha_min,alpha_plus,omega_m,omega_p)

final_pz_plots <- Pz_exp_plot(final,rate)
final_pz_plots[1]
final_pz_plots[2]

final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  #filter((type == "Z-" & Z < 3e5) | type == "Z+") %>% 
  dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate=as.numeric(rate[1])), dexp(Z,rate=as.numeric(rate[2])))) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(density), fill=type), color="black", bins=180) +
  scale_fill_manual(values=my_palette) +
  geom_line(aes(x=Z,y=rate)) +
  facet_grid(~type) 

# individual plots
final %>% 
  filter(`Z-` < 1e6) %>% 
  ggplot(aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), color = "black", fill = rgb(102,204,102,maxColorValue = 255), bins = 120) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate[1])), size=0.8) +
  theme(axis.title=element_text(size=20),
        axis.text = element_text(size=14)) +
  xlim(0,850000) + ylim(0,8e-6)


ggplot(final,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 120) +
  stat_function(fun = dexp, args = list(rate = as.numeric(rate[2])), size=0.8) +
  theme(axis.title=element_text(size=20),
        axis.text = element_text(size=14)) +
  ylim(0,0.0002) +
  scale_x_continuous(breaks=c(30000,60000,90000),limits = c(0,90000))


# when not exponential, upload analytic distribution for Z-
analytic <- read.csv("./imgs/[10_20_005]/P(Z-)[10_20_005][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/[10_20_005]/P(Z+)[10_20_005][0.8] mathematica.csv") %>% 
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
final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z-") %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=150, color="black",fill=rgb(102,204,102,maxColorValue = 255)) +
  geom_line(data=analytic %>% filter(type=="Z-"), aes(x=Z,y=P), size=1)

# Z+
final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  filter(type == "Z+" & Z < 823432) %>% 
  ggplot() +
  geom_histogram(aes(x=Z,y=after_stat(density)), bins=150, color="black",fill="#D5D139") +
  geom_line(data=analytic %>% filter(type=="Z+" & Z < 823432), aes(x=Z,y=P)) +
  xlim(0,823433) +
  ylim(0,3e-5)


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

