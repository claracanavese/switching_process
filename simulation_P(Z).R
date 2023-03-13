library(LaplacesDemon)
library(ggplot2)
library(magrittr)
library(tidyverse)

simulation <- function(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, index) {
  # initialize Z and t
  Z_minus = 1; Z_plus = 0; t = 0
  Z <- c(Z_minus,Z_plus)
  # define stoichiometric vectors
  o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
  o <- rbind(o1,o2,o3,o4,o5,o6)
  # tibble to store evolution of populations with time
  Zt_output <- tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]) 
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
    Zt_output <- bind_rows(Zt_output,tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]))
  }
  if (index == 2) {
    # store Z(t)
    saveRDS(Zt_output, file = paste0("./simulations_time/output_",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,".rds"))
  }
  # return final populations
  output = tibble("Z-" = Z[1],"Z+" = Z[2], "time" = t)
  return(output)
}
Zt_plots <- function(dat) {
  #add columns with sum and proportions
  dat <- mutate(dat, sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum)
  
  # plot the evolution of both populations on same graph (log scale)
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
    scale_y_continuous(trans = 'log10')
  
  # plot relative proportions
  plot2 <- dat %>% 
    dplyr::select(t, dplyr::contains("rate")) %>% 
    melt(id="t", variable.name="type", value.name="rate") %>% 
    ggplot() +
    geom_line(aes(x=t, y=rate, color=type)) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139"))
  output <- list(plot1,plot2)
  return(output)
}

alpha_min = 20;beta_min = 0;alpha_plus = 10;beta_plus = 0;omega_p = 0.1;omega_m = 0.01

# run simulation 1000 times
output_list = lapply(1:500, function(i){
  x <- simulation(alpha_min, alpha_plus, beta_min, beta_plus, omega_m, omega_p, i)
  x$step <- i
  print(i)
  x
})

final = output_list %>% do.call(rbind, .)

# save file in rds format
saveRDS(final, file = paste0("./simulations/final",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,".rds"))
# open rds file
final <- readRDS("./simulations/final[20_10_01_001] 500:0.5.rds")
final_time <- readRDS("./simulations_time/output_[20_10_01_001][0.5].rds")

# TIME EVOLUTION

final_time_plots <- Zt_plots(final_time)
final_time_plots[1]
final_time_plots[2]

# MARGINAL DISTRIBUTIONS

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
rate <- Pz_exp(final,alpha_min,alpha_plus,omega_m,omega_p)

# plot histograms and, if exponential, the analytic solution
Pz_exp_plot <- function(dat,rate) {
  dat %>% 
    reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
    # dplyr::mutate(rate=ifelse(type=="Z-", exprate_minus, exprate_plus)) %>% 
    dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate=exprate_minus), dexp(Z,rate=exprate_plus))) %>%
    ggplot() +
    geom_histogram(aes(x=Z, y=after_stat(density), fill=type), color="black", bins=100) +
    geom_line(aes(x=Z,y=rate)) +
    # stat_function(aes(x=Z), fun = dexp, args = list(rate = exprate_plus)) +
    facet_grid(type~., scales="free_y") #+ ylim(0,5e-5) + xlim(0,5e5)
}
#aes(y=after_stat(density)),
final %>% 
  reshape2::melt(id=c("time","step"), variable.name="type", value.name="Z") %>%
  # dplyr::mutate(rate=ifelse(type=="Z-", exprate_minus, exprate_plus)) %>% 
  dplyr::mutate(rate=ifelse(type=="Z-", dexp(Z,rate=exprate_minus), dexp(Z,rate=exprate_plus))) %>%
  ggplot() +
  geom_histogram(aes(x=Z, y=after_stat(density), fill=type), color="black", bins=100) +
  geom_line(aes(x=Z,y=rate)) +
  # stat_function(aes(x=Z), fun = dexp, args = list(rate = exprate_plus)) +
  facet_grid(type~., scales="free_y") #+ ylim(0,5e-5) + xlim(0,5e5)
  
ggplot(final,aes(x=`Z-`)) + 
  geom_histogram(aes(y=after_stat(density)), color = "black", fill = rgb(102,204,102,maxColorValue = 255), bins = 40) +
  stat_function(fun = dexp, args = list(rate = exprate_minus))

ggplot(final,aes(x=`Z+`)) + 
  geom_histogram(aes(y=after_stat(density)),color = "black", fill = "#D5D139", bins = 100) + #+ xlim(0,78000) + ylim(0,60)
  stat_function(fun = dexp, args = list(rate = exprate_plus))

# when not exponential, upload analytic distribution for Z-
analytic <- read.csv("./imgs/[10_20_01_001]/P(Z-)[10_20_01_001][0.8] mathematica.csv") %>% 
  tibble::as_tibble() %>% 
  mutate(type="Z-") %>% 
  add_row(read.csv("./imgs/[10_20_01_001]/P(Z+)[10_20_01_001][0.8] mathematica.csv") %>% 
            tibble::as_tibble() %>% 
            mutate(type="Z+"))


# plot analytic distribution over the histograms
final %>% 
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


# final %>% top_n(5,`Z-`)
# final %>% top_n(5,`Z+`)

# plot histograms on same graph
p <- ggplot() +
  geom_histogram(aes(x = final$`Z-`, fill = "Z+"), alpha = 0.5) +
  geom_histogram(aes(x = final$`Z+`, fill = "Z-"), alpha = 0.5) +
  scale_fill_manual(values = c("Z-" = "red", "Z+" = "green")) +
  xlab("Z") +
  ylab("count") +
  xlim(0,80000) +
  ylim(0,250)
p

# plot Z- VS Z+
final <- final %>% mutate(sum = `Z+`+`Z-`) %>% mutate("Z- rate" = `Z-`/sum) %>% mutate("Z+ rate" = `Z+`/sum)
final %>% ggplot(aes(x=`Z-`,y=`Z+`)) + geom_point(color = 'cyan4')
final %>% ggplot(aes(x=`Z- rate`,y=`Z+ rate`)) + geom_point(color = 'cyan4')

rm(list = ls())

