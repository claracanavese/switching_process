library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)
library(ggmuller)
library(patchwork)

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

while (t < 1.2) { 
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
  print(t)
  i <- rcat(1,anorm)
  Z <- Z+o[i,]
  output <- bind_rows(output,tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]))
}

# save data in rds file
saveRDS(output, file = paste0("./simulations_time/10.1_10/sim1_[1_1.1][1.2]"))


Zt_plots <- function(dat) {
  #add columns with sum and proportions
  dat <- mutate(dat, sum = `Z+`+`Z-`, "Z- rate" = `Z-`/sum, "Z+ rate" = `Z+`/sum)
  
  # plot the evolution of both populations on same graph 
  dat1 <- melt(dat[,1:3], value.name = "Z", id = "t")
  plot1 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) #+ 
    #scale_y_continuous(trans = 'log10')

  plot2 <- dat1 %>% ggplot(aes(t,Z, col=variable)) + ylab("Z") + geom_point() +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10')
  
  # plot relative proportions
  plot3 <- dat %>% 
    dplyr::select(t, dplyr::contains("rate")) %>% 
    melt(id="t", variable.name="type", value.name="rate") %>% 
    ggplot() +
    geom_line(aes(x=t, y=rate, color=type)) +
    scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139"))
  output <- list(plot1,plot2,plot3)
  return(output)
}
final <- Zt_plots(final_time)
final[1]
final[2]
final[3]

final_time <- readRDS("./simulations_time/output_[10_20_005][1.05].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot3 <- final_time1 %>% 
  mutate(eexp1=omega_p/alpha_plus*pi/sin(pi*alpha_min/alpha_plus)*omega_m/(alpha_plus-alpha_min)*exp((alpha_plus+omega_m*omega_p/(alpha_plus-alpha_min))*t)+exp(alpha_min*t)) %>% 
  mutate(eexp2=omega_p/alpha_plus*pi/sin(pi*alpha_min/alpha_plus)*exp((alpha_plus+omega_m*omega_p/(alpha_plus-alpha_min))*(t-0.1))) %>% 
  ggplot(aes(t, col=variable)) + 
  ylab("Z") + 
  geom_point(aes(y=Z),size=0.5) +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10') +
  geom_line(aes(y=eexp1), color='black') +
  geom_line(aes(y=eexp2), color='black')
plot3  

final_time <- readRDS("./simulations_time/output_[20_20_01][0.673].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot2 <- final_time1 %>% 
  filter(Z<1e5) %>% 
  mutate(eexp1=cosh(sqrt(omega_m*omega_p)*t)*exp(alpha_min*t)) %>% 
  mutate(eexp2=sinh(sqrt(omega_m*omega_p)*t)*exp(alpha_min*t)*sqrt(omega_p/omega_m)) %>% 
  ggplot(aes(t, col=variable)) + 
  ylab("Z") + 
  geom_point(aes(y=Z),size=0.5) +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10') +
  geom_line(aes(y=eexp1), color='black') +
  geom_line(aes(y=eexp2), color='black') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(bquote(~ lambda['-']==20 ~ lambda['+']==20)) +
  ylim(0,1e5)
plot2

final_time <- readRDS("./simulations_time/output_[15_10_001][0.983].rds")
final_time1 <- melt(final_time[,1:3], value.name = "Z", id = "t")
plot1 <- final_time1 %>% 
  filter(Z<1e5) %>% 
  mutate(eexp1=exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*t)) %>% 
  mutate(eexp2=omega_p/(alpha_min-alpha_plus)*exp((alpha_min+omega_m*omega_p/(alpha_min-alpha_plus))*t)) %>% 
  ggplot(aes(t, col=variable)) + 
  ylab("Z") + 
  geom_point(aes(y=Z),size=0.5) +
  scale_colour_manual(values=c(rgb(102,204,102,maxColorValue = 255),"#D5D139")) + 
  scale_y_continuous(trans = 'log10') +
  geom_line(aes(y=eexp1), color='black') +
  geom_line(aes(y=eexp2), color='black') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(bquote(~ lambda['-']==15 ~ lambda['+']==10)) +
  ylim(0,1e5)
plot1

plot1 + plot2

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
