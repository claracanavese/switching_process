library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)

# define parameters
alpha_min = 10;beta_min = 0;alpha_plus = 20;beta_plus = 0;omega_p = 0.1;omega_m = 0.01

# population starting with 1 cell in state -
Z_minus = 1; Z_plus = 0; t = 0
Z <- c(Z_minus,Z_plus)

# define stoichiometric vectors
o1 <- c(1,0);o2 <- c(-1,0);o3 <- c(-1,1);o4 <- c(0,1);o5 <- c(0,-1);o6 <- c(1,-1)
o <- rbind(o1,o2,o3,o4,o5,o6)

# create tibble to store Z values for each t
output <- tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2])

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
  print(t)
  i <- rcat(1,anorm)
  Z <- Z+o[i,]
  output <- bind_rows(output,tibble("t" = t,"Z-" = Z[1],"Z+" = Z[2]))
}

# save data in rds file
saveRDS(output, file = paste0("./simulations_time/output_",alpha_min,"_",alpha_plus,"_",omega_p,"_",omega_m,".rds"))


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
final <- Zt_plots(output)
final[1]
final[2]




# MULLER
# edges <- data.frame(Parent = paste0(rep("Z.",2), LETTERS[1:2]),Identity = paste0(rep("Z.",2), LETTERS[2:1]))
edges <- data.frame(Parent = c("Z.0","Z.A"), Identity = c("Z.A","Z.B"))
edges
colnames(dat1) <- c("Time","Identity","Population")
dat1 = rbind(dat1, data.frame("Time"=0,"Identity"="Z.0","Population"=1))

# View(dat1)
Muller_df <- get_Muller_df(edges,dat1)
Muller_plot(Muller_df, add_legend = TRUE, xlab = "Time", ylab = "Proportion")
dat1


ggplot(Muller_df, aes_string(x = "Time", y = "Frequency", fill = "Identity", colour = "Identity")) +
  geom_area() + theme(legend.position = "right") + 
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = 25 * (0:4), name = "Percentage")

rm(list=ls())
