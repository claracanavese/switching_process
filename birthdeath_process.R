library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)
library(rainbow)
library(viridis)

# define parameters
alpha = 15;beta = 10;
# population starting with 1 cell
Z = 1; t = 0
# create tibble to store Z values for each t
output <- tibble("t" = t,"Z" = Z)

while (t < 1.0) { 
  a1 = alpha*Z
  a2 = beta*Z
  a <- c(a1,a2)
  a0 <- sum(a)
  anorm <- a/a0
  tau <- rexp(n = 1, rate = a0)
  t <- t + tau
  i <- rcat(1,anorm)
  if (i==1) {
    Z = Z+1
  } else if (i==2) {
    Z = Z-1
  }
  print(t)
  output <- bind_rows(output,tibble("t" = t,"Z" = Z))
}

output <- output %>% 
  mutate(eexp=exp((alpha-beta)*t))

# save data in rds file
saveRDS(output, file = paste0("./simulations_time.rds"))

# open rds
output1 <- readRDS("./simulations_time/birthdeath/birthdeath_1.rds")
output2 <- readRDS("./simulations_time/birthdeath/birthdeath_2.rds")
output3 <- readRDS("./simulations_time/birthdeath/birthdeath_3.rds")
output4 <- readRDS("./simulations_time/birthdeath/birthdeath_4.rds")
output5 <- readRDS("./simulations_time/birthdeath/birthdeath_5.rds")
output6 <- readRDS("./simulations_time/birthdeath/birthdeath_6.rds")
output7 <- readRDS("./simulations_time/birthdeath/birthdeath_7.rds")
output8 <- readRDS("./simulations_time/birthdeath/birthdeath_8.rds")
output9 <- readRDS("./simulations_time/birthdeath/birthdeath_9.rds")
output10 <- readRDS("./simulations_time/birthdeath/birthdeath_10.rds")
output11 <- readRDS("./simulations_time/birthdeath/birthdeath_11.rds")
output12 <- readRDS("./simulations_time/birthdeath/birthdeath_12.rds")
output13 <- readRDS("./simulations_time/birthdeath/birthdeath_13.rds")
output14 <- readRDS("./simulations_time/birthdeath/birthdeath_14.rds")
output15 <- readRDS("./simulations_time/birthdeath/birthdeath_15.rds")
output16 <- readRDS("./simulations_time/birthdeath/birthdeath_16.rds")
output17 <- readRDS("./simulations_time/birthdeath/birthdeath_17.rds")

final_plot <- ggplot() +
  geom_line(data=output1,aes(x=t,y=Z,color="1")) +
  geom_line(data=output2,aes(x=t,y=Z,color="2")) +
  geom_line(data=output3,aes(x=t,y=Z,color="3")) +
  geom_line(data=output4,aes(x=t,y=Z,color="4")) +
  geom_line(data=output5,aes(x=t,y=Z,color="5")) +
  geom_line(data=output6,aes(x=t,y=Z,color="6")) +
  geom_line(data=output7,aes(x=t,y=Z,color="7")) +
  #geom_line(data=output8,aes(x=t,y=Z,color="8")) +
  geom_line(data=output9,aes(x=t,y=Z,color="9")) +
  geom_line(data=output10,aes(x=t,y=Z,color="10")) +
  geom_line(data=output11,aes(x=t,y=Z,color="11")) +
  geom_line(data=output12,aes(x=t,y=Z,color="12")) +
  geom_line(data=output13,aes(x=t,y=Z,color="13")) +
  #geom_line(data=output14,aes(x=t,y=Z,color="14")) +
  geom_line(data=output15,aes(x=t,y=Z,color="15")) +
  geom_line(data=output16,aes(x=t,y=Z,color="16")) +
  geom_line(data=output17,aes(x=t,y=Z,color="17")) +
  geom_line(data = output5, aes(x=t,y=eexp), linewidth = 1.5) +
  scale_fill_viridis() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size=27))

final_plot  
ggsave("./simulations_time/birthdeath/bd_plot.png",dpi=600)

output12 %>% 
  ggplot(aes(x=t)) + 
  geom_point(aes(y=Z, color="simulation")) +
  geom_line(aes(y=eexp, color="exponential")) +
  scale_color_manual(values=c("black","indianred"),labels = c(expression(e^{(alpha - beta)*t}),"simulation")) +
  labs(color=NULL) +
  theme(legend.text = element_text(size=18),
        axis.title = element_text(size=16),
        axis.text = element_text(size=12)) +
  theme(legend.text.align = 0)
  

rm(list = ls())
  