library(reshape2)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)

# define parameters
alpha = 20;beta = 10;
# population starting with 1 cell
Z = 1; t = 0
# create tibble to store Z values for each t
output <- tibble("t" = t,"Z" = Z)

while (t < 0.8) { 
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
saveRDS(output, file = paste0("./simulations_time/birthdeath_12.rds"))
output <- readRDS("./simulations_time/birthdeath/birthdeath_12.rds")

output %>% 
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
  