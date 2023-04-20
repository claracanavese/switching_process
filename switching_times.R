library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(patchwork)

df1_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.288,0.702,0.334,0.641,0.441,0.702))
df1_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.33,0.69,0.38,0.66,0.49,0.69), errors = c(0.11,0.08,0.12,0.09,0.11,0.08))

plot1 <- ggplot() +
  geom_point(data=df1_theory, aes(x=case,y=times,color="thoery"),size=4, shape = 17) + 
  geom_point(data=df1_sim, aes(x=case,y=times,color="simulation"),size=3) +
  geom_errorbar(data=df1_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_color_manual(values=c("black","chartreuse3")) +
  labs(col=NULL) +
  guides(colour = guide_legend(override.aes = list(shape = c(16,17))))
plot1

df2_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.265,0.500,0.300,0.459,0.380,0.500))
df2_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.27,0.51,0.30,0.50,0.38,0.53), errors = c(0.09,0.06,0.09,0.07,0.08,0.06))

plot2 <- ggplot() +
  geom_point(data=df2_theory, aes(x=case,y=times),size=4, shape = 17,color="chartreuse3") + 
  geom_point(data=df2_sim, aes(x=case,y=times),size=3) +
  geom_errorbar(data=df2_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) 
plot2

df3_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.462,0.762,0.530,0.750,0.691,0.876))
df3_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.46,0.74,0.53,0.72,0.68,0.83), errors = c(0.14,0.11,0.15,0.12,0.16,0.13))

plot3 <- ggplot() +
  geom_point(data=df3_theory, aes(x=case,y=times),size=4, shape = 17,color="chartreuse3") + 
  geom_point(data=df3_sim, aes(x=case,y=times),size=3) +
  geom_errorbar(data=df3_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) 
plot3
plot1 + plot2 + plot3
