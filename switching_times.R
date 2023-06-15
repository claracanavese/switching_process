library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(patchwork)

df1_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.288,0.702,0.334,0.641,0.441,0.702))
df1_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.33,0.69,0.38,0.66,0.49,0.69), errors = c(0.11,0.08,0.12,0.09,0.11,0.08))

plot1 <- ggplot() +
  geom_jitter(data=df1_theory, aes(x=case,y=times,color="thoery"),size=4, shape = 17, height = 0, width = 0.03) + 
  geom_point(data=df1_sim, aes(x=case,y=times,color="simulation"),size=3) +
  geom_errorbar(data=df1_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_color_manual(values=c("black","chartreuse3")) +
  labs(col=NULL) +
  guides(colour = guide_legend(override.aes = list(shape = c(16,17)))) +
  theme(legend.position = "none")
plot1

df2_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.265,0.500,0.300,0.459,0.380,0.500))
df2_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.27,0.51,0.30,0.50,0.38,0.53), errors = c(0.09,0.06,0.09,0.07,0.08,0.06))

plot2 <- ggplot() +
  geom_jitter(data=df2_theory, aes(x=case,y=times),size=4, shape = 17,color="chartreuse3",height = 0, width = 0.1) + 
  geom_point(data=df2_sim, aes(x=case,y=times),size=3) +
  geom_errorbar(data=df2_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) 
plot2

df3_theory <- data.frame(case = c(1,1,2,2,3,3),times = c(0.462,0.762,0.530,0.750,0.691,0.876))
df3_sim <- data.frame(case = c(1,1,2,2,3,3), times = c(0.46,0.74,0.53,0.72,0.68,0.83), errors = c(0.14,0.11,0.15,0.12,0.16,0.13))

plot3 <- ggplot() +
  geom_jitter(data=df3_theory, aes(x=case,y=times),size=4, shape = 17,color="chartreuse3",height = 0, width = 0.1) + 
  geom_point(data=df3_sim, aes(x=case,y=times),size=3) +
  geom_errorbar(data=df3_sim, aes(x=case, ymin=times-errors, ymax=times+errors),width=.08,color="black") +
  scale_x_continuous(breaks = c(1,2,3)) 
plot3
plot1 + plot2 + plot3

ggsave("./imgs/switching_times/jitter.png",dpi=600)

# cumulative distribution plot
omega_m = 0.01
omega_p = 0.1
s <- 1; v <- 0.5
x <- seq(-1.1, 1.1, length = 201)

x_list = lapply(seq(0, 1.8, length = 100), function(t){
  return(-omega_m*2*exp(20*t)/20)
})
t_list = seq(0, 1.8, length = 100)
x_num <- as.numeric(unlist(x_list))

plot(t_list, lerch(x_num, s = s, v = v), type = "l", col = "blue",
     las = 1, main = paste0("lerch(x, s = ", s,", v = ", v, ")"))
# abline(v = 0, h = 1, lty = "dashed", col = "gray")

C2 <- function(t) {
  lerchphi = lerch(-omega_m*2*exp(20*t)/20, s = 1, v = 0.5)
  return(omega_p*omega_m/400*2*exp(20*t)*lerchphi)
}
C2(0.1)
C2(0.5)
plot(t_list, time_cumulative3(t_list), type = "l", col = "blue",
     las = 1, main = paste0("P(T2)>t"), xlab="t", ylab="C2")

lerch(-omega_m*2*exp(20*0.5)/20, s = 1, v = 0.5, tolerance = 1.0e-10, iter = 1000)
