
###### data visualization ######

sub.corn <- subset(lasrosas.corn, year==2001)

ggplot(sub.corn) + geom_point(aes(long,lat,colour=yield)) + 
  scale_colour_gradient(low = "yellow",high = "green4") + 
  xlab("Longitude") + ylab("Latitude") +
  labs(color = "Yield") + thm1
ggplot(sub.corn) + geom_point(aes(long,lat,colour=topo)) + 
  xlab("Longitude") + ylab("Latitude") + 
  labs(color = "Topographic Factor") + thm1



ggplot(sub.corn) + geom_density_ridges_gradient(aes(yield,topo,fill=topo)) + 
  xlab("Yield") + ylab("Topographic Factor") + 
  thm1 + theme(legend.position = "na")


ggplot(sub.corn) + 
  geom_histogram(aes(x=yield,y=..density..), colour="black", fill="white")+
  geom_density(aes(yield,colour=yield),size=2,alpha=.2, fill="#56B4E9") + 
  xlab("Yield") + ylab("Density") + 
  xlim(c(10,130)) + thm1


ggplot(sub.corn) + geom_point(aes(long,lat,colour=rep)) + 
  xlab("Longitude") + ylab("") + 
  labs(color = "Reps") + thm1

ggplot(sub.corn) + geom_density_ridges_gradient(aes(yield,rep,fill=rep)) + 
  xlab("Yield") + ylab("Replicate") + 
  thm1 + theme(legend.position = "na")


ggplot(data=gridded_df_ordered,
       aes(as.numeric(gridId), yield, color = topo)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = lm,size=2) +
  thm1 + xlab("Index") + ylab("Yield") + 
  labs(color = "Topographic Factor")

ggplot(data=gridded_df_ordered,
       aes(as.numeric(gridId), yield)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = lm) +
  thm1

ggplot(data=gridded_df_ordered,
       aes(as.numeric(gridId), yield, color = rep)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = lm) +
  thm1


ggplot(data=gridded_df_ordered,
       aes(as.numeric(gridId), yield, color = row)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = lm) +
  theme(legend.position = "na")


ggplot(data=gridded_df_ordered,
       aes(as.numeric(gridId), yield, color = bv)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = lm) +
  thm1






########### prior ############

set.seed(1837273)
R1 <- rlkjcorr(1e4,K=2,eta=0.1)
R2 <- rlkjcorr(1e4,K=2,eta=1)
R3 <- rlkjcorr(1e4,K=2,eta=5)
R4 <- rlkjcorr(1e4,K=2,eta=15)

lkj.data <- data.frame(rhos=c(R1[,1,2],R2[,1,2],R3[,1,2],R4[,1,2]),
                       eta =c(rep(0.1,1e4),rep(1,1e4),rep(5,1e4),rep(15,1e4)))
lkj.data$eta <- factor(lkj.data$eta)

ggplot(lkj.data) + 
  stat_density(aes(rhos,y=..density..,group=eta,colour=eta),
               geom="line",position="identity", size = 1) + 
  xlab(expression(rho)) + ylab("Density") + 
  labs(color=expression(paste(epsilon," "))) + thm1


R5 <- rlkjcorr(1e4,K=3,eta=1)
plot(R5[,1,2], R5[,1,3] , col=col.alpha("black",0.2) , pch=16 )

lkj.data2 <- data.frame(r12=R5[,1,2],r13=R5[,1,3])
ggplot(lkj.data2) + geom_point(aes(r12,r13),alpha = 0.3) + 
  xlab(expression(rho[12])) + ylab(expression(rho[13])) + 
  thm1


############ prior checking ############

## t distribution ## scale parameter 
s  <- seq(-5,5,length=1e4)
sigs<- c(1,2,5)
ds <- matrix(0,length(sigs),length(s))
for(i in 1:length(sigs))
  ds[i,] <- dstudent_t(s,3,0,sigs[i])

t.data <- data.frame(s = rep(s,3),dens=c(ds[1,],ds[2,],ds[3,]),
                     sigs = rep(sigs,each=length(s)))
t.data$sigs <- as.factor(t.data$sigs)
ggplot(t.data) + 
  geom_line(aes(s,dens,group=sigs,colour=sigs),size = 1) + 
  ylab("Density") + xlab("") + 
  labs(color="scale parameter")+ thm1



set.seed(seed = 19284)

N  <- datAll$N
M <- 100

y_prior <- matrix(0,M,N)
for(i in 1:M){
  tau0  <- 1 / sqrt(rgamma(1, 1, rate = 100))
  sigma <- 1 / sqrt(rgamma(1, 1, rate = 100))
  beta0 <- rnorm(1,median(datAll$Y),100)
  beta12 <- rnorm(2, 0, 100)
  beta0i <- rnorm(3, 0, tau0)
  y_prior[i,] <- beta0 + datAll$Xc%*%beta12+
    beta0i[datAll$J_1]+ rnorm(N, mean = 0, sd = sigma)
}
pp_check(datAll$Y,y_prior,fun = "ppc_dens_overlay") + thm1


y_prior1 <- matrix(0,M,N)
for(i in 1:M){
  tau0  <- abs(rstudent(1,3, 0, 1))
  sigma <- abs(rstudent(1,3, 0, 1))
  beta0 <- rstudent(1,3,median(gridded_df_ordered$yield),30)
  beta12 <- rnorm(2,0, 0.01)
  beta0i <- rstudent(3,3, 0, tau0)

  y_prior1[i,] <- beta0 + datAll$Xc%*%beta12 + beta0i[datAll$J_1]+
    rnorm(N, mean = 0, sd = sigma)
}
pp_check(datAll$Y,y_prior1,fun = "ppc_dens_overlay") + thm1 



y_prior2 <- matrix(0,M,N)
for(i in 1:M){
  tau0  <- abs(rstudent(1,3, 0, 30))
  sigma <- abs(rstudent(1,3, 0, 30))
  beta0 <- rstudent(1,3,median(gridded_df_ordered$yield),5)
  beta1 <- rnorm(1,0, 0.1)
  beta2 <- rnorm(1,0, 0.01)
  beta0i <- rstudent(3,3, 0, tau0)
  
  y_prior2[i,] <-beta0 + datAll$Xc%*%c(beta1,beta2)+ 
    beta0i[datAll$J_1]+
    rnorm(N, mean = 0, sd = sigma)
}
pp_check(datAll$Y,y_prior2,fun = "ppc_dens_overlay") + thm1


# ppc_stat(datAll$Y,y_prior2)


##### model checking #####
pp_check(mod.null) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_modelnull.pdf")

pp_check(mod.1) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model1.pdf")



pp_check(mod.2) + thm1 + 
  coord_cartesian(xlim = c(-10, 150))
ggsave("Plots/ppcheck_model2.pdf")



pp_check(mod.4) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model3.pdf")


pp_check(mod.4) + thm1 + 
  coord_cartesian(xlim = c(-10, 150))
ggsave("Plots/ppcheck_model4.pdf")




pp_check(mod.5) + thm1 + ## with priors
  coord_cartesian(xlim = c(-10, 160), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model5.pdf")






##### model checking #####
pp_check(mod.brm.1) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model1.pdf")

pp_check(mod.brm.2) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model2.pdf")


pp_check(mod.brm.3) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model3.pdf")

pp_check(mod.brm.4) + thm1 + 
  coord_cartesian(xlim = c(-10, 150), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model4.pdf")

pp_check(mod.brm.5) + thm1 + ## with priors
  coord_cartesian(xlim = c(-10, 160), ylim = c(0, 0.025))
ggsave("Plots/ppcheck_model5.pdf")
