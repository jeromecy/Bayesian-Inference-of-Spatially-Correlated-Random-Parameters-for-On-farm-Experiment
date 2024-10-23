source("Initialization.R")
summary(gridded_df_ordered)

####### for paper: stan and brms models #########

######## model 1 #############
### replace cnitro, cnitro.sq with nitro, nitro.sq
form.1 <- brmsformula(yield ~ 1 + nitro + nitro.sq + (1|rep))

prior.1 <- get_prior(form.1,data=gridded_df_ordered)

prior.1$prior[2] <- "normal(0,0.1)"
prior.1$prior[3] <- "normal(0,0.01)"
prior.1$prior[4]   <- "student_t(3, 84.7, 5)"
prior.1$prior[5:8] <- "student_t(3,0, 10)"
prior.1$bound[5:8] <- 0

mod.1 <- brm(form.1,data=gridded_df_ordered,
             prior = prior.1,
             sample_prior = "only",
             seed = 10403,
             cores = 4)

pp_check(mod.1) + thm1
summary(mod.1)


# make_stancode(form.1,data=gridded_df_ordered)


mod.1.update <- update(mod.1, sample_prior="no",
                       seed = 202,
                       cores=4)

pp_check(mod.1.update) + thm1
summary(mod.1.update)


dat.1 <- as.data.frame(extract(mod.1$fit))
dat.1.update <- as.data.frame(extract(mod.1.update$fit))

mod1.prior <- data.frame(mean=apply(dat.1, 2,mean),sd=apply(dat.1, 2,sd),
                    par = names(dat.1),pt="prior")
mod1.post <- data.frame(mean=apply(dat.1.update, 2,mean),sd=apply(dat.1.update, 2,sd),
                    par = names(dat.1.update),pt="posterior")

mod1.comb <- rbind(mod1.prior[c(2:8),],mod1.post[c(2:8),])

ggplot(mod1.comb, aes(x=par, y=mean, colour=par)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  facet_wrap(~pt) + 
  geom_point(size=3) +
  scale_x_discrete(labels=c(expression(beta[1]),expression(beta[2]),
                            expression(b[1]),expression(b[2]),
                            expression(b[3]),expression(sigma[b]),
                            expression(sigma[e]))) +
  xlab("") + ylab("") + coord_flip() + thm1 +
  theme(legend.position = "na")

model1.Inter <- data.frame(Inter = c(dat.1$b_Intercept,dat.1.update$b_Intercept),
                           pt = c(rep("prior",4e3),rep("post",4e3)))

ggplot(model1.Inter) + geom_density_ridges(aes(Inter,pt,fill=pt),alpha=0.5) + 
  coord_cartesian(xlim=c(-50,250),ylim=c(1,2.5)) + xlab("Intercept") + ylab("")+
  thm1 + theme(legend.position = "na")




######## model 2 #############
form.2 <- brmsformula(yield ~ 1 + nitro + nitro.sq + (1|rep) +
                        (1|row+col))

prior.2 <- get_prior(form.2,data=gridded_df_ordered)

prior.2$prior[2]  <- "normal(0,0.1)"
prior.2$prior[3]  <- "normal(0,0.01)"
prior.2$prior[4]    <- "student_t(3, 84.7, 5)"
prior.2$prior[5:11] <- "student_t(3,0, 30)"
prior.2$prior[12] <- "student_t(3,0, 30)"
prior.2$bound[5:12] <- 0

mod.2 <- brm(form.2,data=gridded_df_ordered,
             prior = prior.2,
             sample_prior = "only",
             seed =  515190382,
             cores = 4)

pp_check(mod.2) + thm1
summary(mod.2)


# make_stancode(form.2,data=gridded_df_ordered,prior = prior.2)
data.2 <- make_standata(form.2,data=gridded_df_ordered,prior = prior.2)

mod.2.update <- update(mod.2, sample_prior="no",
                       seed = 2133433928,
                       cores=4)

pp_check(mod.2.update) + thm1
summary(mod.2.update)



dat.2 <- as.data.frame(extract(mod.2$fit))
dat.2.update <- as.data.frame(extract(mod.2.update$fit))

mod2.prior <- data.frame(mean=apply(dat.2, 2,mean),sd=apply(dat.2, 2,sd),
                         par = names(dat.2),pt="prior")
mod2.post <- data.frame(mean=apply(dat.2.update, 2,mean),sd=apply(dat.2.update, 2,sd),
                        par = names(dat.2.update),pt="posterior")

mod2.comb <- rbind(mod2.prior[c(2:7,26:28),],mod2.post[c(2:7,26:28),])

ggplot(mod2.comb, aes(x=par, y=mean, colour=par)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  facet_wrap(~pt) + 
  geom_point(size=3) +
  scale_x_discrete(labels=c(expression(beta[1]),expression(beta[2]),
                            expression(b[1]),expression(b[2]),
                            expression(b[3]),expression(sigma[c]),
                            expression(sigma[b]),expression(sigma[r]),
                            expression(sigma[e]))) +
  xlab("") + ylab("") + coord_flip() + thm1 +
  theme(legend.position = "na")

model2.Inter <- data.frame(Inter = c(dat.2$b_Intercept,dat.2.update$b_Intercept),
                           pt = c(rep("prior",4e3),rep("post",4e3)))

ggplot(model2.Inter) + geom_density_ridges(aes(Inter,pt,fill=pt),alpha=0.5) + 
  coord_cartesian(xlim=c(-50,250),ylim=c(1,2.5)) + xlab("Intercept") + ylab("")+
  thm1 + theme(legend.position = "na")



form.21 <- brmsformula(yield ~ 1 + nitro + nitro.sq + (1|topo) +
                        (1|row+col))


mod.21 <- brm(form.21,data=gridded_df_ordered,
              seed =  515190382,
              cores = 4)

pp_check(mod.21) + thm1
summary(mod.21)

loo(mod.21)
loo(mod.2.update)



######## model 3 AR1xAR1 #############


form.3 <- brmsformula(yield ~ 1 + nitro + nitro.sq + (1|rep))

data.3 <- make_standata(form.3,data=gridded_df_ordered)
data.3$nrow <- 93
data.3$ncol <- 18
# saveRDS(data.3,"Data/data3.rds")

# stan.3 <- stan("StanModel/spatialAR1.stan",
#                data=data.3,
#                chains = 1,iter = 10)

# stan.3 <- readRDS("Outcome/StanSpatialAR1_PC.rds")
stan.3 <- readRDS("Outcome/StanSpatialAR1prior_Pawsey.rds")

print(stan.3)

# parNms <- c("b[1]","b[2]","sigma","sigma_ar","sd_1[1]","rho_r","rho_c")

dat.3.update <- as.data.frame(extract(stan.3))[,-c(3,7:9,15)]

mod3.post <- data.frame(mean=apply(dat.3.update, 2,mean),
                        sd=apply(dat.3.update, 2,sd),
                        par = names(dat.3.update),pt="posterior")

ggplot(mod3.post, aes(x=par, y=mean, colour=par)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=3) +
  scale_x_discrete(labels=c(expression(beta[1]),expression(beta[2]),
                            expression(b[1]),expression(b[2]),
                            expression(b[3]),expression(rho[c]),
                            expression(rho[r]),expression(sigma[b]),
                            expression(sigma[e]),expression(sigma[xi]))) +
  xlab("") + ylab("") + coord_flip() + thm1 +
  theme(legend.position = "na")

# ggsave("Plots/Model3PostPoints.pdf",width=11.96,height=8.27)


# dat.3.melt <- melt(dat.3.update)
# dat.3.melt$variable<- factor(dat.3.melt$variable,levels = levels(mod3.post$par))
# ggplot(dat.3.melt) + 
#   geom_density_ridges(aes(value,variable,fill=variable)) +
#   coord_cartesian(xlim=c(-5,20))+
#   scale_y_discrete(labels=c(expression(beta[1]),expression(beta[2]),
#                             expression(b[1]),expression(b[2]),
#                             expression(b[3]),expression(rho[c]),
#                             expression(rho[r]),expression(sigma[b]),
#                             expression(sigma[e]),expression(sigma[xi])))+
#   xlab("") + ylab("") + thm1 +
#   theme(legend.position = "na")

var_names <- c(
  "b.1" = "beta[1]",
  "b.2" = "beta[2]",
  "r_1_1.1" = "b[1]",
  "r_1_1.2" = "b[2]",
  "r_1_1.3" = "b[3]",
  "rho_c"   = "rho[c]",
  "rho_r" = "rho[r]",
  "sd_1" = "sigma[b]",
  "sigma" = "sigma[e]",
  "sigma_ar" = "sigma[ar]"
)


ggplot(dat.3.melt) + 
  stat_density(aes(value,fill=variable)) + 
  facet_wrap(~variable,scales="free",
             labeller = labeller(variable = as_labeller(var_names,label_parsed)))+
  xlab("") + ylab("") + thm1 +
  theme(legend.position = "na",axis.text = element_text(size=10))

ggsave("Plots/Model3PostDensity02.pdf",width=11.96,height=8.27)






dat.3.pp <- as.array(rstan::extract(stan.3))





N <- data.3$N
i <- 1
j <- 1
M <- 100 # length(dat.3.pp$b)

Xc<- matrix(0,N, 2)
means_X <- numeric(2)
for (i in 2:3) {
  means_X[i - 1] = mean(data.3$X[, i]);
  Xc[, i - 1] = data.3$X[, i] - means_X[i - 1];
}

y_rep_ar1  <- matrix(0,M,N)
# loglik_ar1 <- matrix(0,M,N)

mus.ar1 <- matrix(0,M,N)
set.seed(128374)

for(i in 1:M){
  mus.ar1[i,] = dat.3.pp$Intercept[i] + Xc %*% dat.3.pp$b[i,] +
    dat.3.pp$r_1_1[i,data.3$J_1]
  
  BigAR1 = dat.3.pp$sigma_ar[i]^2*kronecker(AR_matrix(dat.3.pp$rho_r[i],data.3$nrow),
                                            AR_matrix(dat.3.pp$rho_c[i],data.3$ncol))
  BigSig = BigAR1+ diag(dat.3.pp$sigma[i]^2,data.3$N)
  chol_BigSig <- chol(BigSig)
  y_rep_ar1[i,] <- mus.ar1[i,] + chol_BigSig%*%rnorm(nrow(chol_BigSig))
  
  # loglik_ar1[i,] <- dmvnorm(x=as.vector(data.3$Y),mean=mus.ar1[i,],
  #                           sigma=BigSig,
  #                           log=TRUE)
  cat(paste0(round(i/M*100,2),"% "))
}
pp_check(data.3$Y,y_rep_ar1,ppc_dens_overlay) + thm1

# ggsave("Plots/ppcheck_model3.pdf",width=11.96,height=8.27)





######## model 4 TPRS #############


form.4 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1|rep) + s(x,y))

prior.4 <- get_prior(form.4,data=gridded_df_ordered)

prior.4$prior[1:5] <- "normal(0,0.1)"
prior.4$prior[3] <- "normal(0,0.01)"
prior.4$prior[6] <- "student_t(3, 84.7, 5)"
prior.4$prior[c(7,10,12)] <- "student_t(3,0, 30)"
prior.4$bound[7:12] <- 0

mod.4 <- brm(form.4,data=gridded_df_ordered,
             prior = prior.4,
             sample_prior = "only",
             seed = 1283494313,
             cores = 4)

pp_check(mod.4) + thm1
summary(mod.4)


mod.4.update <- update(mod.4, sample_prior="no",
                       seed = 1101294840,
                       cores=4)

pp_check(mod.4.update) + thm1
ggsave("Plots/ppcheck_model4.pdf",width=11.96,height=8.27)

summary(mod.4.update)



dat.4 <- as.data.frame(rstan::extract(mod.4$fit))
dat.4.update <- as.data.frame(rstan::extract(mod.4.update$fit))

mod4.prior <- data.frame(mean=apply(dat.4, 2,mean),sd=apply(dat.4, 2,sd),
                         par = names(dat.4),pt="prior")
mod4.post <- data.frame(mean=apply(dat.4.update, 2,mean),sd=apply(dat.4.update, 2,sd),
                        par = names(dat.4.update),pt="posterior")

mod4.comb <- rbind(mod4.prior[c(2:11),],mod4.post[c(2:11),])

new_order <- with(mod4.comb, reorder(par, mean, median))
mod4.comb$par <- factor(mod4.comb$par,levels = levels(new_order))
ggplot(mod4.comb, aes(x=par, y=mean, colour=par)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  facet_wrap(~pt) + 
  geom_point(size=3) +
  # scale_x_discrete(labels=c(expression(b[2]),expression(b[3]),
  #                           expression(b[2]),expression(beta[2]),
  #                           expression(beta[2]),
  #                           expression(sigma[b]),expression(sigma))) +
  xlab("") + ylab("") + coord_flip() + thm1 +
  theme(legend.position = "na")

model2.Inter <- data.frame(Inter = c(dat.4$b_Intercept,dat.4.update$b_Intercept),
                           pt = c(rep("prior",4e3),rep("post",4e3)))

ggplot(model2.Inter) + geom_density_ridges(aes(Inter,pt,fill=pt),alpha=0.5) + 
  coord_cartesian(xlim=c(-50,250),ylim=c(1.5,2.5)) + xlab("Intercept") + ylab("")+
  thm1 + theme(legend.position = "na")





form.4.2 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                          (1|row+col)+ (1|rep) + s(x,y))

prior.4.2 <- get_prior(form.4.2,data=gridded_df_ordered)

prior.4.2$prior[1:5] <- "normal(0,0.1)"
prior.4.2$prior[3] <- "normal(0,0.01)"
prior.4.2$prior[6] <- "student_t(3, 84.7, 5)"
prior.4.2$prior[c(7,14,16)] <- "student_t(3,0, 30)"
prior.4.2$bound[7:16] <- 0

mod.4.2 <- brm(form.4.2,data=gridded_df_ordered,
              prior = prior.4.2,
              sample_prior = "only",
              seed = 1682577970,
              cores = 4)


mod.4.2.update <- update(mod.4.2, sample_prior="no",
                         seed = 905136004,
                         cores=4)

pp_check(mod.4.2.update) + thm1

bayes_R2(mod.4.2.update)

bayes_R2(mod.2.update)
bayes_R2(mod.4.update)





########## random U ~ ar1xar1 ###############


## seeds  1072473474
stanRandU.rep <- readRDS("Outcome/StanRandUAR1_PC.rds")

check_hmc_diagnostics(stanRandU.rep)
get_num_leapfrog_per_iteration(stanRandU.rep)

check_energy(stanRandU.rep)

names(stanRandU.rep)

post.rand.rep <- as.array(rstan::extract(stanRandU.rep))
names(post.rand.rep)

dim(post.rand.rep$r_1_1)

traceplot(stanRandU.rep,c("rho_r","rho_c"),chain=1)


form.rand.rep <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep)+
                           (1 + cnitro + cnitro.sq|gridId))

datarand.rep <- make_standata(form.rand.rep,data=gridded_df_ordered)
datarand.rep$ncol <- 18
datarand.rep$nrow <- 93

N <- datarand.rep$N
i <- 1
j <- 1
M <- length(post.rand.rep$b_Intercept)

y_rep_randU.rep <- matrix(0,M,N)
loglik_randU.rep <- matrix(0,M,N)
mus.rep <- matrix(0,M,N)
set.seed(4885514)
for(i in 1:M){

  mus.rep[i,] <- post.rand.rep$b_Intercept[i] + 
    datarand.rep$X[,2:3]%*%post.rand.rep$b[i,] +
    as.numeric(post.rand.rep$r_1_1[i,]*datarand.rep$Z_1_1) + 
    as.numeric(post.rand.rep$r_1_2[i,]*datarand.rep$Z_1_2) +
    as.numeric(post.rand.rep$r_1_3[i,]*datarand.rep$Z_1_3) +
    post.rand.rep$r_2_1[i,datarand.rep$J_2] 
  
  y_rep_randU.rep[i,] <- rnorm(N,mus.rep[i,],post.rand.rep$sigma[i])
  loglik_randU.rep[i,] <- dnorm(datarand.rep$Y,
                            mus.rep[i,],
                            abs(post.rand.rep$sigma[i]),
                            log = TRUE)
  cat(paste0(round(i/M*100,2),"% "))
}

pp_check(datarand.rep$Y, y_rep_randU.rep[sample(1:M,100),], ppc_dens_overlay) + thm1
# ggsave("Plots/ppcheck_ModelRandU.pdf",width=11.96,height=8.27)
ppc_stat(datarand.rep$Y, y_rep_randU.rep[sample(1:M,100),])


var_mu.rep <- apply(mus.rep, 1, var)
sigma2.rep <- post.rand.rep$sigma
rsq_bayes.rep <- var_mu.rep / (var_mu.rep + sigma2.rep)

summary(rsq_bayes.rep)


psis_result.rep <- psis(-loglik_randU.rep)
plot(psis_result.rep, label_points = TRUE)

loo(loglik_randU.rep)

wh <- which(psis_result.rep$diagnostics$pareto_k>0.7)

plot(datarand.rep$Y)
points(wh,datarand.rep$Y[wh],col="red",pch=20)


loo_moment_match(stanRandU)



hist(post.rand.rep$Cor_1[,1,2],30)

post.rand.rep$cor_1

prdX <- matrix(c(10.5,-445.266666666667),N,2,byrow = TRUE)

i <- 1
j <- 1
M <- 50

y_rep_randU_prd <- matrix(0,M,N)

ct <- 1
for(i in sample(1:N,M)){
  mu <- post.rand.rep$b_Intercept[i] + prdX%*%post.rand.rep$b[i,]
  for(j in 1:N){
    mu[j] <- mu[j] + post.rand.rep$r_1_1[i,datarand$J_1[j]] + 
      post.rand.rep$r_1_2[i,datarand$J_1[j]] + 
      post.rand.rep$r_1_3[i,datarand$J_1[j]] + 
      post.rand.rep$r_2_1[i,datarand$J_2[j]] 
  }
  y_rep_randU_prd[ct,] <- rnorm(N,mu,post.rand.rep$sigma[i])
  cat(ct)
  ct <- ct + 1
}

prd <- apply(y_rep_randU_prd,2,mean)

prd.dat <- data.frame(gridId = gridded_df_ordered$gridId)

prd.dat$x <- as.data.frame(gridded_df_ordered$x)
prd.dat$y <- as.data.frame(gridded_df_ordered$y)
prd.dat$yield <- prd

sp::coordinates(prd.dat) <- c("x","y")

spplot(prd.dat, 'yield')



mod.yield <- as(prd.dat["yield"], "ppp")
yield.smooth <- Smooth(mod.yield)
plot(yield.smooth, col=rev(rainbow(255)),  ribsep=0.01, 
     ribside="bottom", main="")
options(scipen = 999)

plot(yield.smooth, col=matlab.like2,  ribsep=0.001, 
     ribside="bottom", main="")
contour(yield.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
## ggsave("Plots/ModelRandU_prd.pdf",width=11.96,height=8.27)

i <- 1
j <- 1
M <- 50
ct <- 1

beta0 <- matrix(0,M,N)
beta1 <- matrix(0,M,N)
beta2 <- matrix(0,M,N)
for(i in sample(1:N,M)){
  beta0[ct,] <- post.rand$b_Intercept[i] + post.rand$r_1_1[i,] - 
    post.rand$r_1_2[i,]*means_X[1]-post.rand$r_1_3[i,]*means_X[2]
  beta1[ct,] <- post.rand$b[i,1] + post.rand$r_1_2[i,]
  beta2[ct,] <- post.rand$b[i,2] + post.rand$r_1_3[i,]
  cat(ct)
  ct <- ct + 1
}

prd.dat$b0 <- apply(beta0,2,mean)
prd.dat$b1 <- apply(beta1,2,mean)
prd.dat$b2 <- apply(beta2,2,mean)



mod.b0 <- as(prd.dat["b0"], "ppp")
b0.smooth <- Smooth(mod.b0)
plot(b0.smooth, col=matlab.like2,  ribsep=0.001, 
     ribside="bottom", main="")
contour(b0.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)

mod.b1 <- as(prd.dat["b1"], "ppp")
b1.smooth <- Smooth(mod.b1)
plot(b1.smooth, col=matlab.like2,  ribsep=0.001, 
     ribside="bottom", main="")
contour(b1.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)


mod.b2 <- as(prd.dat["b2"], "ppp")
b2.smooth <- Smooth(mod.b2)
plot(b2.smooth, col=matlab.like2,  ribsep=0.001, 
     ribside="bottom", main="")
contour(b2.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)


post.rand$b_Intercept[i]

r1 <- apply(post.rand$r_1_1,2,mean) + mean(post.rand$b_Intercept)
r2 <- apply(post.rand$r_1_2,2,mean) + 0.0502569140 
r3 <- apply(post.rand$r_1_3,2,mean) - 0.0001612903





########## FULL model long time ###############


## seeds  1072473474
stanFULL <- readRDS("Outcome/StanFullModel_PC.rds")

check_hmc_diagnostics(stanFULL)

post.full <- as.array(rstan::extract(stanFULL))
names(post.full)


traceplot(stanFULL,c("rho_r","rho_c"))


form.full <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep)+
                           (1 + cnitro + cnitro.sq|gridId)+
                           s(x,y),
                         sigma~(1|gridId))

datafull <- make_standata(form.full,data=gridded_df_ordered)
datafull$ncol <- 18
datafull$nrow <- 93



N <- datafull$N
i <- 1
j <- 1
M <- length(post.full$Intercept)

y_rep_full <- matrix(0,M,N)
loglik_full<- matrix(0,M,N)
mus.full   <- matrix(0,M,N)
sigma.full <- matrix(0,M,N)
set.seed(18473012)
for(i in 1:M){
  mus.full[i,] <- post.full$b_Intercept[i] + 
    datafull$X[,2:3]%*%post.full$b[i,] + 
    datafull$Xs %*% post.full$bs[i,] + 
    datafull$Zs_1_1 %*% post.full$s_1_1[i,] +
    as.numeric(post.full$r_1_1[i,]*datafull$Z_1_1) + 
    as.numeric(post.full$r_1_2[i,]*datafull$Z_1_2) + 
    as.numeric(post.full$r_1_3[i,]*datafull$Z_1_3) + 
    as.numeric(post.full$r_2_1[i,datafull$J_2])
  
  sigma.full[i,] <- post.full$Intercept_sigma[i]+
    post.full$r_3_sigma_1[i,]
  
  y_rep_full[i,] <- rnorm(N,mus.full[i,],exp(sigma.full[i,]))
  loglik_full[i,] <- dnorm(datafull$Y,
                           mus.full[i,],
                           exp(sigma.full[i,]),log = TRUE)
  
  cat(paste0(round(i/M*100,2),"% "))
}

pp_check(datafull$Y, y_rep_full[sample(1:M,100),], ppc_dens_overlay) + thm1

# ggsave("Plots/ppcheck_ModelLongFull.pdf",width=11.96,height=8.27)


var_mu.full <- apply(mus.full, 1, var)
sigma2.full <- apply(sigma.full, 1, var)
rsq_bayes.full <- var_mu.full / (var_mu.full + sigma2.full)

summary(rsq_bayes.full)
quantile(rsq_bayes.full,c(0.025,0.5,0.975))

loo(loglik_full)

psis_result.full <- psis(-loglik_full)
plot(psis_result.full, label_points = TRUE)

pareto_k_table(psis_result.full)


M <- nrow(post.full$b)
loglik.full <- yloo.full <- sdloo.full <- matrix(nrow = M, ncol = N)
for (i in 1:M) {
  g.full <- (as.numeric(datafull$Y) - mus.full[i,])/exp(sigma.full[i,])^2
  yloo.full[i, ]   <- as.numeric(datafull$Y) - g.full*exp(sigma.full[i,])^2
  sdloo.full[i, ]  <- exp(sigma.full[i,])
  loglik.full[i, ] <- dnorm(as.numeric(datafull$Y), 
                            yloo.full[i, ], 
                            sdloo.full[i, ], log = TRUE)
  cat(paste0(round(i/nrow(post.full$b)*100,2),"% "));
}

psis_result.full <- psis(-loglik.full)
plot(psis_result.full, label_points = TRUE)


loo1.full  <- loo(loglik.full)
psis1.full <- psis(-loglik.full)

pareto_k_table(psis1.full)


hmc_diagnostics_full <- nuts_params(stanFULL)
lp_full <- log_posterior(stanFULL)


mcmc_nuts_divergence(hmc_diagnostics_full,
                     log_posterior(stanFULL),
                     chain = 4)

mcmc_nuts_energy(hmc_diagnostics_full)


wh<- which(hmc_diagnostics_full$Parameter=="divergent__")
table(hmc_diagnostics_full[wh,]$Value)




i <- 1
j <- 1
M <- 50

y_rep_full_prd <- matrix(0,M,N)

ct <- 1
for(i in sample(1:N,M)){
  mu <- post.full$Intercept[i] + prdX%*%post.full$b[i,] + 
    datafull$Xs %*% post.full$bs[i,] + datafull$Zs_1_1 %*% post.full$zs_1_1[i,]
  
  for(j in 1:N){
    mu[j] <- mu[j] + post.full$r_1_1[i,datafull$J_1[j]] + 
      post.full$r_1_2[i,datafull$J_1[j]] + 
      post.full$r_1_3[i,datafull$J_1[j]] + 
      post.full$r_2_1[i,datafull$J_2[j]] 
  }
  
  sigma = post.full$Intercept_sigma[i] + post.full$r_3_sigma_1[i,]
  
  y_rep_full_prd[ct,] <- mu + sigma*rnorm(N)
  cat(ct)
  ct <- ct + 1
}

pp_check(datarand$Y,y_rep_full_prd,ppc_dens_overlay)

prd <- apply(y_rep_full_prd,2,mean)

prd.dat <- data.frame(gridId = gridded_df_ordered$gridId)

prd.dat$x <- as.data.frame(gridded_df_ordered$x)
prd.dat$y <- as.data.frame(gridded_df_ordered$y)
prd.dat$yield <- prd

sp::coordinates(prd.dat) <- c("x","y")

spplot(prd.dat, 'yield')



mod.yield <- as(prd.dat["yield"], "ppp")
yield.smooth <- Smooth(mod.yield)
plot(yield.smooth, col=rev(rainbow(255)),  ribsep=0.01, 
     ribside="bottom", main="")
options(scipen = 999)

plot(yield.smooth, col=matlab.like2,  ribsep=0.001, 
     ribside="bottom", main="")
contour(yield.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
ggsave("Plots/ModelFULL_prd.pdf",width=11.96,height=8.27)












