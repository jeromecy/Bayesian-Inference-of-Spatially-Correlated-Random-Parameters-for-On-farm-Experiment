
M <- nrow(gridded_df_ordered)

set.seed(2020)
V_u   <- rlkjcorr(n=1,K=3,eta=4)
# rho_c <- runif(1,0,1)
# rho_r <- runif(1,0,1)
b0 <- rnorm(M,75,10)
b1 <- rnorm(M,0,1)
b2 <- rnorm(M,0,0.01)

# sig_u <- c(abs(rnorm(1,11,1)),abs(rnorm(2, 0, 1)))
sig_u <- abs(rnorm(3, 0, 1))

# matC <- chol_AR_matrix(rho_c,18)
# matR <- chol_AR_matrix(rho_r,93)
# matAR1 <- chol_kronecker_product(matR,matC)
# Sigma_u <- chol_kronecker_product(t(diag(sig_u)%*%chol(V_u)),matAR1)

Sigma_u <- chol_kronecker_product(t(diag(sig_u)%*%chol(V_u)),diag(M))
us   <- Sigma_u%*%rnorm(M*3,0,1)
matU <- matrix(us,nrow=M,ncol=3,byrow = TRUE)

rand2 <- matU[,1]+matU[,2]*(gridded_df_ordered$nitro-mean(gridded_df_ordered$nitro)) + 
  matU[,3]*(gridded_df_ordered$nitro.sq-mean(gridded_df_ordered$nitro.sq))

rand1 <- b0+b1*(gridded_df_ordered$nitro-mean(gridded_df_ordered$nitro)) +
  b2*(gridded_df_ordered$nitro.sq-mean(gridded_df_ordered$nitro.sq))


sig_e <- rnorm(M, 0, 1)

sim <- rand1+rand2+sig_e

plot(sim)
plot(gridded_df_ordered$yield,sim)


plot(gridded_df_ordered$yield,rand1)


plot(gridded_df_ordered$yield,log(abs(sim)))

fit.gam <- gam(yield ~ s(x,y,k=30),data=gridded_df_ordered)

plot(gridded_df_ordered$yield-fit.gam$fitted.values,sim)

plot(gridded_df_ordered$yield,sim+fit.gam$fitted.values)


plot(gridded_df_ordered$yield-fit.gam$fitted.values,rand1)


ggplot(gridded_df_ordered) + 
  geom_histogram(aes(yield,fill=topo),bins=30)

ggplot(gridded_df_ordered) + 
  geom_density_ridges(aes(yield,y=topo))


subdata <- subset(gridded_df_ordered,topo=="HT")
M <- nrow(subdata)

set.seed(2020)
V_u   <- rlkjcorr(n=1,K=3,eta=4)
# rho_c <- runif(1,0,1)
# rho_r <- runif(1,0,1)
b0 <- rnorm(M,0,1)
b1 <- rnorm(M,0,1)
b2 <- rnorm(M,0,0.01)

sig_u <- c(abs(rnorm(1, 0, 1)),
           abs(rnorm(1, 0, 1)),
           abs(rnorm(1, 0, 0.01)))

Sigma_u <- chol_kronecker_product(t(diag(sig_u)%*%chol(V_u)),diag(M))
us   <- Sigma_u%*%rnorm(M*3,0,1)
matU <- matrix(us,nrow=M,ncol=3,byrow = TRUE)

rand2 <- matU[,1]+matU[,2]*(subdata$nitro-mean(subdata$nitro)) + 
  matU[,3]*(subdata$nitro.sq-mean(subdata$nitro.sq))

rand1 <- b0+b1*(subdata$nitro-mean(subdata$nitro)) +
  b2*(subdata$nitro.sq-mean(subdata$nitro.sq))

sig_e <- rnorm(M, 0, 1)

sim <- rand1+rand2+sig_e

plot(sim)
plot(subdata$yield,sim)

plot(subdata$yield,log(abs(sim)))

plot(subdata$yield,rand1)
plot(subdata$yield,rand2)







keep_pars <- c("rho_c","rho_r")

posterior2 <- as.array(stan.test, pars = keep_pars)

color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.15, div_alpha = 0.4)

mcmc_parcoord(
  posterior2,
  size = 0.15,
  alpha = 0.2,
  np_style = div_style
) 
mcmc_scatter(
  posterior2,
  size = 3,
  alpha = 0.2,
  np_style = div_style
)



plot(gridded_df_ordered$yield)



# Intercept    75.92      0.91    74.15    77.70 1.00     3330     2813
# cnitro        0.05      0.02     0.02     0.08 1.00     2470     1955
# cnitro.sq    -0.00      0.00    -0.00     0.00 1.00     4783     2643



# Plot: prior predictive with weakly informative priors
set.seed(seed = 1923840479)
tau0 <- abs(rnorm(1, 0, 1))
tau1 <- abs(rnorm(1, 0, 1))
sigma <- abs(rnorm(1, 0, 1))
beta0i <- rnorm(8, 0, tau0)
beta1i <- rnorm(8, 0, tau1)
beta0 <- rnorm(1, 0, 1)
beta1 <- rnorm(1, 1, 1)

data2 <- data.frame(
  log_pm25 = GM$log_pm25,
  sim = beta0 + beta0i[GM$super_region] +
    (beta1 + beta1i[GM$super_region]) * GM$log_sat +
    rnorm(Nsim, mean = 0, sd = sigma)
)





prior.1$prior[1:3] <- "normal(0,1)"
prior.brm.1 <- brm(bf.1,data=gridded_df_ordered,
                   cores = 4,
                   prior = prior.1,
                   sample_prior = "only")

pp_check(prior.brm.1)+thm1


prior.2 <- prior.1
prior.2$prior[1:3] <- "cauchy(0,1)"
prior.brm.2 <- brm(bf.1,data=gridded_df_ordered,
                   cores = 4,
                   prior = prior.2,
                   sample_prior = "only")

pp_check(prior.brm.2)+thm1




prior.3 <- prior.1
prior.3$prior[1] <- "normal(70,10)"
prior.3$prior[2:3] <- "cauchy(0,1)"
prior.brm.3 <- brm(bf.1,data=gridded_df_ordered,
                   cores = 4,
                   prior = prior.3,
                   sample_prior = "only")

pp_check(prior.brm.3)+thm1





gridded_df_ordered$repID <- as.numeric(gsub("R","",gridded_df_ordered$rep))


Nsim <- nrow(gridded_df_ordered)

set.seed(seed = 1923840479)
tau0  <- abs(rnorm(1, 0, 1))
sigma <- abs(rnorm(1, 0, 1))

beta0 <- rnorm(median(gridded_df_ordered$yield), 0, 1)
beta1 <- rnorm(1, 0, 1)
beta2 <- rnorm(1, 0, 1)
beta0i <- rnorm(3, 0, tau0)


data2 <- data.frame(
  raw = gridded_df_ordered$yield,
  sim = beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
)


ggplot(data2, aes(x = raw, y = sim)) +
  geom_point(alpha = 0.1)


set.seed(seed = 1923840479)
tau0  <- 1 / sqrt(rgamma(1, 1, rate = 100))
sigma <- 1 / sqrt(rgamma(1, 1, rate = 100))

beta0 <- rnorm(median(gridded_df_ordered$yield), 0, 1)
beta1 <- rnorm(1, 0, 1)
beta2 <- rnorm(1, 0, 1)
beta0i <- rnorm(3, 0, tau0)


data1 <- data.frame(
  raw = gridded_df_ordered$yield,
  sim = beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
)

ggplot(data1, aes(x = raw, y = sim)) +
  geom_point(alpha = 0.1)


set.seed(seed = 1923840479)

y_prior <- matrix(0,50,N)
for(i in 1:50){
  
  tau0  <- 1 / sqrt(rgamma(1, 1, rate = 100))
  sigma <- 1 / sqrt(rgamma(1, 1, rate = 100))
  beta0 <- rnorm(median(gridded_df_ordered$yield), 0, 100)
  beta1 <- rnorm(1, 0, 100)
  beta2 <- rnorm(1, 0, 100)
  beta0i <- rnorm(3, 0, tau0)
  y_prior[i,] <- beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
}
pp_check(standata.model$Y,y_prior,fun = "ppc_dens_overlay") +
  coord_cartesian(xlim = c(-1000, 2000))


y_prior1 <- matrix(0,50,N)
for(i in 1:50){
  
  tau0  <- abs(rstudent(1,3, 0, 30))
  sigma <- abs(rstudent(1,3, 0, 30))
  beta0 <- rstudent(1,3,median(gridded_df_ordered$yield),30)
  beta1 <- rnorm(1,0, 0.01)
  beta2 <- rnorm(1,0, 0.01)
  beta0i <- rstudent(3,3, 0, tau0)
  
  y_prior1[i,] <-beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
}
pp_check(standata.model$Y,y_prior1,fun = "ppc_dens_overlay") +
  coord_cartesian(xlim = c(-1000, 2000))

plot(standata.modelar1$Y,y_prior1[1,],ylim=c(-6000,8000))
for(i in 2:50) points(standata.modelar1$Y,y_prior1[i,])




y_prior2 <- matrix(0,50,N)
for(i in 1:50){
  
  tau0  <- abs(rnorm(1, 0, 10))
  sigma <- abs(rnorm(1, 0, 10))
  beta0 <- rnorm(median(gridded_df_ordered$yield), 0, 10)
  beta1 <- rnorm(1, 0, 10)
  beta2 <- rnorm(1, 0, 10)
  beta0i <- rnorm(3, 0, tau0)
  
  y_prior2[i,] <-beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
}
pp_check(standata.model$Y,y_prior2,fun = "ppc_dens_overlay") +
  coord_cartesian(xlim = c(-1000, 2000))

plot(standata.modelar1$Y,y_prior2[1,],ylim=c(-6000,8000))
for(i in 2:50) points(standata.modelar1$Y,y_prior2[i,])



plot(y_prior[1,])
plot(y_prior1[1,])
plot(y_prior2[1,])


y_prior3 <- matrix(0,50,N)
for(i in 1:50){
  
  tau0  <- abs(rstudent_t(1,3, 0, 31.3))
  sigma <- abs(rstudent_t(1,3, 0, 31.3))
  beta0 <- rstudent_t(1,3,median(gridded_df_ordered$yield), 31.3)
  beta1 <- rnorm(1, 0, 1)
  beta2 <- rnorm(1, 0, 1)
  beta0i <- rstudent_t(3,3, 0, tau0)
  
  y_prior3[i,] <-beta0 + beta1*gridded_df_ordered$cnitro+
    beta2*gridded_df_ordered$cnitro.sq+
    beta0i[gridded_df_ordered$repID]+
    rnorm(Nsim, mean = 0, sd = sigma)
}
pp_check(standata.model$Y,y_prior3,fun = "ppc_dens_overlay") +
  coord_cartesian(xlim = c(-1000, 2000))

plot(standata.modelar1$Y,y_prior3[1,],ylim=c(-6000,8000))
for(i in 2:50) points(standata.modelar1$Y,y_prior3[i,])



data5 <- data.frame(y=c(gridded_df_ordered$yield,
                        y_prior,y_prior1,y_prior2,y_prior3),
                    x=c(rep("Raw",N),rep(c("P","P1","P2","P3"),each=83700)))

ggplot(data5) + geom_density_ridges_gradient(aes(y,x,fill=x))

ggplot(data5) + geom_density(aes(y,colour=x,fill=x),size=1,alpha=0.2) + 
  coord_cartesian(xlim=c(-1e3,1e3))


rs <- rstudent_t(1e5,3, 0, 31.3)
dens(rs)

rs <- rstudent_t(1e5,3, 0, 10)
dens(rs)


# Plot: prior predictive with weakly informative priors
set.seed(seed = 192384)

N <- standata.model$N
Xc <- matrix(0,N,2)
means_X <- numeric(2)
for (i in 2:3) {
  means_X[i - 1] = mean(standata.model$X[, i]);
  Xc[, i - 1] = standata.model$X[, i] - means_X[i - 1];
}

i <- 1
j <- 1
M <- 50


y_prior <- rnorm(rstudent_t(1, 3, 84.7, 31.3) + Xc%*%rnorm(2) + 
                   rstudent_t(3, 3, 0, 31.3)[standata.model$J_1],
                 abs(rstudent_t(1, 3, 0, 31.3)))

plot(standata.model$Y,y_prior)


y_prior2 <- rnorm(rstudent_t(1, 3, 84.7, 1) + Xc%*%rnorm(2) + 
                    rstudent_t(3, 3, 0, 1)[standata.model$J_1],
                  abs(rstudent_t(1, 3, 0, 1)))

plot(standata.model$Y,y_prior2)

plot(standata.model$Y,y_prior2,ylim=c(-5,30))
points(standata.model$Y,y_prior,col="red")



y_prior3 <- rnorm(rnorm(1,84.7, 100) + Xc%*%rnorm(2) + 
                    rnorm(3, 0, 100)[standata.model$J_1],
                  abs(rnorm(1,0, 100)))

plot(standata.model$Y,y_prior3)


y_prior <- matrix(0,50,N)
for(i in 1:50){
  y_prior[i,] <- rnorm(rstudent_t(1, 3, 84.7, 30) + Xc%*%rnorm(2,0,0.01)  + 
                         rstudent_t(3, 3, 0, 1)[standata.model$J_1],
                       abs(rstudent_t(1, 3, 0, 30)))
}
pp_check(standata.model$Y,y_prior,fun = "ppc_dens_overlay")


y_prior2 <- matrix(0,50,N)
for(i in 1:50){
  y_prior2[i,] <- rnorm(rstudent_t(1, 3, 84.7, 30) + Xc%*%rnorm(2,0,0.01)  + 
                          rstudent_t(3, 3, 0, 1)[standata.model$J_1],
                        abs(rstudent_t(1, 3, 0, 30)))
}
pp_check(standata.model$Y,y_prior2,fun = "ppc_dens_overlay")


y_prior3 <- matrix(0,50,N)
for(i in 1:50){
  y_prior3[i,] <- rnorm(rnorm(1,84.7, 1) + Xc%*%rnorm(2,0,0.01) + 
                          rnorm(3, 0, 100)[standata.model$J_1],
                        abs(rnorm(1,0, 100)))
}
pp_check(standata.model$Y,y_prior3,fun = "ppc_dens_overlay")





tau0 <- abs(rnorm(1, 0, 1))
tau1 <- abs(rnorm(1, 0, 1))
sigma <- abs(rnorm(1, 0, 1))
beta0i <- rnorm(8, 0, tau0)
beta1i <- rnorm(8, 0, tau1)
beta0 <- rnorm(1, 0, 1)
beta1 <- rnorm(1, 1, 1)

data2 <- data.frame(
  log_pm25 = GM$log_pm25,
  sim = beta0 + beta0i[GM$super_region] +
    (beta1 + beta1i[GM$super_region]) * GM$log_sat +
    rnorm(Nsim, mean = 0, sd = sigma)
)

ggplot(data2, aes(x = log_pm25, y = sim)) +
  geom_point(alpha = 0.1) + 
  xysim_labs


data3 <- data.frame(
  x = standata.modelar1$Y, 
  wip = y_prior[1,], 
  vague = y_prior[2,]
)
ggplot(data3, aes(x=x, y=wip)) + 
  geom_point(alpha = 0.1) + 
  geom_point(
    aes(y = vague), 
    color = "red", 
    alpha = 0.1
  )








########## prior check ############


form.null <- brmsformula(yield ~ 1)

# make_stancode(form.null,data=gridded_df_ordered)
prior.1 <- get_prior(form.null,data=gridded_df_ordered)

prior.1 <- c(set_prior("normal(0,1)", class = "b"),
             set_prior("student_t(3, 84.7, 10)", class = "Intercept"),
             set_prior("student_t(3, 0, 10)", class = "sigma"))

mod.null <- brm(form.null,data=gridded_df_ordered,
                sample_prior = "only",
                prior = prior.1,
                cores = 4)

pp_check(mod.null)

summary(mod.null)

prior.2 <- c(set_prior("normal(0,1)", class = "b"),
             set_prior("student_t(3, 84.7, 100)", class = "Intercept"),
             set_prior("student_t(3, 0, 100)", class = "sigma"))

mod.null.2 <- brm(form.null,data=gridded_df_ordered,
                  sample_prior = "only",
                  prior = prior.2,
                  cores = 4)

pp_check(mod.null.2)



form.brm <- brmsformula(yield ~ 1+cnitro + cnitro.sq + (1|rep))

standata.model <- make_standata(form.brm,data=gridded_df_ordered)
N <- nrow(gridded_df_ordered)
y_prior <- matrix(0,50,N)
for(i in 1:50){
  y_prior[i,] <- rnorm(rstudent_t(1, 3, 84.7, 10) + standata.model$X,
                       abs(rstudent_t(1, 3, 0, 5)))
}
pp_check(standata.model$Y,y_prior,fun = "ppc_dens_overlay")



post.1 <- as.data.frame(extract(mod.null$fit))
names(post.1)

rs <- abs(rstudent_t(1e3,3,0, 1))

DF <- data.frame(x = c(post.1$sigma,rs),
                 y = c(rep("post",length(post.1$sigma)),
                       rep("pr",length(rs))))

ggplot(DF, aes(x = x, fill = y)) +
  geom_density(alpha = 0.3)



rs <- abs(rstudent_t(1e3,10,84.7, 10))

DF <- data.frame(x = c(post.1$b_Intercept,rs),
                 y = c(rep("post",length(post.1$b_Intercept)),
                       rep("pr",length(rs))))

ggplot(DF, aes(x = x, fill = y)) +
  geom_density(alpha = 0.3)



prior.2 <- c(set_prior("normal(0,1)", class = "b"),
             set_prior("student_t(3, 0, 10)", class = "sd"),
             set_prior("student_t(3, 84.7, 10)", class = "Intercept"),
             set_prior("student_t(3, 0, 10)", class = "sigma"))

# prior.1$prior[1:3] <- "normal(0,1)"

mod.null.2 <- brm(form.null,data=gridded_df_ordered,
                  # sample_prior = "only",
                  prior = prior.2,
                  cores = 4)

pp_check(mod.null.2)
pp_check(mod.null)

bayes_R2(mod.null)
bayes_R2(mod.null.2)





pp_check(mod.null) + thm1 + coord_cartesian(xlim=c(-200,500))


make_stancode(form.null,data=gridded_df_ordered,
              prior = prior.2)


form.1 <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep))
mod.1 <- brm(form.1,data=gridded_df_ordered,
             cores = 4)



prior.1 <- get_prior(form.1,data=gridded_df_ordered)
prior.1$prior[1:3] <- "normal(0,0.01)"
mod.1 <- brm(form.1,data=gridded_df_ordered,
             prior = prior.1,
             sample_prior = "only", 
             cores = 4)

me_loss_prior1 <- marginal_effects(
  mod.1, conditions = data.frame(y=gridded_df_ordered$topo), 
  re_formula = NULL, method = "predict"
)
plot(me_loss_prior1,points =TRUE)


pp_check(mod.1) + thm1
summary(mod.1)


prior.2 <- prior.1
prior.2$prior[5:7] <- "student_t(3, 0, 1)"
mod.2 <- brm(form.1,data=gridded_df_ordered,
             prior = prior.2,
             sample_prior = "only", 
             cores = 4)

pp_check(mod.2) + thm1
summary(mod.2)




prior.3 <- prior.1

prior.3$prior[1:3] <- "normal(0,0.01)"
prior.3$prior[4]   <- "student_t(3, 84.7, 30)"
prior.3$prior[5:8] <- "student_t(3, 0, 30)"
prior.3$bound[5:8] <- "0"

mod.1.2 <- brm(form.1,data=gridded_df_ordered,
               prior = prior.3,
               sample_prior = "only", 
               cores = 4)

pp_check(mod.1.2) + thm1

me_loss_prior2 <- marginal_effects(
  mod.1.2, conditions = data.frame(y=gridded_df_ordered$topo), 
  re_formula = NULL, method = "predict"
)
plot(me_loss_prior2,points =TRUE)


make_stancode(form.1,data=gridded_df_ordered)


mcmc_areas(
  as.array(mod.1.2), 
  pars = c("b_Intercept"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)

mcmc_areas(
  as.array(mod.1.2), 
  pars = c("b_cnitro",
           "b_cnitro.sq"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)


mcmc_areas(
  as.array(mod.1.2), 
  pars = c("sd_rep__Intercept",
           "r_rep[R1,Intercept]", "r_rep[R2,Intercept]",
           "r_rep[R3,Intercept]","sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)

mod.1.update <- update(mod.1.2, sample_prior="no", seed = 1234,
                       cores=4)

summary(mod.1.update)

pp_check(mod.1.update)
pp_check(mod.1)

summary(mod.1)

bayes_R2(mod.1.update)
bayes_R2(mod.1)




form.3 <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep)+
                        (1 + cnitro + cnitro.sq|gridId))

prior.3<- get_prior(form.3,data=gridded_df_ordered)

prior.3$prior[1:3] <- "normal(0,0.01)"
prior.3$prior[4:5] <- "lkj(3)"
prior.3$prior[9:13] <- "normal(0,0.01)"


mod.3 <- brm(form.3,data=gridded_df_ordered,
             prior = prior.3,
             sample_prior = "only",
             cores = 4)
pp_check(mod.3)




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
  ylab("Density")










form.4 <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep)+s(x,y))

prior.4<- get_prior(form.4,data=gridded_df_ordered)

mod.4 <- brm(form.4,data=gridded_df_ordered,
             #prior = prior.1,
             cores = 4)

pp_check(mod.4) + thm1
summary(mod.4)




form.5 <- brmsformula(yield ~ 1 + cnitro + cnitro.sq + (1|rep)+
                        (1 + cnitro + cnitro.sq |gridId) + 
                        s(x,y))

prior.5<- get_prior(form.5,data=gridded_df_ordered)

prior.5$prior[1:5] <- "normal(1,10)"

mod.5 <- brm(form.5,data=gridded_df_ordered,
             prior = prior.5,
             # sample_prior = "only",
             chains = 1)

pp_check(mod.5)




fit.brm.5 <- mod.5$fit

traceplot(fit.brm.5,"sigma")


sigs <- as.data.frame(extract(fit.brm.5,"sigma"))
hist(sigs$sigma,30)

dens(sigs$sigma)

rs <- rstudent_t(1e3,3, 5, 5) 
dens(rs)

DF <- data.frame(x = c(sigs$sigma,rs),
                 y = c(rep("post",1000),rep("pr",1e3)))

ggplot(DF, aes(x = x, fill = y)) +
  geom_density(alpha = 0.3) + 
  coord_cartesian(xlim = c(-5,10))


coord_cartesian(xlim = c(-5,10)) +
  labs(x = "beta")




fit.as <- asreml(yield ~ 1+cnitro+cnitro.sq,
                 random = ~ rep,
                 residual = ~ar1(row):ar1(col),
                 data=gridded_df_ordered)
wald(fit.as)
summary(fit.as)
plot(varioGram(fit.as))


fit.as.1 <- asreml(yield ~ 1+cnitro+cnitro.sq,
                   random = ~ topo,
                   residual = ~ar1(row):ar1(col),
                   data=gridded_df_ordered)
wald(fit.as.1)

summary(fit.as.1)
plot(varioGram(fit.as.1))

lrt(fit.as.1,fit.as)




fit.as.2 <- asreml(yield ~ 1+cnitro+cnitro.sq,
                   random = ~ topo+rep,
                   residual = ~ar1(row):ar1(col),
                   data=gridded_df_ordered)
summary(fit.as.2)

lrt(fit.as.1,fit.as,fit.as.2)







ggplot(gridded_df_ordered, aes(y = yield, x = cnitro)) +
  geom_point(
    aes(color = rep), 
    alpha = 0.4,
    size = rel(0.8)
  ) + 
  geom_smooth(
    # formula = y ~ x,
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  )



ggplot(gridded_df_ordered, aes(y = yield, x = cnitro)) +
  geom_point(aes(colour = rep), alpha = 0.2, size = rel(0.75)) + 
  geom_smooth(
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  ) + 
  scale_color_manual(
    values = 
      c("R1" = "#00C094",
        "R2" = "#FB61D7",
        "R3" = "#53B400")
  ) +
  geom_smooth(aes(colour = rep), method = lm) + 
  coord_equal() +
  legend_none() 

ggplot(gridded_df_ordered, aes(y = yield, x = cnitro)) +
  geom_point(aes(colour = topo), alpha = 0.2, size = rel(0.75)) + 
  geom_smooth(
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  ) + 
  scale_color_manual(
    values = 
      c("E" = "#00C094",
        "HT" = "#FB61D7",
        "LO" = "#53B400",
        "W" = "#F8766D")
  ) +
  geom_smooth(aes(colour = topo), method = lm) + 
  coord_equal() +
  legend_none() 
















