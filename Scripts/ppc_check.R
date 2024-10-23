
########## RRC without rep #########
form.randU <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                            (1 + nitro + nitro.sq | gridId))

data.randU <- make_standata(form.randU,data=gridded_df_ordered)
data.randU$nrow <- 93
data.randU$ncol <- 18

stanRandU <- readRDS("Outcome/stanRandU_norep_pawsey.rds")
post.rand <- as.array(rstan::extract(stanRandU))

N <- data.randU$N
i <- 1
j <- 1
M <- length(post.rand$Intercept)
y_rep_randU <- matrix(0,M,N)
loglik_randU <- matrix(0,M,N)

mus.rand <- matrix(0,M,N)
set.seed(128374)
for(i in 1:M){
  mus.rand[i,] <- post.rand$b_Intercept[i] + 
    data.randU$X[,2:3]%*%post.rand$b[i,] +
    as.numeric(post.rand$r_1_1[i,]*data.randU$Z_1_1) + 
    as.numeric(post.rand$r_1_2[i,]*data.randU$Z_1_2) +
    as.numeric(post.rand$r_1_3[i,]*data.randU$Z_1_3)
  
  y_rep_randU[i,] <- rnorm(N,mus.rand[i,],abs(post.rand$sigma[i]))
  
  loglik_randU[i,] <- dnorm(data.randU$Y,
                            mus.rand[i,],abs(post.rand$sigma[i]),
                            log = TRUE)
  
  cat(paste0(round(i/M*100,2),"% "))
}
pp_check(data.randU$Y, y_rep_randU[sample(1:M,100),], 
         ppc_dens_overlay) + thm1


var_mu.rand <- apply(mus.rand, 1, var)
sigma2.rand <- post.rand$sigma^2
rsq_bayes.rand <- var_mu.rand / (var_mu.rand + sigma2.rand)

summary(rsq_bayes.rand)
quantile(rsq_bayes.rand,c(0.025,0.5,0.975))


loo1.randU <- loo(loglik_randU)
psis.randU <- psis(-loglik_randU)
plot(loo1.randU, label_points = TRUE)

pareto_k_table(loo1.randU)


pit1.randU <- loo_pit(object = y_rep_randU,
                      y = data.randU$Y, 
                      lw = weights(psis.randU))

unifs.randU <- matrix(runif(length(pit1.randU) * 100), nrow = 100)
ppc_dens_overlay(pit1.randU, unifs.randU)


ppc_loo_pit_qq(data.randU$Y, 
               y_rep_randU, lw = weights(psis.randU))


loo1.randU <- loo(loglik_randU)
psis_result <- psis(-loglik_randU)
plot(psis_result, label_points = TRUE)

pareto_k_table(psis_result)


pp_check(gridded_df_ordered$yield,y_rep_randU[sample(1:M,100),],
         ppc_dens_overlay)

ppc_loo_pit_overlay(
  y = data.randU$Y,
  yrep = y_rep_randU,
  lw = weights(psis_result)
)

ppc_loo_pit_qq(
  y = data.randU$Y,
  yrep = y_rep_randU,
  lw = weights(psis_result)
)

ppc_loo_pit_qq(
  y = data.randU$Y,
  yrep = y_rep_randU,
  lw = weights(psis_result),
  compare = "normal"
)




psis_result.2 <- psislw(-loglik_randU)

pit1.randU <- rstantools::loo_pit(object = y_rep_randU,
                      y = data.randU$Y, 
                      lw = psis_result.2$lw_smooth)
unifs.randU <- matrix(runif(length(pit1.randU) * 100), nrow = 100)
ppc_dens_overlay(pit1.randU, unifs.randU)

ppc_loo_pit_qq(data.randU$Y, 
               y_rep_randU, lw = psis_result.2$lw_smooth)





loo1.rep <- loo(loglik_randU.rep)
psis.rep <- psis(-loglik_randU.rep)
plot(loo1.rep, label_points = TRUE)

pareto_k_table(loo1.rep)

pit1.rep <- loo_pit(object = y_rep_randU.rep,
                    y = datarand.rep$Y, 
                    lw = weights(psis.rep))
unifs.rep <- matrix(runif(length(pit1.rep) * 100), nrow = 100)
ppc_dens_overlay(pit1.rep, unifs.rep)

ppc_loo_pit_qq(datarand.rep$Y,
               y_rep_randU.rep, lw = weights(psis.rep))

ppc_loo_pit_overlay(
  y = datarand.rep$Y,
  yrep = y_rep_randU.rep,
  lw = weights(psis.rep)
)




loo1.full <- loo(loglik_full)
psis_result.full <- psis(-loglik_full)
plot(psis_result.full, label_points = TRUE)

pareto_k_table(psis_result.full)

pp_check(gridded_df_ordered$yield,y_rep_full[sample(1:M,100),],
         ppc_dens_overlay)

ppc_loo_pit_overlay(
  y = datafull$Y,
  yrep = y_rep_full,
  lw = weights(psis_result.full)
)




########### outers <- c(74,1116,1490) ############## 

outers <- c(74,1116,1490)
gridded_out3 <-gridded_df_ordered
gridded_out3$outer <- 0
gridded_out3$outer[outers] <- 1
gridded_out3[outers,]

table(gridded_out3$outer)

form.randU.out <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                                (1 + nitro + nitro.sq | gridId) + 
                                (1|outer))
data.randU.out3 <- make_standata(form.randU.out,
                                 data=gridded_out3)
data.randU.out3$nrow <- 93
data.randU.out3$ncol <- 18


stanRandUout3 <- readRDS("Outcome/stanRandUout3_pawsey.rds")

post.rand.out3 <- as.array(rstan::extract(stanRandUout3))

N <- data.randU.out3$N
i <- 1
j <- 1
M <- length(post.rand.out3$Intercept)
y_rep_randU_out3 <- matrix(0,M,N)
loglik_randU_out3 <- matrix(0,M,N)

mus.out3 <- matrix(0,M,N)
set.seed(128374)
for(i in 1:M){
  mus.out3[i,] <- post.rand.out3$b_Intercept[i] + 
    data.randU.out3$X[,2:3]%*%post.rand.out3$b[i,] +
    as.numeric(post.rand.out3$r_1_1[i,]*data.randU.out3$Z_1_1) + 
    as.numeric(post.rand.out3$r_1_2[i,]*data.randU.out3$Z_1_2) +
    as.numeric(post.rand.out3$r_1_3[i,]*data.randU.out3$Z_1_3)
  
  for (j in 1:N) {
    mus.out3[i,j]<- mus.out3[i,j]+ 
      post.rand.out3$r_2_1[i,data.randU.out3$J_2[j]]
  }
  y_rep_randU_out3[i,] <- rnorm(N,mus.out3[i,],abs(post.rand.out3$sigma[i]))
  loglik_randU_out3[i,] <- dnorm(data.randU.out3$Y,
                                 mus.out3[i,],
                                 abs(post.rand.out3$sigma[i]),
                                 log=TRUE)
  
  cat(paste0(round(i/M*100,2),"% "))
}
pp_check(data.randU.out3$Y, y_rep_randU_out3[sample(1:M,100),], 
         ppc_dens_overlay) + thm1



var_mu.out3 <- apply(mus.out3, 1, var)
sigma2.out3 <- post.rand.out3$sigma^2
rsq_bayes.out3 <- var_mu.out3 / (var_mu.out3 + sigma2.out3)

summary(rsq_bayes.out3)
quantile(rsq_bayes.out3,c(0.025,0.5,0.975))


loo1.out3 <- loo(loglik_randU_out3)
psis.out3 <- psis(-loglik_randU_out3)
plot(loo1.out3, label_points = TRUE)

pareto_k_table(loo1.out3)

pit1.out3 <- loo_pit(object = y_rep_randU_out3,
                     y = data.randU.out3$Y, 
                     lw = weights(psis.out3))
unifs.out3 <- matrix(runif(length(pit1.out3) * 100), nrow = 100)
ppc_dens_overlay(pit1.out3, unifs.out3) + legend_none()

ppc_loo_pit_qq(data.randU.out3$Y, 
               y_rep_randU_out3, lw = weights(psis.out3))



########## full model ##############

stanFULL <- readRDS("Outcome/StanFullModel_PC.rds")
post.full <- as.array(rstan::extract(stanFULL))

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


var_mu.full <- apply(mus.full, 1, var)
sigma2.full <- apply(sigma.full, 1, var)
rsq_bayes.full <- var_mu.full / (var_mu.full + sigma2.full)

summary(rsq_bayes.full)
quantile(rsq_bayes.full,c(0.025,0.5,0.975))


loo1.full <- loo(loglik_full)
psis.full <- psis(-loglik_full)
plot(loo1.full, label_points = TRUE)

pareto_k_table(loo1.full)

pit1.full <- loo_pit(object = y_rep_full,
                     y = datafull$Y, 
                     lw = weights(psis.full))
unifs.full <- matrix(runif(length(pit1.full) * 100), nrow = 100)
ppc_dens_overlay(pit1.full, unifs.full) + legend_none()


ppc_loo_pit_qq(datafull$Y, 
               y_rep_full, lw = weights(psis.full))



############ rand U with rep #############
## seeds  1072473474
stanRandU.rep <- readRDS("Outcome/StanRandUAR1_PC.rds")
post.rand.rep <- as.array(rstan::extract(stanRandU.rep))

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



var_mu.rep <- apply(mus.rep, 1, var)
sigma2.rep <- post.rand.rep$sigma^2
rsq_bayes.rep <- var_mu.rep / (var_mu.rep + sigma2.rep)

summary(rsq_bayes.rep)
quantile(rsq_bayes.rep,c(0.025,0.5,0.975))


loo1.rep <- loo(loglik_randU.rep)
psis.rep <- psis(-loglik_randU.rep)
plot(loo1.rep, label_points = TRUE)

pareto_k_table(loo1.rep)

pit1.rep <- loo_pit(object = y_rep_randU.rep,
                     y = datarand.rep$Y, 
                     lw = weights(psis.rep))
unifs.rep <- matrix(runif(length(pit1.rep) * 100), nrow = 100)
ppc_dens_overlay(pit1.rep, unifs.rep) + legend_none()

ppc_loo_pit_qq(datarand.rep$Y,
               y_rep_randU.rep, lw = weights(psis.rep))




##################

ggplot(gridded_df_ordered) + 
  geom_histogram(aes(yield,fill=topo))

ggplot(gridded_df_ordered) + 
  geom_histogram(aes(yield,fill=rep))

ggplot(gridded_df_ordered) + 
  geom_histogram(aes(yield)) + 
  facet_wrap(~rep)

ggplot(gridded_df_ordered) + 
  geom_histogram(aes(yield)) + 
  facet_wrap(~topo)


plot(gridded_df_ordered$bv,gridded_df_ordered$yield)


ggplot(gridded_df_ordered) + 
  geom_point(aes(x,y,colour=bv))

ggplot(sub.corn) + geom_point(aes(long,lat,colour=bv)) + 
  xlab("Longitude") + ylab("Latitude")
  





