source("Initialization.R")
form.rand <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                           (1 + nitro + nitro.sq | gridId))
data.rand <- make_standata(form.rand,
                           data=gridded_df_ordered)
data.rand$nrow <- 93
data.rand$ncol <- 18

nuts_controls   <- list(max_treedepth = 15, adapt_delta = 0.99)

stan_ST_Rand_NS <- stan("StanModel/ST_rand_NoSpatial.stan",
                        data  = data.rand,
                        control = nuts_controls,
                        cores = 4)
# saveRDS(stan_ST_Rand_NS,"stan_ST_Rand_NS_lap1599.rds")

stan_GS_Rand_NS <- stan("StanModel/GS_rand_NoSpatial.stan",
                        data  = data.rand,
                        control = nuts_controls,
                        cores = 4)
# saveRDS(stan_GS_Rand_NS,"stan_GS_Rand_NS_lap1599.rds")

stan_ST_Rand <- stan("StanModel/ST_rand_faster_try2.stan",
                     data = data.rand,
                     control = nuts_controls,
                     cores = 4)
# saveRDS(stan_ST_Rand,"stan_ST_Rand_lap1599.rds")

stan_GS_Rand <- stan("StanModel/GS_rand_faster_try2.stan",
                     data = data.rand,
                     control = nuts_controls,
                     cores = 4)
# saveRDS(stan_GS_Rand,"Outcome/stan_GS_Rand_lap1599.rds")

########### AR1 x AR1 residuals ###########

nuts_controls   <- list(max_treedepth = 10, adapt_delta = 0.90)

stan_AR1 <- stan("StanModel/GS_AR1.stan",
                 data = data.rand,
                 control = nuts_controls,
                 cores = 4)
# saveRDS(stan_AR1,"Outcome/stan_AR1_1090.rds")





############ further ################

source("Initialization.R")
# nuts_controls   <- list(max_treedepth = 15, adapt_delta = 0.99)
nuts_controls   <- list(max_treedepth = 12, adapt_delta = 0.95)

gridded_df_ordered$gridId <- as.factor(gridded_df_ordered$gridId)

form.rand <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                           (1 + nitro + nitro.sq | gridId))
data.rand <- make_standata(form.rand,
                           data=gridded_df_ordered)
data.rand$nrow <- 93
data.rand$ncol <- 18

## order columns 18 nested within rows 93
stan_GS_SP_1295_try2 <- stan("StanModel/GS_rand_faster_try2.stan",
                        data = data.rand,
                        control = nuts_controls,
                        cores = 4)
saveRDS(stan_GS_SP_1295_try2,"stan_GS_SP_1295_try2_PC.rds")

stan_ST_SP_1295_try2 <- stan("StanModel/ST_rand_faster_try2.stan",
                        data = data.rand,
                        control = nuts_controls,
                        cores = 4)
saveRDS(stan_ST_SP_1295_try2,"stan_ST_SP_1295_try2_PC.rds")





DistMat <- as.matrix(dist(gridded_df_ordered[c("row","col")], 
                method = "euclidean", diag = FALSE, upper = FALSE, p = 2))

data.rand$Distmat <- DistMat
stan_Matern <- stan("StanModel/GS_simu_matern.stan",
                    data = data.rand,
                    iter = 10,
                    # control = nuts_controls,
                    cores = 4)
# saveRDS(stan_Matern,"stan_Matern.rds")

stan_Matern <- stan("StanModel/GS_simu_matern_V2.stan",
                    data = data.rand,
                    iter = 10,
                    # control = nuts_controls,
                    cores = 4)



## order rows 93 nested within columns 18
dat_neworder <- gridded_df_ordered[order(gridded_df_ordered$col,gridded_df_ordered$row),]

dat_neworder$gridId <- as.factor(1:nrow(dat_neworder))

form.rand <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                           (1 + nitro + nitro.sq | gridId))
data.rand <- make_standata(form.rand,
                           data=dat_neworder)
data.rand$nrow <- 93
data.rand$ncol <- 18

stan_GS_SP_try1 <- stan("StanModel/GS_rand_faster_try1.stan",
                        data = data.rand,
                        # control = nuts_controls,
                        iter = 10,
                        cores = 4)



dat_neworder <- gridded_df_ordered[order(gridded_df_ordered$col,gridded_df_ordered$row),]

dat_neworder$gridId <- as.factor(1:nrow(dat_neworder))

form.rand <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                           (1 + nitro + nitro.sq | gridId))
data.rand <- make_standata(form.rand,
                           data=dat_neworder)
data.rand$nrow <- 93
data.rand$ncol <- 18


# nuts_controls   <- list(max_treedepth = 15, adapt_delta = 0.99)


# stan_ST_Rand_NS <- stan("StanModel/ST_rand_NoSpatial.stan",
#                         data  = data.rand,
#                         control = nuts_controls,
#                         cores = 4)
# # saveRDS(stan_ST_Rand_NS,"stan_ST_Rand_NS_lap1599.rds")
# 
# stan_GS_Rand_NS <- stan("StanModel/GS_rand_NoSpatial.stan",
#                         data  = data.rand,
#                         control = nuts_controls,
#                         cores = 4)
# # saveRDS(stan_GS_Rand_NS,"stan_GS_Rand_NS_lap1599.rds")



stan_ST_Rand <- stan("StanModel/ST_rand_faster_further.stan",
                     data = data.rand,
                     # control = nuts_controls,
                     iter = 10,
                     cores = 4)
# saveRDS(stan_ST_Rand,"stan_ST_Rand_lap1599.rds")



stan_ST_Rand <- stan("StanModel/ar1Urand_NoRep_ST.stan",
                     data = data.rand,
                     # control = nuts_controls,
                     iter = 10,
                     cores = 4)

stan_GS_SP_try <- stan(model_code = gsmodel,
                     data = data.rand,
                     # control = nuts_controls,
                     # iter = 10,
                     cores = 4)

stan_ST_SP_try <- stan(model_code = stmodel,
                       data = data.rand,
                       # control = nuts_controls,
                       # iter = 10,
                       cores = 4)


# saveRDS(stan_GS_SP_try,"Outcome/stan_GS_SP_try.rds")

# stan_GS_SP_try <- readRDS("Outcome/stan_GS_SP_try")

check_hmc_diagnostics(stan_GS_SP_try)

post.gs.try <- extract(stan_GS_SP_try)

M <- length(post.gs.try$Intercept)

pp_check(gridded_df_ordered$yield,
         post.gs.try$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_GSRand.pdf",width=11.96,height=8.27)


loo.gs.try  <- loo(post.gs.try$log_lik)
psis.gs.try<- psis(-post.gs.try$log_lik)

rsq_bayes(post.gs.try)


stan_GS_SP_try1 <- stan("StanModel/GS_rand_faster_try1.stan",
                       data = data.rand,
                       # control = nuts_controls,
                       # iter = 10,
                       cores = 4)

# stan_GS_SP_try1 <- stan(model_code = stantry1,
#                         data = data.rand,
#                         # control = nuts_controls,
#                         # iter = 10,
#                         cores = 4)
saveRDS(stan_GS_SP_try1,"stan_GS_SP_try1.rds")
# stan_GS_Rand <- stan("StanModel/GS_rand_faster.stan",
#                      data = data.rand,
#                      control = nuts_controls,
#                      cores = 4)
# # saveRDS(stan_GS_Rand,"Outcome/stan_GS_Rand_lap1599.rds")

stan_ST_SP_try1 <- stan("StanModel/ST_rand_faster_try1.stan",
                        data = data.rand,
                        # control = nuts_controls,
                        cores = 4)
saveRDS(stan_ST_SP_try1,"stan_ST_SP_try1.rds")


library(cmdstanr)
library(posterior)
library(bayesplot)

file <- file.path("StanModel/ar1Urand_NoRep_ST.stan")
mod <- cmdstan_model(file)

mod$print()


dat_neworder <- gridded_df_ordered[order(gridded_df_ordered$col,gridded_df_ordered$row),]

dat_neworder$gridId <- as.factor(1:nrow(dat_neworder))

form.rand <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                           (1 + nitro + nitro.sq | gridId))



data.rand <- make_standata(form.rand,
                           data=dat_neworder)

make_stancode(form.rand,data=dat_neworder,
              family = student)

data.rand$nrow <- 93
data.rand$ncol <- 18

fit <- mod$sample(
  data = data.rand,
  seed = 2021,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)



######### topo testing ##########

form.topo <- brmsformula(yield ~ 1 + nitro + nitro.sq + (1|topo) + 
                           (1 + nitro + nitro.sq | gridId))
data.topo <- make_standata(form.topo,
                           data=gridded_df_ordered)
data.topo$nrow <- 93
data.topo$ncol <- 18

nuts_controls   <- list(max_treedepth = 12, adapt_delta = 0.95)


stan_GS_Rand_Topo <- stan("StanModel/GS_rand_topo_try2.stan",
                          data  = data.topo,
                          control = nuts_controls,
                          cores = 4)


stan_GS_Rand_Topo <- readRDS("Outcome/stan_GS_Rand_Topo.rds")

check_hmc_diagnostics(stan_GS_Rand_Topo)

traceplot(stan_GS_Rand_Topo,c("b"))

post.gs.rand.topo <- extract(stan_GS_Rand_Topo)

M <- length(post.gs.rand.topo$Intercept)

pp_check(gridded_df_ordered$yield,
         post.gs.rand.topo$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_GSRNS.pdf",width=11.96,height=8.27)


loo.gs.rand.topo  <- loo(post.gs.rand.topo$log_lik)
psis.gs.rand.topo <- psis(-post.gs.rand.topo$log_lik)

rsq_bayes(post.gs.rand.topo)



sumstan(post.gs.rand.topo$b)
sumstan(post.gs.rand.topo$rho_r)
sumstan(post.gs.rand.topo$rho_c)




gridded_df_ordered$topo <- as.factor(gridded_df_ordered$topo)
form.topo2 <- brmsformula(yield ~ 1 + nitro + nitro.sq + topo + 
                           (1 + nitro + nitro.sq | gridId))
data.topo2 <- make_standata(form.topo2,
                           data=gridded_df_ordered)
data.topo2$nrow <- 93
data.topo2$ncol <- 18

nuts_controls   <- list(max_treedepth = 12, adapt_delta = 0.95)


stan_GS_fix_Topo <- stan("StanModel/GS_rand_faster_try2.stan",
                          data  = data.topo2,
                          control = nuts_controls,
                          cores = 4)
# saveRDS(stan_GS_fix_Topo,"stan_GS_fix_Topo.rds")


stan_ST_fix_Topo <- readRDS("Outcome/stan_GS_fix_Topo.rds")

check_hmc_diagnostics(stan_ST_fix_Topo)

traceplot(stan_ST_fix_Topo,c("b"))

post.gs.fix.topo <- extract(stan_ST_fix_Topo)

M <- length(post.gs.fix.topo$Intercept)

pp_check(gridded_df_ordered$yield,
         post.gs.fix.topo$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_GSRNS.pdf",width=11.96,height=8.27)


loo.gs.fix.topo  <- loo(post.gs.fix.topo$log_lik)
psis.gs.fix.topo <- psis(-post.gs.fix.topo$log_lik)

rsq_bayes(post.gs.fix.topo)


sumstan(post.gs.fix.topo$b)
sumstan(post.gs.fix.topo$rho_r)
sumstan(post.gs.fix.topo$rho_c)


output <- post.gs.rand.topo
output <- post.gs.fix.topo

N <- nrow(gridded_df_ordered)
beta0.topo <- matrix(0,M,N)
beta1.topo <- matrix(0,M,N)
beta2.topo <- matrix(0,M,N)
for(i in 1:M){
  beta0.topo[i,] <- output$b_Intercept[i] + output$r_1_1[i,]
  beta1.topo[i,] <- output$b[i,1] + output$r_1_2[i,]
  beta2.topo[i,] <- output$b[i,2] + output$r_1_3[i,]
  cat(paste0(round(i/M*100,2),"% "))
}

out_b0 <- apply(beta0.topo,2,median)
out_b1 <- apply(beta1.topo,2,median)
out_b2 <- apply(beta2.topo,2,median)



nitr <- numeric(N)
wh <- which(out_b2<0)
nitr[wh] <- pmin(-out_b1[wh]/(2*out_b2[wh]),124.6)
nitr[-wh] <- 124.6 

table(round(nitr,0))

temp <- data.frame(nitr = nitr, x = gridded_df_ordered$x,
                   y = gridded_df_ordered$y)

sp::coordinates(temp) <- c("x","y")

spplot(temp, 'nitr')

temp <- as(temp["nitr"], "ppp")
temp <- Smooth(temp)

plot(temp, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(temp,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)









