addrew <- "C:/Users/279302d/OneDrive - Curtin/Documents/Curtin/Presentation/SAGI_Sym_2020"


ggplot(data=gridded_df_ordered,
       aes(nitro, yield)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = "gam",
              formula=y~s(x,k=3)) +
  theme_bw() + xlab("")+ylab("")

# ggsave(paste0(addrew,"/plots/LasrosasInadequateGlobal.pdf"),width=11.96,height=8.27)


ggplot(data=gridded_df_ordered,
       aes(nitro, yield, color = topo)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = "gam",
              formula=y~s(x,k=3)) +
  facet_wrap(~topo) +
  theme_bw() + xlab("")+ylab("") +
  theme(legend.position = "na")
# ggsave(paste0(addrew,"/plots/LasrosasInadequateLocal.pdf"),width=11.96,height=8.27)


ggplot(data=subset(gridded_df_ordered,row %in% (1:30)),
       aes(nitro, yield, color = row)) + 
  geom_point() + 
  geom_smooth(se = TRUE, method = "gam",
              formula=y~s(x,k=3)) +
  facet_wrap(~row) +
  theme_bw() + xlab("")+ylab("") +
  theme(legend.position = "na")
# ggsave(paste0(addrew,"/plots/LasrosasRows.pdf"),width=11.96,height=8.27)


##### RRC I ######


form.0 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId))

mod.0.update <- brm(form.0, 
                    data = gridded_df_ordered,
                    iter = 1000,
                    cores= 4)
color_scheme_set("blue")
pp_check(mod.0.update,nsamples=100) + thm1 + xlim(c(-10,160))
# ggsave(paste0(addrew,"/plots/pp_model1.pdf"),width=11.96,height=8.27)

loo.0 <-loo(mod.0.update, save_psis = TRUE)

yrep.0 <- posterior_predict(mod.0.update)  

ppc_loo_pit_overlay(
  y = gridded_df_ordered$yield,
  yrep = yrep.0,
  lw = weights(loo.0$psis_object)
) + thm1 
# ggsave(paste0(addrew,"/plots/pitloo_m0.pdf"),width=11.96,height=8.27)



bayes_R2(mod.0.update)


###### RRC AR1 ######
form.randU <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                            (1 + nitro + nitro.sq | gridId))

# prior check 
data.randU <- make_standata(form.randU,data=gridded_df_ordered)
data.randU$nrow <- 93
data.randU$ncol <- 18


## seeds  1072473474
# stanRandU <- readRDS("Outcome/StanRandUAR1_PC.rds")
# stanRandU <- readRDS("Outcome/StanRandUAR1_Pawsey.rds")
# stanRandU <- readRDS("Outcome/stan_randU_norep_cnitro_PC.rds")
# stanRandU <- readRDS("Outcome/stanRandU_norep_pawsey.rds")

# stanRandU@sim$warmup  <- 500
# stanRandU@sim$warmup2 <- c(500,500,500,500)
# stanRandU@stan_args[[1]]$warmup <- 500
# stanRandU@stan_args[[2]]$warmup <- 500
# stanRandU@stan_args[[3]]$warmup <- 500
# stanRandU@stan_args[[4]]$warmup <- 500
# 
# for(i in 1:4){
#   per1 <- stanRandU@sim$permutation[[i]]
#   wh <- order(per1)[c(1:100)]
#   stanRandU@sim$permutation[[i]] <- per1[-wh]-100
# }

# data.randU$rho_c <- 0.95
# data.randU$rho_r <- 0.95
# 
# stanRandU <- stan("StanModel/ar1Urand_NoRep_Prior_withRho.stan",
#                   data = data.randU,iter = 100,
#                   control = nuts_controls, 
#                   cores = 4)

stanRandU <- readRDS("Outcome/stanRandU_pawsey_warmup.rds")

check_hmc_diagnostics(stanRandU)
# get_num_leapfrog_per_iteration(stanRandU)

hmc_diagnostics <- nuts_params(stanRandU)

wh<- which(hmc_diagnostics$Parameter=="divergent__")
table(hmc_diagnostics[wh,]$Value)



post.rand <- as.array(rstan::extract(stanRandU))
names(post.rand)


summary(post.rand$b_Intercept)
summary(post.rand$b)
quantile(post.rand$b[,1],c(0.025,0.975))
quantile(post.rand$b[,2],c(0.025,0.975))


traceplot(stanRandU,c("rho_r","rho_c","b_Intercept"))

traceplot(stanRandU,c("sigma"))


pars <- c("b","b_Intercept","sigma","sd_1",
          "rho_r","rho_c")

N <- data.randU$N
i <- 1
j <- 1
M <- length(post.rand$Intercept)
y_rep_randU <- matrix(0,M,N)
loglik_randU <- matrix(0,M,N)

mus.rand <- matrix(0,M,N)
set.seed(128374)
for(i in 1:M){
  mus.rand[i,] <- post.rand$b_Intercept[i] + data.randU$X[,2:3]%*%post.rand$b[i,] +
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
         ppc_dens_overlay) + thm1 + xlim(c(-10,160))
# ggsave(paste0(addrew,"/plots/ppcheckRandU.pdf"),width=11.96,height=8.27)


pp_check(data.randU$Y, y_rep_randU[1:5,],  ppc_error_hist)
pp_check(data.randU$Y, y_rep_randU[1:5,],  ppc_hist) 

skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

ppc_stat(data.randU$Y, y_rep_randU[sample(1:M,1000),],
         stat = "skew", stat_bin=30) + thm1
# ggsave("Plots/RRC/ppcheckRandUHist.pdf",width=11.96,height=8.27)

div_style <- parcoord_style_np(div_color = "green", div_size = 0.15, div_alpha = 0.4)

mcmc_scatter(
  stanRandU,
  regex_pars = c("rho_r", "rho_c"),
  size = 3.5,
  alpha = 0.25,
  np = hmc_diagnostics,
  np_style = div_style
) + stat_density_2d(color = "black", size = .5) + 
  xlab(expression(rho[r]))+ ylab(expression(rho[c]))+
  thm1
# ggsave(paste0(addrew,"/plots/biscatter.pdf"),width=11.96,height=8.27)

mcmc_pairs(
  stanRandU,
  regex_pars = c("sd"),
  np = hmc_diagnostics,
  np_style = div_style,
  diag_fun = c("hist", "dens")
)

mcmc_nuts_divergence(hmc_diagnostics,log_posterior(stanRandU),
                     chain = 4)
# ggsave(paste0(addrew,"/plots/mcmc_nuts_divergence.pdf"),
#        width=11.96,height=8.27)

mcmc_nuts_energy(hmc_diagnostics)

# corMat <- matrix(0,3,3)
# corMat[1,] <- apply(post.rand$Cor_1[,1,], 2, mean)
# corMat[2,] <- apply(post.rand$Cor_1[,2,], 2, mean)
# corMat[3,] <- apply(post.rand$Cor_1[,3,], 2, mean)
# rownames(corMat) <- paste("sd",1:3)
# colnames(corMat) <- c("$sigma[1]","$sigma[2]","$sigma[3]")
# 
# corrplot(corMat)
# corrplot.mixed(corMat,tl.cex=3)


var_mu.rand <- apply(mus.rand, 1, var)
sigma2.rand <- post.rand$sigma^2
rsq_bayes.rand <- var_mu.rand / (var_mu.rand + sigma2.rand)

summary(rsq_bayes.rand)
quantile(rsq_bayes.rand,c(0.025,0.5,0.975))



# use loo for psis smoothing
loo.rand <- loo(post.rand$log_lik)
psis_result.rand <- psis(-post.rand$log_lik)
plot(psis_result.rand, label_points = TRUE)

color_scheme_set("red")
ppc_loo_pit_overlay(
  y = gridded_df_ordered$yield,
  yrep = y_rep_randU,
  lw = weights(psis_result.rand)
) + thm1
# ggsave(paste0(addrew,"/plots/pitloo.pdf"),width=11.96,height=8.27)



prdX <- matrix(c(75.4,75.4^2),N,2,byrow = TRUE)
prdX <- matrix(c(124.6,124.6^2),N,2,byrow = TRUE)

i <- 1
j <- 1
M <- length(post.rand$Intercept)
y_rep_randU_prd <- matrix(0,M,N)
ct <- 1
musprd <- matrix(0,M,N)

for(i in 1:M){
  musprd[i,] <- post.rand$b_Intercept[i] + 
    prdX%*%post.rand$b[i,] +
    as.numeric(post.rand$r_1_1[i,]) + 
    as.numeric(post.rand$r_1_2[i,]*prdX[,1]) +
    as.numeric(post.rand$r_1_3[i,]*prdX[,2])
  
  y_rep_randU_prd[i,] <- rnorm(N,musprd[i,],post.rand$sigma[i])
  cat(paste0(round(i/M*100,2),"% "));
}

prd     <- apply(y_rep_randU_prd,2,median)
prd.dat <- data.frame(gridId = gridded_df_ordered$gridId)
prd.dat$x <- as.data.frame(gridded_df_ordered$x)
prd.dat$y <- as.data.frame(gridded_df_ordered$y)
prd.dat$yield <- prd

sp::coordinates(prd.dat) <- c("x","y")

spplot(prd.dat, 'yield')

mod.yield <- as(prd.dat["yield"], "ppp")
yield.smooth <- Smooth(mod.yield)
# plot(yield.smooth, col=rev(rainbow(255)),  ribsep=0.01, 
#      ribside="bottom", main="")
# options(scipen = 999)

# pdf(paste0(addrew,"/plots/ModelRandU_prd_01.pdf"),
#     width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(yield.smooth, col=terrain.colors(n=100,rev=FALSE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(yield.smooth,"marks.nitro", add=TRUE, lwd=1, nlevels = 20,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()



beta0 <- matrix(0,M,N)
beta1 <- matrix(0,M,N)
beta2 <- matrix(0,M,N)
for(i in 1:M){
  beta0[i,] <- post.rand$b_Intercept[i] + post.rand$r_1_1[i,]
  beta1[i,] <- post.rand$b[i,1] + post.rand$r_1_2[i,]
  beta2[i,] <- post.rand$b[i,2] + post.rand$r_1_3[i,]
  cat(paste0(round(i/M*100,2),"% "))
}

prd.dat$b0 <- apply(beta0,2,median)
prd.dat$b1 <- apply(beta1,2,median)
prd.dat$b2 <- apply(beta2,2,median)

mod.b0 <- as(prd.dat["b0"], "ppp")
b0.smooth <- Smooth(mod.b0)
pdf(paste0(addrew,"/plots/ModelRandU_b0_01.pdf"),
    width=11.96, height=4)
par(mar = c(0,0,0,0))
plot(b0.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b0.smooth,"marks.nitro", add=TRUE, lwd=1, nlevels = 20,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
dev.off()



mod.b1 <- as(prd.dat["b1"], "ppp")
b1.smooth <- Smooth(mod.b1)
pdf(paste0(addrew,"/plots/ModelRandU_b1_01.pdf"),
    width=11.96, height=4)
par(mar = c(0,0,0,0))
plot(b1.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b1.smooth,"marks.nitro", add=TRUE, lwd=1, nlevels = 20,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
dev.off()


mod.b2 <- as(prd.dat["b2"], "ppp")
b2.smooth <- Smooth(mod.b2)
pdf(paste0(addrew,"/plots/ModelRandU_b2_01.pdf"),
    width=11.96, height=4)
par(mar = c(0,0,0,0))
plot(b2.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b2.smooth,"marks.nitro", add=TRUE, lwd=1, nlevels = 15,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
dev.off()






# loo1.full   <- loo(loglik.full)
# elpdi1.full <- loo1.full$pointwise[, "elpd_loo"]
# 
# loo1.rep <- loo(loglik_randU.rep)
# 
# loo_compare(loo1.rand,loo1.out1,loo1.full,loo1.rep)
# 
# corrplot(post.rand$Cor_1[7,,])
# 
# sumstan(post.rand$b_Intercept)
# sumstan(post.rand$b[,1])
# sumstan(post.rand$sigma)
# sumstan(post.rand$rho_r)
# sumstan(post.rand$rho_c)
# 
# sumstan(post.rand$sd_1[,1])
# 
# names(post.rand)
# 
# sumstan(post.rand$r_1_1)


# prd.dat$bi1 <- apply(post.rand$r_1_2,2,mean)
# prd.dat$bi2 <- apply(post.rand$r_1_3,2,mean)
# 
# 
# prd.dat$bi1 <- post.rand$r_1_2[5,]
# prd.dat$bi2 <- post.rand$r_1_3[20,]
# 
# mod.bi1 <- as(prd.dat["bi1"], "ppp")
# bi1.smooth <- Smooth(mod.bi1)
# plot(bi1.smooth, col=matlab.like2,  ribsep=0.01, 
#      ribside="bottom", main="")
# contour(bi1.smooth,"marks.nitro", add=TRUE, lwd=1, 
#         vfont=c("sans serif", "bold italic"), labcex=0.9)
# 
# mod.bi2 <- as(prd.dat["bi2"], "ppp")
# bi2.smooth <- Smooth(mod.bi2)
# plot(bi2.smooth, col=matlab.like2,  ribsep=0.01, 
#      ribside="bottom", main="")
# contour(bi2.smooth,"marks.nitro", add=TRUE, lwd=1, 
#         vfont=c("sans serif", "bold italic"), labcex=0.9)





prd.dat$nitro <- gridded_df_ordered$nitro
prd.dat$N <- 0
wh <- which(prd.dat$b2<0)
prd.dat$N[wh] <- pmin(-prd.dat$b1[wh]/(2*prd.dat$b2[wh]),124.6)

prd.dat$N[-wh] <- 124.6 


# fixef(mod.1.update)
# 
# brmpost <- as.data.frame(extract(mod.1.update$fit))
# head(brmpost)






mod.N <- as(prd.dat["N"], "ppp")
N.smooth <- Smooth(mod.N)
# pdf(paste0(addrew,"/plots/Optimal_N_01.pdf"),
#     width=11.96, height=4)
# par(mar = c(0,0,0,0))
plot(N.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(N.smooth,"marks.nitro", add=TRUE, lwd=1, nlevels=10,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
# dev.off()


prd.dat$prdYield <-prd.dat$b0+prd.dat$b1*prd.dat$N + 
  prd.dat$b2*prd.dat$N^2

mod.prdYield <- as(prd.dat["prdYield"], "ppp")
prdYieldsmooth <- Smooth(mod.prdYield)
#pdf(paste0(addrew,"/plots/Optimal_yield_01.pdf"),
#    width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(prdYieldsmooth, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(prdYieldsmooth,"marks.nitro", add=TRUE, lwd=1, nlevels=20,
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()

### colour terrain.colors

sum(prd.dat$prdYield - gridded_df_ordered$yield)
sum(prd.dat$prdYield - prd)


