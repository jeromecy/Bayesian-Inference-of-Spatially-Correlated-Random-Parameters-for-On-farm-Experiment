source("Initialization.R")

summary(gridded_df_ordered)

skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

rsq_bayes <- function(post){
  var_mu <- apply(post$mu, 1, var)
  sigma2 <- post$sigma^2
  rsq_bayes <- var_mu / (var_mu + sigma2)
  print(summary(rsq_bayes))
  print(quantile(rsq_bayes,c(0.025,0.5,0.975)))
}

sumstan <- function(fit){
  
  return(c(mean(fit),
           sd(fit),
           quantile(fit,c(0.025,0.5,0.975))))
  
}

##################  GS NS #####################
color_scheme_set("blue")

# stan_GS_Rand_NS <- readRDS("Outcome/stan_GS_Rand_NS_pawsey1295.rds")
stan_GS_Rand_NS <- readRDS("Outcome/stan_GS_Rand_NS_1295_laptop.rds")

check_hmc_diagnostics(stan_GS_Rand_NS)

traceplot(stan_GS_Rand_NS,c("b"))

post.gs.rand.ns <- extract(stan_GS_Rand_NS)

M <- length(post.gs.rand.ns$Intercept)

pp_check(gridded_df_ordered$yield,
         post.gs.rand.ns$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_GSRNS.pdf",width=11.96,height=8.27)


loo.gs.rand.ns  <- loo(post.gs.rand.ns$log_lik)
psis.gs.rand.ns <- psis(-post.gs.rand.ns$log_lik)

rsq_bayes(post.gs.rand.ns)


################## GS Rand  ####################

# stan_GS_Rand <- readRDS("Outcome/stan_GS_Rand_pawsey1295.rds") ## faster code
# stan_GS_Rand <- readRDS("Outcome/stan_GS_SP_pawsey_try2.rds")
# stan_GS_Rand <- readRDS("Outcome/stan_GS_SP_try2_PC1095.rds")
# stan_GS_Rand <- readRDS("Outcome/stan_GS_SP_1599_try2_PC.rds")

stan_GS_Rand <- readRDS("Outcome/stan_GS_SP_1295_try2_PC.rds")
# stan_GS_Rand <- readRDS("Outcome/stan_GS_Rand_pawsey1599_rhos.rds")

check_hmc_diagnostics(stan_GS_Rand)

post.gs.rand <- extract(stan_GS_Rand)

M <- length(post.gs.rand$Intercept)

pp_check(gridded_df_ordered$yield,
         post.gs.rand$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_GSRand.pdf",width=11.96,height=8.27)


loo.gs.rand  <- loo(post.gs.rand$log_lik)
psis.gs.rand <- psis(-post.gs.rand$log_lik)

rsq_bayes(post.gs.rand)


################## ST NS  ####################

stan_ST_Rand_NS <- readRDS("Outcome/stan_ST_Rand_NS_pawsey1295.rds")

check_hmc_diagnostics(stan_ST_Rand_NS)

post.st.rand.ns <- extract(stan_ST_Rand_NS)

M <- length(post.st.rand.ns$Intercept)

pp_check(gridded_df_ordered$yield,
         post.st.rand.ns$y_rep[sample(1:M,100),],
         ppc_dens_overlay) +
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_STRNS.pdf",width=11.96,height=8.27)

loo.st.rand.ns  <- loo(post.st.rand.ns$log_lik)
psis.st.rand.ns <- psis(-post.st.rand.ns$log_lik)

rsq_bayes(post.st.rand.ns)


################## ST Rand  ####################

# stan_ST_Rand <- readRDS("Outcome/stan_ST_Rand_pawsey1295.rds")

# stan_ST_Rand <- readRDS("Outcome/stan_ST_SP_pawsey_try2.rds")
# stan_ST_Rand <- readRDS("Outcome/stan_ST_SP_try2_PC1095.rds")
# stan_ST_Rand <- readRDS("Outcome/stan_ST_SP_1599_try2_pawsey.rds")

stan_ST_Rand <- readRDS("Outcome/stan_ST_SP_1295_try2_PC.rds")
# stan_ST_Rand <- readRDS("Outcome/stan_ST_Rand_pawsey1599_rhos.rds")

check_hmc_diagnostics(stan_ST_Rand)

post.st.rand <- extract(stan_ST_Rand)

M <- length(post.st.rand$Intercept)

pp_check(gridded_df_ordered$yield,
         post.st.rand$y_rep[sample(1:M,100),],
         ppc_dens_overlay) +
  xlim(c(-20,160)) + thm1
# ggsave("SomeWriting/Images/ppcheck_STRand.pdf",width=11.96,height=8.27)

loo.st.rand  <- loo(post.st.rand$log_lik)
psis.st.rand <- psis(-post.st.rand$log_lik)

rsq_bayes(post.st.rand)


#### skew theme colour #####
color_scheme_set("blue")

ppc_stat(gridded_df_ordered$yield, 
         post.gs.rand.ns$y_rep,
         stat = "skew",binwidth=0.01) + thm1 +
  xlim(-0.5,0.3)
# ggsave("SomeWriting/Images/skew_GSRNS.pdf",width=11.96,height=8.27)


ppc_stat(gridded_df_ordered$yield, 
         post.gs.rand$y_rep,
         stat = "skew",binwidth=0.002) + thm1 +
  xlim(-0.5,0.3)
# ggsave("SomeWriting/Images/skew_GSRand.pdf",width=11.96,height=8.27)


ppc_stat(gridded_df_ordered$yield, 
         post.st.rand.ns$y_rep,
         stat = "skew",binwidth=0.01) + thm1 +
  xlim(-0.5,0.3)
# ggsave("SomeWriting/Images/skew_STRNS.pdf",width=11.96,height=8.27)


ppc_stat(gridded_df_ordered$yield, 
         post.st.rand$y_rep,
         stat = "skew",binwidth=0.01) + thm1 +
  # xlim(c(-0.5,0))
  xlim(-0.5,0.3)
# ggsave("SomeWriting/Images/skew_STRand.pdf",width=11.96,height=8.27)


#### PIT LOO theme colour #####
color_scheme_set("red")

ppc_loo_pit_overlay(gridded_df_ordered$yield,
                    yrep = post.gs.rand.ns$y_rep,
                    lw = weights(psis.gs.rand.ns))+
  thm1 +
  ylim(0,2.5)
# ggsave("SomeWriting/Images/pit_GSRNS_scale.pdf",width=11.96,height=8.27)


ppc_loo_pit_overlay(gridded_df_ordered$yield,
                    yrep = post.gs.rand$y_rep,
                    lw = weights(psis.gs.rand))+
  thm1 +
  ylim(0,2.5)
# ggsave("SomeWriting/Images/pit_GSRand_scale.pdf",width=11.96,height=8.27)


ppc_loo_pit_overlay(gridded_df_ordered$yield,
                    yrep = post.st.rand.ns$y_rep,
                    lw = weights(psis.st.rand.ns)) +
  thm1 +
  ylim(0,2.5)
# ggsave("SomeWriting/Images/pit_STRNS_scale.pdf",width=11.96,height=8.27)


ppc_loo_pit_overlay(gridded_df_ordered$yield,
                    yrep = post.st.rand$y_rep,
                    lw = weights(psis.st.rand)) +
  thm1 + 
  ylim(0,2.5)
# ggsave("SomeWriting/Images/pit_STRand_scale.pdf",width=11.96,height=8.27)



################ Results prd X ######################

sumstan(post.st.rand$Intercept)
sumstan(post.st.rand$b_Intercept)

sumstan(post.st.rand$b[,1])
sumstan(post.st.rand$b[,2])
sumstan(post.st.rand$sd_1[,1])
sumstan(post.st.rand$sd_1[,2])
sumstan(post.st.rand$sd_1[,3])
sumstan(post.st.rand$sigma)

# sumstan(post.st.rand$L_1[,2,1])
# sumstan(post.st.rand$L_1[,3,1])
# sumstan(post.st.rand$L_1[,3,2])

sumstan(post.st.rand$Cor_1[,1,2])
sumstan(post.st.rand$Cor_1[,1,3])
sumstan(post.st.rand$Cor_1[,2,3])


sumstan(post.st.rand$rho_c)
sumstan(post.st.rand$rho_r)
sumstan(post.st.rand$nu)


M <- length(post.st.rand$Intercept)
N <- nrow(gridded_df_ordered)

prdX.st <- as.data.frame(data.rand$X[,-1])
prdX.st$nitro <- 124.6 # 75.4 # 124.6
prdX.st$nitro.sq <- 124.6^2 # 75.4^2

y_prd.st.rand <- matrix(0,M,N)
y_st_prd <- matrix(0,M,N)

mus_prd.st.rand<- matrix(0,M,N)

for(i in 1:M){
  mus_prd.st.rand[i,] <- post.st.rand$b_Intercept[i] + 
    as.matrix(prdX.st)%*%post.st.rand$b[i,] +
    as.numeric(post.st.rand$r_1_1[i,]) + 
    as.numeric(post.st.rand$r_1_2[i,]*prdX.st[,1]) +
    as.numeric(post.st.rand$r_1_3[i,]*prdX.st[,2])
  
  y_st_prd[i,] <- rstudent_t(N,post.st.rand$nu[i],mus_prd.st.rand[i,],post.st.rand$sigma[i])
  cat(paste0(round(i/M*100,2),"% "));
}

prd.st     <- apply(y_st_prd,2,median)
prd.dat.st.rand <- data.frame(gridId = gridded_df_ordered$gridId)
prd.dat.st.rand$x <- as.data.frame(gridded_df_ordered$x)
prd.dat.st.rand$y <- as.data.frame(gridded_df_ordered$y)
prd.dat.st.rand$yield <- prd.st

sp::coordinates(prd.dat.st.rand) <- c("x","y")

spplot(prd.dat.st.rand, 'yield')

mod.yield.st.rand <- as(prd.dat.st.rand["yield"], "ppp")
yield.smooth.st.rand <- Smooth(mod.yield.st.rand)

#pdf("SomeWriting/Images/STRand_prd.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(yield.smooth.st.rand, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(yield.smooth.st.rand,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()



beta0 <- matrix(0,M,N)
beta1 <- matrix(0,M,N)
beta2 <- matrix(0,M,N)
for(i in 1:M){
  beta0[i,] <- post.st.rand$b_Intercept[i] + post.st.rand$r_1_1[i,]
  beta1[i,] <- post.st.rand$b[i,1] + post.st.rand$r_1_2[i,]
  beta2[i,] <- post.st.rand$b[i,2] + post.st.rand$r_1_3[i,]
  cat(paste0(round(i/M*100,2),"% "))
}

prd.dat.st.rand$b0 <- apply(beta0,2,median)
prd.dat.st.rand$b1 <- apply(beta1,2,median)
prd.dat.st.rand$b2 <- apply(beta2,2,median)

mod.b0 <- as(prd.dat.st.rand["b0"], "ppp")
b0.smooth <- Smooth(mod.b0)
#pdf("SomeWriting/Images/STRand_b0.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(b0.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b0.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()

mod.b1 <- as(prd.dat.st.rand["b1"], "ppp")
b1.smooth <- Smooth(mod.b1)
#pdf("SomeWriting/Images/STRand_b1.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(b1.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b1.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()

mod.b2 <- as(prd.dat.st.rand["b2"], "ppp")
b2.smooth <- Smooth(mod.b2)
#pdf("SomeWriting/Images/STRand_b2.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(b2.smooth, col=matlab.like2,  ribsep=0.01, 
     ribside="bottom", main="")
contour(b2.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()




############## optimal Nitrogen ###################

prd.dat.st.rand$N <- 0
wh <- which(prd.dat.st.rand$b2<0)
prd.dat.st.rand$N[wh] <- pmin(-prd.dat.st.rand$b1[wh]/(2*prd.dat.st.rand$b2[wh]),124.6)
prd.dat.st.rand$N[-wh] <- 124.6 


# optN <- matrix(0,M,N)
# for(i in 1:M){
#   beta0 <- post.st.rand$b_Intercept[i] + post.st.rand$r_1_1[i,]
#   beta1 <- post.st.rand$b[i,1] + post.st.rand$r_1_2[i,]
#   beta2 <- post.st.rand$b[i,2] + post.st.rand$r_1_3[i,]
# 
#   wh <- which(beta2<0)
#   optN[i,wh] <- pmin(-beta1[wh]/(2*beta2[wh]),124.6)
#   optN[i,-wh] <- 124.6
#   cat(paste0(round(i/M*100,2),"% "))
# }
# prd.dat.st.rand$N <- apply(optN,2,median)

mod.N <- as(prd.dat.st.rand["N"], "ppp")
N.smooth <- Smooth(mod.N)
#pdf("SomeWriting/Images/ST_opNitrogen.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(N.smooth, ribsep=0.01,col= heat.colors(100),  ## heat.colors, topo.colors
     ribside="bottom", main="")
contour(N.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()

# pdf("SomeWriting/Images/ST_opNitrogen_V2.pdf",width=11.96, height=4)
# par(mar = c(0,0,0,0))
# pal <- colorRampPalette(c("red","white"))
# plot(N.smooth, ribsep=0.01,col= pal(100),  ## heat.colors, topo.colors
#      ribside="bottom", main="")
# contour(N.smooth,"marks.nitro", add=TRUE, lwd=1, 
#         vfont=c("sans serif", "bold italic"), labcex=0.9)
# dev.off()

prd.dat.st.rand$prd <-prd.dat.st.rand$b0+prd.dat.st.rand$b1*prd.dat.st.rand$N + 
  prd.dat.st.rand$b2*prd.dat.st.rand$N^2

mod.prd <- as(prd.dat.st.rand["prd"], "ppp")
prd.smooth <- Smooth(mod.prd)
#pdf("SomeWriting/Images/ST_opYield.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(prd.smooth, col=terrain.colors(100,rev=TRUE),ribsep=0.01, 
     ribside="bottom", main="")
contour(prd.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()


sum(prd.dat.st.rand$prd-gridded_df_ordered$yield)



prd.dat.st.rand$diff <- prd.dat.st.rand$prd - prd.dat.st.rand$yield

# prd.dat.st.rand$diff <- prd.dat.st.rand$prd - gridded_df_ordered$yield

mod.diff <- as(prd.dat.st.rand["diff"], "ppp")
diff.smooth <- Smooth(mod.diff)

#pdf("SomeWriting/Images/ST_DiffYield.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(diff.smooth, col=terrain.colors(100,rev=TRUE),ribsep=0.01, 
     ribside="bottom", main="")
contour(diff.smooth,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()

# sum(prd.dat.st.rand$diff)
# 3273.801




sumstan(post.st.rand$rho_r)
sumstan(post.st.rand$rho_c)


sumstan(post.gs.rand$rho_r)
sumstan(post.gs.rand$rho_c)



sp::coordinates(gridded_df_ordered) <- c("x","y")

spplot(gridded_df_ordered, 'yield')

mod.gridded_df_ordered <- as(gridded_df_ordered["yield"], "ppp")
mod.gridded_df_ordered <- Smooth(mod.gridded_df_ordered)

#pdf("SomeWriting/Images/STRand_prd.pdf",width=11.96, height=4)
#par(mar = c(0,0,0,0))
plot(mod.gridded_df_ordered, col=terrain.colors(100,rev=TRUE),  ribsep=0.01, 
     ribside="bottom", main="")
contour(mod.gridded_df_ordered,"marks.nitro", add=TRUE, lwd=1, 
        vfont=c("sans serif", "bold italic"), labcex=0.9)
#dev.off()





############# corMat #################

corMat <- matrix(0,3,3)
corMat[1,] <- apply(post.st.rand$Cor_1[,1,], 2, median)
corMat[2,] <- apply(post.st.rand$Cor_1[,2,], 2, median)
corMat[3,] <- apply(post.st.rand$Cor_1[,3,], 2, median)
rownames(corMat) <- paste("sd",1:3)
colnames(corMat) <- c("$sigma[u[0]]","$sigma[u[1]]","$sigma[u[2]]")

corrplot(corMat)
corrplot.mixed(corMat,tl.cex=3)

sumstan(post.st.rand$Cor_1[,1,2])
sumstan(post.st.rand$Cor_1[,1,3])
sumstan(post.st.rand$Cor_1[,2,3])






############ AR1xAR1 residuals ############

color_scheme_set("blue")

stan_AR1 <- readRDS("Outcome/stan_AR1_pawsey_1599.rds")

check_hmc_diagnostics(stan_AR1)

traceplot(stan_AR1,c("b"))
traceplot(stan_AR1,c("rho_c","rho_r","sigma"))

post.ar1 <- extract(stan_AR1)

M <- length(post.ar1$Intercept)

pp_check(gridded_df_ordered$yield,
         post.ar1$y_rep[sample(1:M,100),],
         ppc_dens_overlay) + 
  xlim(c(-20,160)) + thm1

loo.ar1  <- loo(post.ar1$log_lik)
psis.ar1 <- psis(-post.ar1$log_lik)

rsq_bayes(post.ar1)


ppc_loo_pit_overlay(gridded_df_ordered$yield,
                    yrep = post.ar1$y_rep,
                    lw = weights(psis.ar1))












