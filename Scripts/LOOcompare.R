source("Initialization.R")
summary(gridded_df_ordered)

thm2 <-  xlim(c(-50,250))

loo_pit <-function(model,loo1){
  nms <- substitute(model)
  yrep.model <- posterior_predict(model)
  ppc_loo_pit_overlay(gridded_df_ordered$yield,
                      yrep = yrep.model,
                      lw = weights(loo1$psis_object)) + 
    thm1 + ggtitle(nms)
  ggsave(paste0("Plots/LOO/loopit_",nms,".png"),
         width=11.96,height=8.27)
}

pp_plot <-function(model){
  nms <- substitute(model)
  pp_check(model,nsamples=50) + thm2 + ggtitle(nms)
  ggsave(paste0("Plots/LOO/ppcheck_",nms,".png"),
         width=11.96,height=8.27)
}

cont_prd <-function(model,forDat){
  
  # sp::coordinates(forDat) <- c("x","y")
  nms <- deparse(substitute(model))
  est <- predict(model,forDat)
  forDat$prd <- est[,1]
  mod.prd    <- as(forDat["prd"], "ppp")
  prd.smooth <- Smooth(mod.prd)
  
  pdf(paste0("Plots/LOO/estyield_",nms,".pdf"),width=11.96, height=4)
  par(mar = c(0,0,1,0))
  plot(prd.smooth, col=terrain.colors(n=100,rev=TRUE), 
       ribsep=0.01, ribside="bottom", main=nms)
  contour(prd.smooth,"marks.nitro", add=TRUE, lwd=1, 
          vfont=c("sans serif", "bold italic"), labcex=0.9)
  dev.off()
  
}

####### models comparison --- computation ########
## form 0: 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId)

form.0 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId))

mod.0.gs <- brm(form.0, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.0.st <- brm(form.0, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)


## form 1: topo + 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId) + (1|row+col) 

form.1 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId)+
                        (1|row+col))

mod.1.gs <- brm(form.1, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.1.st <- brm(form.1,data = gridded_df_ordered,
                iter = 1000, cores = 4, 
                family = student)


## form 2: topo + 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId)

form.2 <- brmsformula(yield ~ topo + 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId))

mod.2.gs <- brm(form.2, data = gridded_df_ordered,
                iter = 1000, cores= 4)


mod.2.st <- brm(form.2, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)


## form 3: 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId) + (1|topo)

form.3 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId) +
                        (1|topo))

mod.3.gs <- brm(form.3, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.3.st <- brm(form.3, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)


## form 4: 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId) + (1|row)

form.4 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|gridId)+
                        (1|row))

mod.4.gs <- brm(form.4, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.4.st <- brm(form.4, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)


## form 5: 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|row+col)

form.5 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq|row+col))

mod.5.gs <- brm(form.5, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.5.st <- brm(form.5, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)


## form 6: topo + 1 + nitro + nitro.sq 
## + (1 + nitro + nitro.sq|gridId)

form.6 <- brmsformula(yield ~ topo +
                        (1 + nitro + nitro.sq|gridId))

mod.6.gs <- brm(form.6, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.6.st <- brm(form.6, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)

## form 7;  1 + nitro + nitro.sq 
## + (1 | row + col) 

form.7 <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 | row + col))

mod.7.gs <- brm(form.7, data = gridded_df_ordered,
                iter = 1000, cores= 4)

mod.7.st <- brm(form.7, data = gridded_df_ordered,
                iter = 1000, cores= 4,
                family = student)

###### model checking and evaluation ######

### pp check, loo

loo.0.gs <- loo(mod.0.gs, save_psis = TRUE)
loo.0.st <- loo(mod.0.st, save_psis = TRUE)
loo.1.gs <- loo(mod.1.gs, save_psis = TRUE)
loo.1.st <- loo(mod.1.st, save_psis = TRUE)
loo.2.gs <- loo(mod.2.gs, save_psis = TRUE)
loo.2.st <- loo(mod.2.st, save_psis = TRUE)
loo.3.gs <- loo(mod.3.gs, save_psis = TRUE)
loo.3.st <- loo(mod.3.st, save_psis = TRUE)
loo.4.gs <- loo(mod.4.gs, save_psis = TRUE)
loo.4.st <- loo(mod.4.st, save_psis = TRUE)
loo.5.gs <- loo(mod.5.gs, save_psis = TRUE)
loo.5.st <- loo(mod.5.st, save_psis = TRUE)
loo.6.gs <- loo(mod.6.gs, save_psis = TRUE)
loo.6.st <- loo(mod.6.st, save_psis = TRUE)
loo.7.gs <- loo(mod.7.gs, save_psis = TRUE)
loo.7.st <- loo(mod.7.st, save_psis = TRUE)


loo_compare(loo.0.gs,loo.0.st,loo.1.gs,loo.1.st,
            loo.2.gs,loo.2.st,loo.3.gs,loo.3.st,
            loo.4.gs,loo.4.st,loo.5.gs,loo.5.st,
            loo.6.gs,loo.6.st,loo.7.gs,loo.7.st)

bayes_R2(mod.0.gs)
bayes_R2(mod.0.st)
bayes_R2(mod.1.gs)
bayes_R2(mod.1.st)
bayes_R2(mod.2.gs)
bayes_R2(mod.2.st)
bayes_R2(mod.3.gs)
bayes_R2(mod.3.st)
bayes_R2(mod.4.gs)
bayes_R2(mod.4.st)
bayes_R2(mod.5.gs)
bayes_R2(mod.5.st)
bayes_R2(mod.6.gs)
bayes_R2(mod.6.st)
bayes_R2(mod.7.gs)
bayes_R2(mod.7.st)



pp_plot(mod.0.gs)
pp_plot(mod.0.st)
pp_plot(mod.1.gs)
pp_plot(mod.1.st)
pp_plot(mod.2.gs)
pp_plot(mod.2.st)
pp_plot(mod.3.gs)
pp_plot(mod.3.st)
pp_plot(mod.4.gs)
pp_plot(mod.4.st)
pp_plot(mod.5.gs)
pp_plot(mod.5.st)
pp_plot(mod.6.gs)
pp_plot(mod.6.st)
pp_plot(mod.7.gs)
pp_plot(mod.7.st)


### loo_pit check

loo_pit(mod.0.gs,loo.0.gs)
loo_pit(mod.0.st,loo.0.st)
loo_pit(mod.1.gs,loo.1.gs)
loo_pit(mod.1.st,loo.1.st)
loo_pit(mod.2.gs,loo.2.gs)
loo_pit(mod.2.st,loo.2.st)
loo_pit(mod.3.gs,loo.3.gs)
loo_pit(mod.3.st,loo.3.st)
loo_pit(mod.4.gs,loo.4.gs)
loo_pit(mod.4.st,loo.4.st)
loo_pit(mod.5.gs,loo.5.gs)
loo_pit(mod.5.st,loo.5.st)
loo_pit(mod.6.gs,loo.6.gs)
loo_pit(mod.6.st,loo.6.st)
loo_pit(mod.7.gs,loo.7.gs)
loo_pit(mod.7.st,loo.7.st)




### predict yield 

forEst <- data.frame(gridId = gridded_df_ordered$gridId,
                     x = gridded_df_ordered$x,
                     y = gridded_df_ordered$y,
                     row = gridded_df_ordered$row,
                     col = gridded_df_ordered$col,
                     topo = gridded_df_ordered$topo)

forEst$nitro <- 74.5
forEst$nitro.sq <- 74.5^2
sp::coordinates(forEst) <- c("x","y")


cont_prd(mod.0.gs,forEst)
cont_prd(mod.0.st,forEst)
cont_prd(mod.1.gs,forEst)
cont_prd(mod.1.st,forEst)
cont_prd(mod.2.gs,forEst)
cont_prd(mod.2.st,forEst)
cont_prd(mod.3.gs,forEst)
cont_prd(mod.3.st,forEst)
cont_prd(mod.4.gs,forEst)
cont_prd(mod.4.st,forEst)
cont_prd(mod.5.gs,forEst)
cont_prd(mod.5.st,forEst)
cont_prd(mod.6.gs,forEst)
cont_prd(mod.6.st,forEst)
cont_prd(mod.7.gs,forEst)
cont_prd(mod.7.st,forEst)


############


form.xy <- brmsformula(yield ~ 1 + nitro + nitro.sq + 
                        (1 + nitro + nitro.sq | gridId) + 
                        s(x,y))

mod.xy <- brm(form.xy, data = gridded_df_ordered,
              cores = 4)

pp_plot(mod.xy)


loo.xy <- loo(mod.xy, save_psis = TRUE)
loo_pit(mod.xy,loo.xy)


cont_prd(mod.xy,forEst)












