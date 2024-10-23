rm(list=ls())
library(dplyr)
library(ggplot2)
library(rgdal)
# library(sf)
library(raster)
library(lattice)
library(latticeExtra)
library(classInt)
library(GWmodel)
library(maptools)
library(spatstat)
library(tidyr)
library(agridat)
# library(asreml)
library(stringr)
library(mvtnorm)
library(matrixcalc)
library(rstan)
library(brms)
library(metR)
library(colorRamps)
library(corrplot)
library(ggridges)
library(mgcViz)
library(mgcv)
library(rethinking)
library(fields)
library(bayesplot)
library(reshape2)
library(SpatialTools)
library(loo)
library(MASS)
library(convoSPAT)
library(ggpubr)


MSE <- function(a,b) mean((a-b)^2)

thm1 <- theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.key.width = unit(1,"cm"),
        legend.position = "top",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size=20))

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
AR1CorMat <- function(rho,n) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}


# 
# chol_AR_matrix<- function(rho,d){
#   MatAR <- matrix(0,d,d)
#   for(i in 1:d){
#     for(j in i:d){
#       if(j>=i & i==1) MatAR[j,i] <- rho^(j-1)
#       else if(i>=2 & j>=i) MatAR[j,i] <- rho^(j-i)*sqrt(1-rho^2)
#     }
#   }
#   return (MatAR)
# }
# chol_kronecker_product<- function(matA, matB) {
#   matC <- matrix(0,nrow(matA)*nrow(matB),nrow(matA)*nrow(matB))
#   for (k in 1:nrow(matA)){
#     for (l in 1:k){
#       for (m in 1:nrow(matB)){
#         for (n in 1:m){
#           matC[nrow(matB)*(k-1)+m, nrow(matB)*(l-1)+n] = matA[k,l] * matB[m,n]
#         }
#       }
#     }
#   }
#   return(matC)
# }
# 
# AR_matrix<- function(rho,d){
#   MatAR = matrix(0,d,d)
#   for(i in 1:d){
#     for(j in i:d){
#       MatAR[i,j] = rho^abs(i-j)
#       MatAR[j,i] = MatAR[i,j]
#     }
#   }
#   return(MatAR)
# }
# 
# chol_kronecker_product_three_mat<- function(LA,LB,LC,d) {
#   
#   matD <- matrix(0,d,d)
#   for(iA in 1:nrow(LA)){
#     for(jA in 1:iA){
#       for(iB in 1:nrow(LB)){
#         for(jB in 1:iB){
#           for(iC in 1:nrow(LC)){
#             for(jC in 1:iC){
#               matD[nrow(LC)*(nrow(LB)*(iA-1)+iB-1)+iC,nrow(LC)*(nrow(LB)*(jA-1)+jB-1)+jC] = LA[iA,jA]*LB[iB,jB]*LC[iC,jC]
#             }
#           }
#         }
#       }
#     }
#   }
#   return(matD)
# }

# chol_kronecker_product_three<- function(LA,LB,LC,d) {
#   new_d = numeric(length(d))
#   for(iA in 1:nrow(LA)){
#     for(jA in 1:iA){
#       for(iB in 1:nrow(LB)){
#         for(jB in 1:iB){
#           for(iC in 1:nrow(LC)){
#             for(jC in 1:iC){
#               new_d[nrow(LC)*(nrow(LB)*(iA-1)+iB-1)+iC] = new_d[nrow(LC)*(nrow(LB)*(iA-1)+iB-1)+iC] + LA[iA,jA]*LB[iB,jB]*LC[iC,jC]*d[nrow(LC)*(nrow(LB)*(jA-1)+jB-1)+jC];
#             }
#           }
#         }
#       }
#     }
#   }
#   return(new_d)
# }


# ######## loading data ###########


gridded_df_ordered <- read.table("Data/gridded_df_ordered.txt", header = TRUE,
                                 stringsAsFactors = TRUE)


# gridded_df <- readRDS("Data/lasrosasGridded.rds")
# no_col <- 18
# no_row <- 93
# nel <- no_col*no_row
# gridded_df$gridobject <- factor(gridded_df$gridobject)
# gridded_df_ordered <- gridded_df[order(gridded_df$row, gridded_df$col), ]
# gridded_df_ordered$gridId <- factor(1:nrow(gridded_df))
# # Ctreating the spatial covariance matrix
# sigma_sq <- 1.0
# col_corr <- 0.95
# row_corr <- 0.70
# x <- 1:no_col #cols
# y <- 1:no_row #rows
# 
# d_col <- abs(outer(x, x, "-")) 
# sig_col <- col_corr^d_col
# d_row <- abs(outer(y, y, "-"))
# sig_row <- row_corr^d_row
# Cov_mat <- sigma_sq * direct.prod(sig_col, sig_row)
# dim(Cov_mat)
# attr(Cov_mat, "rowNames") <- levels(gridded_df_ordered$gridId)
# Cov_mat_inv <- solve(Cov_mat)
# attr(Cov_mat_inv, "rowNames") <- levels(gridded_df_ordered$gridId)
# 
# ######### stan fun ############
# 
# stanfun <- 
#   "
# functions{
#   matrix chol_AR_matrix(real rho,int d){
#     matrix[d,d] MatAR;
#     MatAR = rep_matrix(0,d,d);
#     for(i in 1:d){
#       for(j in i:d){
#         if(j>=i && i==1) MatAR[j,i] = rho^(j-1);
#         else if(i>=2 && j>=i) MatAR[j,i] = rho^(j-i)*sqrt(1-rho^2);
#       }
#     }
#     return MatAR;
#   }
#   
#   matrix chol_kronecker_product(matrix matA, matrix matB) {
#     matrix[rows(matA)*rows(matB),rows(matA)*rows(matB)] matC;
#     matC = rep_matrix(0,rows(matA)*rows(matB),rows(matA)*rows(matB));
#     for (k in 1:rows(matA)){
#       for (l in 1:k){
#         for (m in 1:rows(matB)){
#           for (n in 1:m){
#             matC[rows(matB)*(k-1)+m, cols(matB)*(l-1)+n] = matA[k,l] * matB[m,n];
#           }
#         }
#       }
#     }
#     return matC;
#   }
#   
#   matrix as_matrix(vector X, int N, int K) { 
#     matrix[N, K] Y; 
#     for (i in 1:N) {
#       Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
#     }
#     return Y; 
#   }
#   vector chol_kronecker_product_man(matrix LA, matrix LG, vector a) {
#     vector[num_elements(a)] new_a;
#     new_a = rep_vector(0, num_elements(a));
#     for(iA in 1:cols(LA)){
#       for(jA in 1:iA){
#         if(LA[iA, jA] > 1e-10){ // avoid calculating products between unrelated individuals
#           for(iG in 1:cols(LG)){
#             for(jG in 1:iG){
#               new_a[(cols(LG)*(iA-1))+iG] = new_a[(cols(LG)*(iA-1))+iG] + 
#                                             LA[iA, jA] * LG[iG, jG] * a[(cols(LG)*(jA-1))+jG];
#             }
#           }
#         }
#       }
#     }
#     return new_a;
#   }
#   
#   
#   vector chol_kronecker_product_three(matrix LA,matrix LB,matrix LC, vector d) {
#     vector[num_elements(d)] new_d;
#     new_d = rep_vector(0, num_elements(d));
#     for(iA in 1:cols(LA)){
#       for(jA in 1:iA){
#         for(iB in 1:cols(LB)){
#           for(jB in 1:iB){
#             for(iC in 1:cols(LC)){
#               for(jC in 1:iC){
#                 new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] = new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] + LA[iA,jA]*LB[iB,jB]*LC[iC,jC]*d[cols(LC)*(cols(LB)*(jA-1)+jB-1)+jC];
#               }
#             }
#           }
#         }
#       }
#     }
#     return new_d;
#   }
#   
#   vector chol_kronecker_two(matrix LA,matrix LB,vector d) {
#     vector[num_elements(d)] new_d;
#     new_d = rep_vector(0, num_elements(d));
#     for(iA in 1:cols(LA)){
#       for(jA in 1:iA){
#         for(iB in 1:cols(LB)){
#           for(jB in 1:iB){
#             new_d[cols(LB)*(iA-1)+iB] = new_d[cols(LB)*(iA-1)+iB] + LA[iA,jA]*LB[iB,jB]*d[cols(LB)*(jA-1)+jB];
#           }
#         }
#       }
#     }
#     return new_d;
#   }
#   
#   matrix kroVec(vector A, vector B){
#     vector[rows(A)*rows(B)] C;
#     for (j in 1:rows(A)) {
#       C[((j - 1)*rows(B) + 1):(j*rows(B))] = B*A[j];
#     }
#     return to_matrix(C,rows(A),rows(B),0);
#   }
#   
# }
# "

#expose_stan_functions(stanc(model_code = stanfun))
