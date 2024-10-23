## data generation ##
datAll <- list()
datAll$N <- nrow(gridded_df_ordered)
datAll$Y <- gridded_df_ordered$yield
# datAll$X <- model.matrix(yield~1 + cnitro + cnitro.sq,
#                          data=gridded_df_ordered)
datAll$X <- model.matrix(yield~1 + nitro + nitro.sq,
                         data=gridded_df_ordered)
datAll$K <- 3 ### group levels rep1,rep2,rep3
datAll$J_1 <- as.numeric(gsub("R","",gridded_df_ordered$rep))   ## rep ID
datAll$J_2 <- gridded_df_ordered$gridId ## grid ID s

# datAll$Zs_1_1 <- data4$Zs_1_1
datAll$nrow <- 93
datAll$ncol <- 18


Xc <- matrix(0,datAll$N,2)
means_X <- numeric(2)
for (i in 2:3) {
  means_X[i - 1] = mean(datAll$X[, i]);
  Xc[, i - 1] = datAll$X[, i] - means_X[i - 1];
}
datAll$Xc <- Xc

















