functions {
  /* compute the logm1 link 
   * Args: 
   *   p: a positive scalar
   * Returns: 
   *   a scalar in (-Inf, Inf)
   */ 
   real logm1(real y) { 
     return log(y - 1);
   }
  /* compute the inverse of the logm1 link 
   * Args: 
   *   y: a scalar in (-Inf, Inf)
   * Returns: 
   *   a positive scalar
   */ 
   real expp1(real y) { 
     return exp(y) + 1;
   }
   
  matrix chol_AR_matrix(real rho,int d){
    matrix[d,d] MatAR;
    MatAR = rep_matrix(0,d,d);
    for(i in 1:d){
      for(j in i:d){
        if(j>=i && i==1) MatAR[j,i] = rho^(j-1);
        else if(i>=2 && j>=i) MatAR[j,i] = rho^(j-i)*sqrt(1-rho^2);
      }
    }
    return MatAR;
  }
  
  matrix chol_kronecker_product(matrix matA, matrix matB) {
    matrix[rows(matA)*rows(matB),rows(matA)*rows(matB)] matC;
    matC = rep_matrix(0,rows(matA)*rows(matB),rows(matA)*rows(matB));
    for (k in 1:rows(matA)){
      for (l in 1:k){
        for (m in 1:rows(matB)){
          for (n in 1:m){
            matC[rows(matB)*(k-1)+m, cols(matB)*(l-1)+n] = matA[k,l] * matB[m,n];
          }
        }
      }
    }
    return matC;
  }
  
  matrix as_matrix(vector X, int N, int K) {
    matrix[N, K] Y;
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]);
    }
    return Y;
  }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  vector[N] Z_1_3;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
  
  int nrow;
  int ncol;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // residual SD
  real<lower=1> nu;  // degrees of freedom or shape
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1*ncol,nrow] z_1;  

  // vector[ncol] z_2;
  // vector[M_1]  z_3;
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  
  real<lower=0,upper=1> rho_r;
  real<lower=0,upper=1> rho_c;
}
transformed parameters {
  matrix[N_1,M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  vector[N_1] r_1_3;
  // compute actual group-level effects

  r_1 =  as_matrix(to_vector(chol_kronecker_product(chol_AR_matrix(rho_c,ncol),diag_pre_multiply(sd_1, L_1))* z_1 * chol_AR_matrix(rho_r,nrow)'),N_1,M_1);

  r_1_1 = r_1[,1];
  r_1_2 = r_1[,2];
  r_1_3 = r_1[,3];
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n];
  }
  // priors including all constants
  
  target += normal_lpdf(Intercept | 80, 10) + normal_lpdf(b[1] | 0,0.01)+normal_lpdf(b[2] | 0,0.001);
  target += normal_lpdf(sd_1[1] | 0, 1) + normal_lpdf(sd_1[2] | 0,0.01)+normal_lpdf(sd_1[3] | 0,0.001);
  target += normal_lpdf(sigma | 0,1);
  target += uniform_lpdf(rho_r | 0,1)+uniform_lpdf(rho_c | 0,1);
  target += gamma_lpdf(nu | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);

  target += std_normal_lpdf(to_vector(z_1));
  // target += std_normal_lpdf(to_vector(z_2));
  // target += std_normal_lpdf(to_vector(z_3));
  
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  // likelihood including all constants
  if (!prior_only) {
    target += student_t_lpdf(Y | nu, mu, sigma);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu = Intercept + Xc * b;
  
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n];
    log_lik[n] = student_t_lpdf(Y[n] | nu, mu[n], sigma);
    y_rep[n]   = student_t_rng(nu,mu[n],sigma);
  }
}
