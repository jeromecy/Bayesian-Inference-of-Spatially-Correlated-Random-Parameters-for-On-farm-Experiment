// generated with brms 2.13.0
functions {
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
  
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
      for (i in 1:N) {
        Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
      }
    return Y; 
  }

  vector chol_kronecker_product_three(matrix LA,matrix LB,matrix LC, vector d) {
    vector[num_elements(d)] new_d;
    new_d = rep_vector(0, num_elements(d));
    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        for(iB in 1:cols(LB)){
          for(jB in 1:iB){
            for(iC in 1:cols(LC)){
              for(jC in 1:iC){
                new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] = new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] + LA[iA,jA]*LB[iB,jB]*LC[iC,jC]*d[cols(LC)*(cols(LB)*(jA-1)+jB-1)+jC];
              }
            }
          }
        }
      }
    }
    return new_d;
  }
  
  vector chol_kronecker_product(matrix LA, matrix LG, vector a) {
    vector[num_elements(a)] new_a;
    new_a = rep_vector(0, num_elements(a));
    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        if(LA[iA, jA] > 1e-10){ // avoid calculating products between unrelated individuals
          for(iG in 1:cols(LG)){
            for(jG in 1:iG){
              new_a[(cols(LG)*(iA-1))+iG] = new_a[(cols(LG)*(iA-1))+iG] + 
                                            LA[iA, jA] * LG[iG, jG] * a[(cols(LG)*(jA-1))+jG];
            }
          }
        }
      }
    }
    return new_a;
  }
  
  matrix chol_kronecker_product_two(matrix matA, matrix matB) {
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
  
  matrix chol_three(matrix LA,matrix LB,matrix LC, vector d) {
    vector[num_elements(d)] new_d;
    new_d = rep_vector(0, num_elements(d));
    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        for(iB in 1:cols(LB)){
          for(jB in 1:iB){
            for(iC in 1:cols(LC)){
              for(jC in 1:iC){
                new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] = new_d[cols(LC)*(cols(LB)*(iA-1)+iB-1)+iC] + LA[iA,jA]*LB[iB,jB]*LC[iC,jC]*d[cols(LC)*(cols(LB)*(jA-1)+jB-1)+jC];
              }
            }
          }
        }
      }
    }
    return to_matrix(new_d,cols(LA)*cols(LB),cols(LC),0);
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
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[M_1*N_1] z_1;  // standardized group-level effects
  vector[N] z_2;
  vector[N] z_3;

  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix

  // real<lower=0,upper=1> rho_r;
  // real<lower=0,upper=1> rho_c;
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  matrix[N_1, M_1] r_2;

  // cholesky_factor_cov[nrow] Chol_Mat_r;
  // cholesky_factor_cov[ncol] Chol_Mat_c;
  // cholesky_factor_cov[3]    Chol_Mat_3;
  cholesky_factor_cov[N]    Chol_Mat_rc;

  // cholesky_factor_cov[3*ncol] Chol_Mat_two;
  // vector[N*3]  Chol_Mat_three;

  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  vector[N_1] r_1_3;

  vector[3] sd_2;
  sd_2[1] = 1;
  sd_2[2] = 0.5;
  sd_2[3] = 0.3;

  // Chol_Mat_r = chol_AR_matrix(rho_r,nrow);
  // Chol_Mat_c = chol_AR_matrix(rho_c,ncol);
  // Chol_Mat_3 = chol_AR_matrix(0.5,3);
  Chol_Mat_rc = chol_AR_matrix(0.5,N);

  // Chol_Mat_two = chol_kronecker_product_two(Chol_Mat_c,Chol_Mat_3);
  // Chol_Mat_three = chol_kronecker_product(Chol_Mat_rc,Chol_Mat_3,z_1);
  
  r_1 = as_matrix(z_1,N,3);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  r_1_3 = r_1[, 3];
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n];
  }
  // priors including all constants
  target += normal_lpdf(b[1] | 0,0.1);
  target += normal_lpdf(b[2] | 0,0.01);
  target += student_t_lpdf(Intercept | 3, 84.7, 5);
  target += student_t_lpdf(sigma | 3, 0, 10)- 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_1 | 3, 0, 5)  - 3 * student_t_lccdf(0 | 3, 0, 5);
  // target += uniform_lpdf(rho_r | 0,1) + uniform_lpdf(rho_c | 0,1);
  
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(Y | mu, sigma);
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
  }
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(Y[n] | mu[n], sigma);
    y_rep[n]   = normal_rng(mu[n],sigma);
  }
  
}