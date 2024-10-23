// generated with brms 2.15.0
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
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
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
  vector[N_1] z_1[M_1];  // standardized group-level effects
  
  matrix[ncol, nrow] z_2;
  real<lower=0,upper=1> rho_r;
  real<lower=0,upper=1> rho_c;
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N] Sig;

  r_1_1 = (sd_1[1] * (z_1[1]));
  
  Sig = sigma * to_vector(chol_AR_matrix(rho_c,ncol) * z_2 * chol_AR_matrix(rho_r,nrow)');
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + Sig[n];
    }
    target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 84.7, 31.3);
  target += student_t_lpdf(sigma | 3, 0, 31.3)
    - 1 * student_t_lccdf(0 | 3, 0, 31.3);
  target += student_t_lpdf(sd_1 | 3, 0, 31.3)
    - 1 * student_t_lccdf(0 | 3, 0, 31.3);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(to_vector(z_2));
  target += uniform_lpdf(rho_r | 0,1) + uniform_lpdf(rho_c | 0,1);

}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu = Intercept + Xc * b;
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + Sig[n];
    log_lik[n] = normal_lpdf(Y[n] | mu[n], 1);
  }
}


