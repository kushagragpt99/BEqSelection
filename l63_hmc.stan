//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
functions {
  
  vector make_tilde(vector X, real t){
    vector[12] X_vec = [1, X[1], X[2], X[3], X[1] ^ 2, X[2] ^ 2, X[3] ^ 2, X[1] * X[2], X[2] * X[3], X[3] * X[1], t, t ^ 2]';
    return(X_vec);
  } 
  
  real X_B_lpdf(vector state, int n_X, int n_theta, int N, int K, matrix y, int[] seq_t, matrix inv_R, matrix inv_lam_0, vector tau_0, real mu, real sigma2, real del_t, real a4, real b4) {
    matrix[3,N+1] X_n = to_matrix(state[1:n_X], 3, N+1);
    matrix[3,12] B = to_matrix(state[(n_X + 1):(n_X + n_theta)], 3, 12);
    matrix[3,K] X_t = X_n[1:3, seq_t];
    real p1 = 0;
    vector[3] y_k;
    vector[3] X_t_k;
    real p2 = 0;
    vector[3] p3_vec_tmp=rep_vector(0,3);
    real p3 = 0;
    vector[3] p3_vec=rep_vector(0,3);
    
    for(k in 1:K){
      y_k = y[1:3,k];
      X_t_k = X_t[1:3,k];
      p1 += (y_k - X_t_k)' * inv_R * (y_k - X_t_k);
    }
    p1 *= -0.5;
    p1 += (-0.5)*(X_n[1:3,1] - tau_0)' * inv_lam_0 * (X_n[1:3,1] - tau_0);
    p2 = (-0.5)*sum((B - mu).*(B - mu)) / sigma2;
    
    for(i in 1:N){
      p3_vec_tmp = (X_n[1:3,i+1] - X_n[1:3,i])/del_t - B * make_tilde(X_n[1:3,i], del_t*(i-1));
      p3_vec += (p3_vec_tmp).*(p3_vec_tmp);
    }
    p3_vec *= (del_t*0.5);
    p3 = 3 * lgamma(a4 + N*0.5) - (a4 + N*0.5) * sum(log(b4 + p3_vec));
    
    return(p1+p2+p3);
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[3,K+1] y;
  int seq_t[K];
  matrix[3,3] R;
  vector[3] tau_0;
  matrix[3,3] lam_0;
  real mu;
  real sigma2;
  real del_t;
  real a4;
  real b4;
  matrix[3,3] inv_R;
  matrix[3,3] inv_lam_0;
  int n_X;
  int n_theta;
  vector[n_X + n_theta] state;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix[3, N+1] X_n;
  matrix[3,12] B;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  //target+= X_B_lpdf(state| n_X, n_theta, N, K, y, seq_t, inv_R, inv_lam_0, tau_0, mu, sigma2, del_t, a4, b4);
    matrix[3,K] X_t = X_n[1:3, seq_t];
    real p1 = 0;
    vector[3] y_k;
    vector[3] X_t_k;
    real p2 = 0;
    vector[3] p3_vec_tmp=rep_vector(0,3);
    real p3 = 0;
    vector[3] p3_vec=rep_vector(0,3);
    
    for(k in 1:K){
      y_k = y[1:3,k];
      X_t_k = X_t[1:3,k];
      p1 += (y_k - X_t_k)' * inv_R * (y_k - X_t_k);
    }
    p1 *= -0.5;
    p1 += (-0.5)*(X_n[1:3,1] - tau_0)' * inv_lam_0 * (X_n[1:3,1] - tau_0);
    p2 = (-0.5)*sum((B - mu).*(B - mu)) / sigma2;
    
    for(i in 1:N){
      p3_vec_tmp = (X_n[1:3,i+1] - X_n[1:3,i])/del_t - B * make_tilde(X_n[1:3,i], del_t*(i-1));
      p3_vec += (p3_vec_tmp).*(p3_vec_tmp);
    }
    p3_vec *= (del_t*0.5);
    p3 = 3 * lgamma(a4 + N*0.5) - (a4 + N*0.5) * sum(log(b4 + p3_vec));
    target += p1+p2+p3;
}

