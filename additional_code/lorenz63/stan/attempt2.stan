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
  vector make_tilde(vector X, real t) {
    vector[11] Xtilde = [X[1], X[2], X[3], X[1]^2, X[2]^2, X[3]^2, X[1]*X[2], X[2]*X[3], X[3]*X[1], t, t^2]';
    return(Xtilde);
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> n_theta;
  int<lower=0> n_X;
  int seq_t[K];
  matrix[3,K] y;
  matrix[3,3] inv_R;
  matrix[3,3] inv_lam_0;
  vector[3] tau_0;
  vector[n_theta] mu;
  real sigma2;
  real del_t;
  real a4;
  real b4;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[n_X] X;
  vector[n_theta] B_vec;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[3, N+1] X_n = to_matrix(X, 3, N+1);
  matrix[3, K] X_t = X_n[1:3, seq_t];
  matrix[3,11] B = to_matrix(B_vec, 3, 11);
  vector[3] X_t_k;
  vector[3] y_k;
  real q1 = 0;
  vector[3] denom_vec = rep_vector(0,3);
  vector[3] del_X;
  vector[11] Xtilde_tmp;
  real denom = 0;
  
  for(k in 1:K){
    y_k = y[1:3,k];
    X_t_k = X_t[1:3,k];
    q1+= (y_k - X_t_k)' * inv_R * (y_k - X_t_k);
  }
  q1+= (X_n[1:3,1] - tau_0)' * inv_lam_0 * (X_n[1:3,1] - tau_0);
  q1+= dot_self(B_vec - mu) / sigma2;
  q1*= (-0.5);
  
  for(i in 1:N){
    del_X = X_n[1:3,i+1] - X_n[1:3,i];
    Xtilde_tmp = make_tilde(X_n[1:3,i], del_t*(i-1));
    denom_vec+= rows_dot_self(del_X/del_t - B*Xtilde_tmp);
  }
  //print(denom_vec);
  denom_vec*= (del_t * 0.5);
  denom_vec+= b4;
  denom = sum(log(denom_vec));
  denom*= (a4 + N*0.5);
  target+= q1 - denom;
}

