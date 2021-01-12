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
  vector drift_fun(vector X, vector theta){
    vector[3] ans;
    ans[1] = theta[1]*(X[2] - X[1]);
    ans[2] = theta[2]*X[1] - X[2] - X[1]*X[3];
    ans[3] = X[1]*X[2] - theta[3]*X[3];
    return(ans);
  }  
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int n_X;
  int n_theta;
  int n_sigma;
  matrix[3,K] y;
  int seq_t[K];
  matrix[3,3] inv_R;
  matrix[3,3] inv_lam_0;
  vector[3] tau_0;
  real del_t;
  real a1;
  real a2;
  real a3;
  real b1;
  real b2;
  real b3;
  real a4;
  real b4;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[n_X] X;
  vector[n_theta] theta;
  vector[n_sigma] sigma_vec;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[3,N+1] X_n = to_matrix(X, 3, N+1);
  matrix[3,3] Sigma = diag_matrix(sigma_vec); //to_matrix(sigma_vec, 3, 3);
  matrix[3,K] X_t = X_n[,seq_t];
  vector[3] y_k;
  vector[3] X_t_k;
  vector[3] del_X;
  vector[3] f;
  matrix[3,3] inv_Sigma;
  real p2=0;
  real p1=0;
  
  for(k in 1:K){
    y_k = y[1:3,k];
    X_t_k = X_t[1:3,k];
    p1+= (y_k - X_t_k)' * inv_R * (y_k - X_t_k);
  }
  p1+= (X_n[1:3,1] - tau_0)' * inv_lam_0 * (X_n[1:3,1] - tau_0);
  p1*= -0.5;
  
  inv_Sigma = inverse(Sigma);
  for(i in 1:N){
    del_X = X_n[1:3,i+1] - X_n[1:3,i];
    f = drift_fun(X_n[1:3,i], theta);
    p2+= (del_X/del_t - f)' * inv_Sigma * (del_X/del_t - f);
  }
  p2*= -0.5*del_t;
  
  p1+= (a1-1)*log(theta[1]) + (a2-1)*log(theta[2]) + (a3-1)*log(theta[3]) - theta[1]/b1 - theta[2]/b2 - theta[3]/b3;
  p1+= -(N*0.5)*log_determinant(Sigma);
  p2+= -(a4+1)*(log(Sigma[1,1]) + log(Sigma[2,2]) + log(Sigma[3,3])) - b4*(1.0/Sigma[1,1] + 1.0/Sigma[2,2] + 1.0/Sigma[3,3]);
  target+= p1+p2;
}

