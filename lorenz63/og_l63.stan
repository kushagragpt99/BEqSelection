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
  
  vector drift_fun(vector X, vector theta) {
    vector[3] ans;
    ans[1] = theta[1] * (X[2] - X[1]);
    ans[2] = theta[2] * X[1] - X[2] - X[1] * X[3];
    ans[3] = X[1] * X[2] - theta[3] * X[3];
    return(ans);
  }
  
}

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[3,K] y;
  int seq_t[K];
  matrix[3,3] R;
  vector[3] tau_0;
  matrix[3,3] lam_0;
  real del_t;
  real a4;
  real b4;
  matrix[3,3] inv_R;
  matrix[3,3] inv_lam_0;
  real alpha1;
  real alpha2;
  real alpha3;
  real beta1;
  real beta2;
  real beta3;
  int n_X;
  int n_theta;
  int n_sigma;
  int n_param;
  //vector[n_X + n_theta] state;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //matrix[3, N+1] X_n;
  //matrix[3,12] B;
  //vector[n_X+n_theta] state;
  vector[n_X] X;
  vector[n_theta] theta;
  vector[n_sigma] sigma_vec;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
    //target+= X_B_lpdf(state| n_X, n_theta, N, K, y, seq_t, inv_R, inv_lam_0, tau_0, mu, sigma2, del_t, a4, b4);
    //matrix[3,N+1] X_n = to_matrix(state[1:n_X], 3, N+1);
    //matrix[3,12] B = to_matrix(state[(n_X + 1):(n_X + n_theta)], 3, 12);

    matrix[3,N+1] X_n = to_matrix(X, 3, N+1);
    //vector[3] theta = state[(n_X + 1):(n.X + n_theta)];
    matrix[3,3] Sigma = diag_matrix(sigma_vec);

    matrix[3,K] X_t = X_n[1:3, seq_t];
    real p1 = 0;
    vector[3] y_k;
    vector[3] X_t_k;
    real p2 = 0;
    vector[3] del_X=rep_vector(0,3);
    real p3 = 0;
    vector[3] p3_vec=rep_vector(0,3);
    matrix[3,3] inv_Sig;
    real p4 = 0;
    vector[3] f_k;
    
    for(k in 1:K){
      y_k = y[1:3,k];
      X_t_k = X_t[1:3,k];
      p1 += (y_k - X_t_k)' * inv_R * (y_k - X_t_k);
    }
    p1 *= -0.5;
    //p1 += (-0.5)*(X_n[1:3,1] - tau_0)' * inv_lam_0 * (X_n[1:3,1] - tau_0);
    inv_Sig = inverse(Sigma);

    for(i in 1:N)   {
        del_X = (X_n[1:3,i+1] - X_n[1:3,i])/del_t;
        f_k = drift_fun(X_n[1:3,i], theta);
        p2 = p2 + (del_X / del_t - f_k)' * inv_Sig * (del_X / del_t - f_k);
    }
    p2 = p2 -0.5 * p2 * del_t;
    p2 = p2 - 0.5 * ( X_n[1:3, 1] - tau_0 )' * inv_lam_0 * ( X_n[1:3, 1] - tau_0 ) - (N *0.5) * log_determinant(Sigma*del_t);

    p3 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3;

    p4 = inv_gamma_lpdf(Sigma[1,1] | a4, b4) + inv_gamma_lpdf(Sigma[2,2] | a4, b4) + inv_gamma_lpdf(Sigma[3,3] | a4, b4);

    target += p1+p2+p3+p4;
}

