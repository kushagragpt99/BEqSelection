set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)
library(rstan)
#setwd("~/Python_Scripts/Bayesian_inference_computer_models")

make_tilde <- function(X,t) {
  X_vec = c(1, X[1], X[2], X[3], X[1] ^ 2, X[2] ^ 2, X[3] ^ 2, X[1] * X[2], X[2] * X[3], X[3] * X[1], t, t ^ 2)
  return(X_vec)
}
# drifet function for Lorenz-63
drift_fun <- function(X, t, B) {
  #print(make_tilde(X,t))
  tildeX = matrix(make_tilde(X, t), nrow = 12, ncol = 1)
  B_mat = matrix(B, nrow = 3)
  #print(B)
  #print(dim(tildeX))
  ans = B_mat %*% tildeX
  return(ans)
}

drift_fun_true <- function(X, theta) {
  ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
  return(t(t(ans)))
}

linchpin <- function(n, init) {
  X_avg = numeric(length = n.X)
  param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
  scale = rep(0.0005, n.X + n.theta)
  scale[(n.X + 1):(n.X + n.theta)] = 0.001
  scale[n.X + non_zero] = 0.01
  #scale[(n.X + 1):(n.X + 3) ] = 0.001
  scale[n.X + c(24,29)] = 0.008
  #scale[n.X + c(3,6,14,17,22,23)] = 0.003
  scale[n.X + 8] = 0.01
  scale[n.X + c(4,5,7)] = 0.08  # 0.05
  scale[n.X + c(7)] = 0.08
  #scale[n.X+c(3)] = 0.0008
  # scale[n.X+4] = 0.5
  scale[n.X+12] = 0.005
  accept.prob = 0
  #chain = metrop(ludfun, init, n, scale = scale)
  #print(chain$accept)
  for (i in 1:n) {
    if(i %% 1e3 == 0) print(c(i, accept.prob/i))
    chain = metrop(ludfun, init, 1, scale = scale)
    state = chain$batch
    accept.prob = accept.prob + chain$accept
    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
    X_avg = X_avg + state[1:n.X]
    param_mat[i, 1:n.theta] = theta
    
    Sigma = numeric(length = 3)
    f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(theta))
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    Sigma[1] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
    Sigma[2] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
    Sigma[3] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
    
    param_mat[i, (n.theta+1):(n.theta + n.sigma)] = Sigma
    init = state
  }
  print(accept.prob/n)
  X_avg = X_avg / n
  final_output = list(param_mat, X_avg)
  return(final_output)
}


# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
  X = matrix(, nrow = 3, ncol = N + 1)
  X[, 1] = X0
  for (i in 2:(N + 1))
    X[, i] = X[, i - 1] + t(drift_fun_true(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
  return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


# hyper-parameters
to = 0 # initial time
tf = 20 # final time
Nobs = 10 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha1 = 20 # Prior for \sigma is Gamma (alpha1, beta1)
alpha2 = 56 # Prior for \rho is Gamma (alpha2, beta2)
alpha3 = 6 # Prior for \beta is Gamma (alpha3, beta3)
beta1 = 0.5
beta2 = 0.5
beta3 = 0.5
a4 = 2
b4 = 6

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 / del_t
R = diag(2, 3) # observational error
inv_R = solve(R)
mu = 0
sigma2 = 1e1
mu_truth = c(rep(0, 3), -10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
mu = matrix(mu_truth, nrow = 3)
n.X = 3 * (N + 1)
n.theta = 36
n.sigma = 3
n.param = n.X + n.theta + n.sigma
seq_t = seq(2, N + 1, N / K)

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)

init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1,mu_truth,sigma=diag(1/50,n.theta))
non_zero = c(4,5,7,8,12,24,29)
named_list1 = list(Xn = X, B = matrix(init[(n.X + 1):(n.X + n.theta)], nrow = 3))
load("l63_linch_reg_bsv_0001_T_20_pv_10_init")
init[(n.X + 1):(n.X + n.theta)] <- head(tail(ans[[1]], 1)[1,], -3)
#ans = linchpin(1e4, init)
mu_truth = matrix(init[(n.X + 1):(n.X + n.theta)], nrow = 3)
model = stan_model('l63_hmc.stan')
named_list2 = list(Xn = X, B = matrix(init[(n.X + 1):(n.X + n.theta)], nrow = 3))
#initial = list(named_list1, named_list2)
initial = list(named_list2)
options(mc.cores = 2)

initf <- function() {
    print('you shall not pass***************************************8')
    return(list(Xn = X, B = matrix(init[(n.X + 1):(n.X + n.theta)], nrow = 3)))
}
mu = matrix(rep(0,36), nrow = 3)
fit <- sampling(model, list(N = N, K = K, y = Y, seq_t = seq_t, R = R, tau_0 = tau_o[,1], lam_0 = lam_o,
                            mu = mu, sigma2 = sigma2, del_t = del_t, a4 = a4, b4 = b4, inv_R = inv_R,
                            inv_lam_0 = inv.lam_o, n_X = n.X, n_theta = n.theta), iter = 1e5, chains = 1,
                            init = initf, control = list(max_treedepth = 3))

p1 = extract(fit, inc_warmup = TRUE, permuted = FALSE)
p2 = p1[, 1, (n.X + 1):(n.param - 3)]

save(fit, file = "L63_HMC_chain_1_mtd3_bp_pv_10")

## compare with n=1e5v metrop runs starting from truth - 7202.980 seconds
