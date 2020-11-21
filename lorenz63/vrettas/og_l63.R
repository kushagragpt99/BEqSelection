set.seed(1)
library(mvtnorm)
#library(matrixcalc)
library(mcmc)
library(invgamma)
library(rstan)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
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
R = diag(2, 3) # observational error
inv_R = solve(R)
n.X = 3 * (N + 1)
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma


X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.param)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- c(10, 28, 8 / 3) # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = 6 # inital \Sigma should also be positive semi definite

scale <- rep(.003, n.param)
scale[(n.X + 1):(n.X + n.theta)] <- .05
scale[(n.X + n.theta + 1):(n.param)] <- .2
#scale[c(6007, 6010, 6012)] <- 100

seq_t = seq(2, N + 1, N / K)
n = 2e4
burn_in_n = n/2

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63

options(mc.cores = 2)
init = numeric(n.param)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- c(10, 28, 8 / 3) # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = 6
initf <- function() {
    print('you shall not pass***************************************8')
    return(list(X = init[(1:n.X)], theta = init[(n.X + 1):(n.X + n.theta)], sigma_vec = init[(n.X + n.theta + 1):(n.param)]))
}
model = stan_model('attempt2og.stan')


#fit <- sampling(model, list(N = N, K = K, n_X = n.X, n_theta = n.theta, n_sigma = n.sigma, y = Y, seq_t = seq_t, inv_R = inv_R,
                #inv_lam_0 = inv.lam_o, tau_0 = tau_o[, 1], del_t = del_t, a1 = alpha1, a2 = alpha2, a3 = alpha3, b1 = beta1, b2 = beta2,
                #b3 = beta3, a4 = a4, b4 = b4), iter = n, warmup = burn_in_n, chains = 1, init = initf, algorithm = "HMC",
                #control = list(stepsize = 0.004, int_time = 0.2), pars = c("theta", "sigma_vec"))
#fit <- sampling(model, list(N = N, K = K, y = Y, seq_t = seq_t, R = R, tau_0 = tau_o[, 1], lam_0 = lam_o,
                #del_t = del_t, a4 = a4, b4 = b4, inv_R = inv_R, inv_lam_0 = inv.lam_o, n_X = n.X,
                #alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta1 = beta1, beta2 = beta2, beta3 = beta3,
                #n_theta = n.theta, n_sigma = n.sigma, n_param = n.param), iter = n, warmup = burn_in_n,
                #chains = 1, init = initf, algorithm = "HMC", control = list(stepsize = 0.04, int_time = 0.2),
                #pars = c("theta", "sigma_vec"))
fit <- sampling(model, list(N = N, K = K, n_X = n.X, n_theta = n.theta, n_sigma = n.sigma, y = Y, seq_t = seq_t, inv_R = inv_R,
                inv_lam_0 = inv.lam_o, tau_0 = tau_o[, 1], del_t = del_t, a1 = alpha1, a2 = alpha2, a3 = alpha3, b1 = beta1,
                b2 = beta2, b3 = beta3, a4 = a4, b4 = b4), iter = n, warmup = burn_in_n, chains = 1, init = initf,
                control = list(max_treedepth = 4), pars = c("theta", "sigma_vec"))

chain_info = capture.output(cat("no of samples from MC is ", n, " \n using warmup ", burn_in_n,
                 "max tree depth is ", 4, " \n starting from truth ", "\n priors centered at truth",
                 " time period ", tf))

print(chain_info)

to_save = list(fit, chain_info)
save(to_save, file = "nuts_vrettas_tf20_lam0_10") ######not 5, 3
p1 = extract(to_save[[1]], inc_warmup = TRUE, permuted = FALSE)
colMeans(p1[,1,])