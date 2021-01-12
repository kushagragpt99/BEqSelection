set.seed(1)

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
sigma2 = 10
mu_truth = c(rep(0, 3), -10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
n.X = 3 * (N + 1)
n.theta = 36
n.sigma = 3
n.param = n.X + n.theta + n.sigma
scale_iter = 1e3
n = 1e3

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
Y0 = Y[, - K]
Y1 = Y[, -1]
SVD_Y0 = svd(Y0)
A = t(SVD_Y0$u) %*% Y1 %*% SVD_Y0$v %*% diag((SVD_Y0$d) ^ -1)
eig = eigen(t(A) %*% A)
print(eig$values)

dynamic_modes = matrix(, nrow = 3, ncol = 3)

for (i in 1:3) {
    dynamic_modes[, i] = eig$values[i] * Y1 %*% SVD_Y0$v %*% diag((SVD_Y0$d) ^ -1) %*% eig$vectors[, i]
}

print(dynamic_modes)