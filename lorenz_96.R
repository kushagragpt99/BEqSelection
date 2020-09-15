set.seed(1)
library(mvtnorm)
library(matrixcalc)
library(mcmc)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = matrix(, nrow = 40, ncol = 1)
    for (i in 0:39) {
        ans[i+1,1] = (X[(i+1)%%40 + 1] - X[(i-2)%%40 + 1]) * X[(i-1)%%40 + 1] - X[i + 1] + theta
    }
    return(ans)
}

make_Sigma <- function(sigma_ar, n) {
    Sigma = matrix(, nrow = n, ncol = n)
    k = 1
    start = 1
    for (i in 1:n) {
        len = n - i + 1
        if (i != 1)
            start = start + n - i + 2
        for (j in 1:len) {
            Sigma[i, (i - 1 + j)] = sigma_ar[start + j - 1]
            Sigma[(i - 1 + j), i] = sigma_ar[start + j - 1]
        }
    }
    return(Sigma)
}
#Sig = make_Sigma(sigma_init, 4)
#make_Sigma(1:10, 4)

# logarithm of the unnormalized posterior
ludfun <- function(state) {
    # State is the vector storing the vectors of length 40*N + 47. The first 40*(N+1) terms are Xs. The next term is the drift partameter \theta. 
    # The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:(40 * (N + 1))], nrow = 40, ncol = N + 1)
    Sigma_vec = state[(40 * N + 42):(40 * N + 41 + 820)]
    Sigma = make_Sigma(Sigma_vec, 40)
    theta = state[40*N + 41] # vector of \sigma, \rho and \beta

    # \Sigma should be positive semi-definite
    if (is.positive.semi.definite(Sigma) == FALSE) {
        #print("mehh")
        return(-Inf)
    }

    X_t = X_n[, seq(1, N + 1, N / K)]
    inv_R = solve(R)

    # pi is the log of likelihood
    p1 = 0
    for (k in 1:K) {
        Y.t = t(t(Y[, k]))
        X_t.t = t(t(X_t[, k]))
        p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    }
    p1 = -0.5 * p1

    # p2 is the log of prior of X conditional on theta
    p2 = 0

    inv_Sig = solve(Sigma)
    for (k in 1:N) {
        del_X = matrix(X_n[, k + 1] - X_n[, k], nrow = 40, ncol = 1)
        f_k = drift_fun(X[, k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + t(del_X / del_t - f_k) %*% inv_Sig %*% (del_X / del_t - f_k)
    }
    p2 = -0.5 * p2
    p2 = p2 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% solve(lam_o) %*% (t(t(X_n[, 1])) - tau_o) - (N / 2) * log(det(Sigma * del_t))

    # p3 is the lof og priors of theta
    p3 = (theta - mu)^2 / (2* sigma^2)

    return(p1 + p2 + p3)

}

# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = 40, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + det(del_t * Sigma) * rmvnorm(1, sigma = diag(1, 40))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

# hyper-parameters
to = 0 # initial time
tf = 10 # final time
Nobs = 8 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 40), nrow = 40, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(1, 40) # prior covariance matrix of X[0]
mu = 8
sigma = 2

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
R = diag(1, 40) # observational error

X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, 8, diag(2, 40)) # generating sample from Lorenz-63
Y = X[, seq(1, N + 1, N / K)] + t(rmvnorm(K + 1, mean = rep(0, 40), sigma = R)) # observations from Lorenz-63
init = runif(40 * N + 41 + 820, 0, 5) # random initial values for MCMC
sigma_init = diag(1, 40)
sigma_init[upper.tri(sigma_init)] <- NA
sigma_init = as.vector(sigma_init)
sigma_init = sigma_init[!is.na(sigma_init)]
init[(40 * N + 42):(40 * N + 41 + 820)] = sigma_init # inital \Sigma should also be positive semi definite

chain = metrop(ludfun, init, nbatch = 1e3, scale = 0.01) # running MH
