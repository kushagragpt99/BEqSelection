set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

make_tilde <- function(X, t) {
    X_vec = c(X[1], X[2], X[3], X[1] ^ 2, X[2] ^ 2, X[3] ^ 2, X[1] * X[2], X[2] * X[3], X[3] * X[1], t, t ^ 2)
    return(X_vec)
}
# drifet function for Lorenz-63
drift_fun <- function(X, t, B) {
    #print(make_tilde(X,t))
    tildeX = matrix(make_tilde(X, t), nrow = 11, ncol = 1)
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

ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 3)

    # all the elements of theta should be positive
    #if (min(theta) <= 0)
    #return(-Inf)

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]


    # pi is the log of likelihood
    # This doesn't need a loop
    p1 = 0
    #print(dim(Y))
    for (k in 1:K) {
        Y.t = t(t(Y[, k]))
        X_t.t = t(t(X_t[, k]))
        p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    }
    p1 = -0.5 * p1
    p1 = p1 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o)

    #######################################################################
    #p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################

    p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = 3 * lgamma(a4 + N / 2) - (a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

linchpin <- function(n, init, scale_vec) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale = rep(0.0005 * 1.5, n.X + n.theta)
    scale[(n.X + 1):(n.X + n.theta)] = 0.001/1.05  # highly susceptibe to this
    scale[n.X + non_zero] = scale_vec[non_zero]
    scale[n.X + param_i] = 1.5 * scale_vec[param_i]
    #scale[n.X + non_zero] = 0.002
    ##scale[(n.X + 1):(n.X + 3) ] = 0.001
    #scale[n.X + c(24,29)] = 0.002
    ##scale[n.X + c(3,6,14,17,22,23)] = 0.003
    #scale[n.X + 8] = 0.002
    #scale[n.X + c(4,5,7)] = 0.002  # 0.05
    #scale[n.X + c(7)] = 0.002
    ##scale[n.X+c(3)] = 0.0008
    ## scale[n.X+4] = 0.5
    #scale[n.X+12] = 0.002
    accept.prob = 0
    #chain = metrop(ludfun, init, n, scale = scale)
    #print(chain$accept)
    for (i in 1:n) {
        if (i %% (n/10) == 0) print(c(i, accept.prob / i))
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

        param_mat[i, (n.theta + 1):(n.theta + n.sigma)] = Sigma
        init = state
    }
    print(accept.prob / n)
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
lam_o = diag(1, 3) # prior covariance matrix of X[0]
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
mu_truth = c(-10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
n.X = 3 * (N + 1)
n.theta = 33
n.sigma = 3
n.param = n.X + n.theta + n.sigma
n <- 5e5

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('../../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)

init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1, mu_truth, sigma = diag(1 / 50, n.theta))
non_zero = c(4, 5, 7, 8, 12, 24, 29) - 3
param_i = c(1, 2, 4, 9)
load("../../l63_linch_reg_bsv_0001_T_20_pv_10_init")
init[(n.X + 1):(n.X + n.theta)] <- head(tail(ans[[1]], 1)[1, - c(1, 2, 3)], -3)

load('../l63_linch_T_20_pv_10_5e5_2')
var1 = cov(to_save[[1]][[1]][, 1:33])
scale_vec = 0.04 * sqrt(diag(var1))
ans = linchpin(n, init, scale_vec)

chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from init ", "\n priors centered at 0 with varuance ",
                            sigma2, " time period ", tf, " lam_0 is 1"))

print(chain_info)
to_save = list(ans, chain_info)
save(to_save, file = "l63_linch_T_20_pv_10_5e5_2_tuned")
pm = ans[[1]]
print(matrix(colMeans(pm), nrow = 3))