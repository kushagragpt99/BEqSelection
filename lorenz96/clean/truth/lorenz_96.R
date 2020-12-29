set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = matrix(, nrow = N.l96, ncol = 1)
    for (i in 0:(N.l96 - 1)) {
        ans[i + 1, 1] = (X[(i + 1) %% N.l96 + 1] - X[(i - 2) %% N.l96 + 1]) * X[(i - 1) %% N.l96 + 1] - X[i + 1] + theta
    }
    return(ans)
}

ludfun <- function(state) {

    X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE)))

    p2 = dnorm(theta, mean = alpha, sd = beta, log = TRUE)

    f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale = rep(0.0012, n.X + n.theta) # 0.003
    scale[(n.X + 1):(n.X + n.theta)] = 0.27 # 0.5
    accept.prob = 0

    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

        chain = metrop(ludfun, init, 1, scale = scale)
        state = chain$batch
        accept.prob = accept.prob + chain$accept
        X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
        theta = state[(n.X + 1)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1] = theta

        Sigma = numeric(length = n.sigma)
        f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
        del_X = t(diff(t(X_n)))
        beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
        Sigma = sapply(1:n.sigma, function(t) rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[t]))

        param_mat[i, 2:(n.sigma + n.theta)] = Sigma
        init = state
    }

    print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg)
    return(final_output)
}


# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = N.l96, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

load('l96_5e4_cwise_spikes_truth_diffuse_try')

# hyper-parameters
to = attr$to # initial time
tf = attr$tf # final time
Nobs = attr$Nobs # no of observations (Y) per time step
N.l96 = attr$N.l96
del_t = attr$del_t # discrete approximation of dt
alpha = 8 # changed later
beta = 2
a4 = attr$a4
b4 = attr$b4


K = attr$K # no of real life observations, i.e. size of Y
N = attr$N # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = attr$burn_in  #/ del_t
R = attr$R # observational error
inv_R = attr$inv_R
seq.Y = attr$seq.Y 
n.X = attr$n.X
n.theta = 1
n.sigma = attr$n.sigma
n.param = attr$n.param
n = 1e6

alpha = mean(colMeans(to_save[[1]][[1]][5e3:1e4, 1:N.l96]))


X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.5, N.l96)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
#X = X[, 1:(N + 1)]

Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, N.l96), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[1:n.X] = as.numeric(X)
init[(n.X + 1):(n.X + n.theta)] = alpha

ans = linchpin(n, init)
pm = ans[[1]]
colMeans(pm)
#plot.ts(pm)
save(ans, file = "l96_linch_1e6")
