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


    # all the elements of theta should be positive
    #if (min(theta) <= 0)
    #return(-Inf)

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE)))
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################
    p2 = dnorm(theta, mean = alpha, sd = beta, log = TRUE)
    #p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale = rep(0.0022, n.X + n.theta) # 0.003
    scale[(n.X + 1):(n.X + n.theta)] = 0.38 # 0.5
    accept.prob = 0

    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

        #if (i == 5e4) {
        #scale[1:n.X] = 0.05
        #}

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
        #for (j in 1:n.sigma) {
            #Sigma[j] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[j])
        #}

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


# hyper-parameters
to = 0 # initial time
tf = 10 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, N.l96), nrow = N.l96, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, N.l96) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha = 8 # changed later
beta = 2
a4 = 2
b4 = .5


K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 #/ del_t
R = diag(.05, N.l96) # observational error
inv_R = diag(1 / (0.05), N.l96)
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
n.X = N.l96 * (N + 1)
n.theta = 1
n.sigma = N.l96
n.param = n.X + n.theta + n.sigma
n = 1e6

load('l96_5e5_cwise_spikes_interp_diffuse_init_theta_try')
alpha = mean(colMeans(to_save[[1]][[1]][5e3:1e4,1:N.l96]))


X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.5, N.l96)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
X = X[, 1:(N + 1)]
#load('../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, N.l96), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
#init[1:n.X] = ans[[3]]
init[(n.X + 1):(n.X + n.theta)] = alpha

# STARTING FROM TRUTH
#init[(1:n.X)] <- as.numeric(X) #+ rnorm(n.X) #runif(n.param, 0, 5)
#init[(n.X + 1):(n.X + n.theta)] <- rnorm(1, alpha, 0.5) # random initial values for MCMC

X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
y.index = 1
#X.interp[,1] = X[,1]
for (i in seq(2, N + 1, N / K)) {
    if (i == 2) {
        X.interp[, 2] = Y[, 1]
    } else {
        for (j in 1:N.l96) {
            X.interp[j, (i - N / K + 1):i] = seq(Y[j, y.index], Y[j, y.index + 1], (Y[j, y.index + 1] - Y[j, y.index]) * K / N)[-1]
        }
        y.index = y.index + 1
    }

}

init[(1:n.X)] <- as.numeric(X.interp)


ans = linchpin(n, init)
pm = ans[[1]]
colMeans(pm)
#plot.ts(pm)
save(ans, file = "l96_linch_1e4_interp")
