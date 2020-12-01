set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

ludfun <- function(state, X2) {

    X1 = matrix(state[1:n.X1], nrow = 3, ncol = K)
    theta = state[(n.X1 + 1):(n.X1 + n.theta)] # vector of \sigma, \rho and \beta    


    # all the elements of theta should be positive
    if (min(theta) <= 0)
        return(-Inf)

    #p1 = 0
    ##print(dim(Y))
    #for (k in 1:K) {
        #Y.t = t(t(Y[, k]))
        #X1.t = t(t(X_t[, k]))
        #p1 = p1 + t(Y.t - X1.t) %*% inv_R %*% (Y.t - X1.t)
    #}
    #p1 = -0.5 * p1

    p1 = sum(dmvnorm(t(Y - X1), sigma = R, log = TRUE))

    X2.index = seq.X1 - 1
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, seq.X1] = X1
    X[, -seq.X1] = X2
    X2.relevant = X[, X2.index]

    p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    f = sapply(split(X2.relevant, rep(1:ncol(X2.relevant), each = nrow(X2.relevant))), drift_fun, theta)
    del_X = X1 - X2.relevant
    beta_tmp = rowSums((del_X / del_t - f) ^ 2) * del_t / 2
    p3 = -(a4 + K / 2) * sum(log(b4 + beta_tmp))

    ans = p1 + p2 + p3
    return(ans)
}

interpolate_X2 <- function(X1, theta) {
    X = matrix(, nrow = 3, ncol = N + 1)
    X[,1] = c(0,0,20)
    X[, seq.X1] = X1
    for (i in seq.X1) {
        for (j in 1:(N / K - 1)) {
            X[,i+j] = X[,i+j-1] + del_t*drift_fun(X[,i+j-1], theta)
        }
    }
    return(X)
}

approx_sampler <- function(n, init, init_X2) {

    param_mat = matrix(nrow = n, ncol = n.theta + n.sigma)
    X_avg = matrix(0, nrow = 3, ncol = N+1)
    scale = rep(0.01, n.X1 + n.theta)
    scale[(n.X1 + 1):(n.X1 + n.theta)] = 0.1
    accept.prob = 0
    X2 = init_X2

    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

        chain = metrop(ludfun, init, 1, scale = scale, X2 = X2)
        state = chain$batch
        accept.prob = accept.prob + chain$accept
        X_n = matrix(state[1:n.X1], nrow = 3, ncol = K)
        theta = state[(n.X1 + 1):(n.X1 + n.theta)] # vector of \sigma, \rho and \beta 
        param_mat[i, 1:n.theta] = theta

        Sigma = numeric(length = n.theta)
        X = matrix(, nrow = 3, ncol = N + 1)
        X[, seq.X1] = X_n
        X[, - seq.X1] = X2
        X2.relevant = X[, seq.X1 - 1]
        f = sapply(split(X2.relevant, rep(1:ncol(X2.relevant), each = nrow(X2.relevant))), drift_fun, theta)
        del_X = X_n - X2.relevant
        beta_tmp = rowSums((del_X / del_t - f) ^ 2) * del_t / 2
        Sigma[1] = rinvgamma(1, shape = K / 2 + a4, rate = b4 + beta_tmp[1])
        Sigma[2] = rinvgamma(1, shape = K / 2 + a4, rate = b4 + beta_tmp[2])
        Sigma[3] = rinvgamma(1, shape = K / 2 + a4, rate = b4 + beta_tmp[3])

        param_mat[i, (n.theta + 1):(n.theta + n.sigma)] = Sigma

        init = state
        X = interpolate_X2(X_n, theta)
        X2 = X[, - seq.X1]
        X_avg = X_avg + X
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
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


# hyper-parameters
to = 0 # initial time
tf = 20 # final time
Nobs = 50 # no of observations (Y) per time step
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


K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 / del_t
R = diag(1 / 2, 3) # observational error
inv_R = solve(R)
n.X = 3 * (N + 1)
n.X1 = 3*K
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma
seq.X1 = seq(2, N + 1, N / K)

n = 5e2*2

try_sigma = c(0.6, 6, 60)


for (r in 1:1) {
    print(r)
    b4 = try_sigma[r]
    #X_total = euler_maruyama(c(0, 0, 25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(try_sigma[r], 3)) # generating sample from Lorenz-63
    #X = X_total[, (burn_in):(N + burn_in)]
    #save(X, file = 'SA_X_6_by_10')
    load('SA_X_6_by_10')
    Y = X[, seq.X1] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
    init = numeric(n.X1 + n.theta)
    init[(1:n.X1)] <- as.numeric(X[, seq.X1]) + rnorm(n.X1, sd = 0.1) #runif(n.param, 0, 5)
    init[(n.X1+ 1):(n.X1 + n.theta)] <- rmvnorm(1, c(10, 28, 8 / 3), sigma = diag(0.5, 3)) # random initial values for MCMC

    X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
    #y.index = 1
    ##X.interp[,1] = X[,1]
    #for (i in seq(2, N + 1, N / K)) {
        #if (i == 2) {
            #X.interp[, 2] = Y[, 1]
        #} else {
            #X.interp[1, (i - N / K + 1):i] = seq(Y[1, y.index], Y[1, y.index + 1], (Y[1, y.index + 1] - Y[1, y.index]) * K / N)[-1]
            #X.interp[2, (i - N / K + 1):i] = seq(Y[2, y.index], Y[2, y.index + 1], (Y[2, y.index + 1] - Y[2, y.index]) * K / N)[-1]
            #X.interp[3, (i - N / K + 1):i] = seq(Y[3, y.index], Y[3, y.index + 1], (Y[3, y.index + 1] - Y[3, y.index]) * K / N)[-1]
            #y.index = y.index + 1
        #}

    #}

    init[(1:n.X1)] <- as.numeric(X.interp[, seq.X1])


    ans = approx_sampler(n, init, X.interp[,-seq.X1])
    pm = ans[[1]]
    print(colMeans(pm))
    plot.ts(pm)

    #chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from interpolate ", "\n priors spike slab ", " time period ",
                            #tf, " lam_0 is 1", "Nobs is ", Nobs, "Sigma is ", try_sigma[r], " R is 0.5"))

    #print(chain_info)
    #to_save = list(ans, chain_info)
    #save_name = paste('l63_linch_1e4_Sigma_x10_R_Nobs_SA', try_sigma[r] * 10, sep = '_')
    #save(to_save, file = save_name)
}