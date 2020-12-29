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

ludfun <- function(state, gamma) {


    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 3)


    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]


    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))) + dmvnorm(t(X_n[, 1] - tau_o), sigma = lam_o, log = TRUE)
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################
    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

ludfun.X <- function(state, gamma, all) {

    all[1:n.X] = state
    X_n = matrix(all[1:n.X], nrow = 3, ncol = N + 1)
    B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 3)


    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))) + dmvnorm(t(X_n[, 1] - tau_o), sigma = lam_o, log = TRUE)
    ######################################################################
    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

sample_gamma <- function(B_vec) {
    gamma = numeric(length = n.theta)
    for (i in 1:n.theta) {
        prob = q[i] * dnorm(B_vec[i], sd = tau1) / (q[i] * dnorm(B_vec[i], sd = tau1) + (1 - q[i]) * dnorm(B_vec[i], sd = tau0))
        gamma[i] = rbinom(1, 1, prob)
    }
    return(gamma)
}

MH.X <- function(init, n, scale, gamma, B_vec) {
    chain = matrix(, nrow = n, ncol = n.X)
    accept.prob = 0
    for (i in 1:n) {
        prop = sapply(init, function(t) rnorm(1, t, scale))
        prop_ludf = c(prop, B_vec)
        init_ludf = c(init, B_vec)
        if (log(runif(1)) < (ludfun(prop_ludf, gamma) - ludfun(init_ludf, gamma))) {
            init = prop
            accept.prob = accept.prob + 1
        }
        chain[i,] = init
    }
    ans = list(chain, accept.prob / n)
    return(ans)
}

MH.B <- function(index, init, n, scale, gamma, state) {
    chain = numeric(length = n)
    accept.prob = 0
    prop_ludf = state
    init_ludf = state
    for (i in 1:n) {
        prop = rnorm(1, init, scale)
        prop_ludf[n.X + index] = prop
        init_ludf[n.X + index] = init
        if (log(runif(1)) < (ludfun(prop_ludf, gamma) - ludfun(init_ludf, gamma))) {
            init = prop
            accept.prob = accept.prob + 1
        }
        chain[i] = init
    }
    ans = list(chain, accept.prob / n)
    return(ans)
}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = 2 * n.theta + n.sigma)
    
    scale.B = scale

    accept.prob = numeric(1 + n.theta)
    state = init

    for (i in 1:n) {
        gamma = sample_gamma(init[(n.X + 1):(n.X + n.theta)])
        param_mat[i, (n.theta + n.sigma + 1):(2 * n.theta + n.sigma)] = gamma

        if (i %% (n / 10) == 0) {
            print(i)
            print(matrix(accept.prob[2:(n.theta + 1)] / i, nrow = 3))
            #print(c(i, accept.prob / i))
        }

        all = init
        chain = metrop(ludfun.X, initial = init[1:n.X], nbatch = 1, scale = scale.X, gamma = gamma, all = all)
        accept.prob[1] = accept.prob[1] + chain$accept
        state[1:n.X] = chain$batch


        for (j in 1:n.theta) {

            ans = MH.B(j, init[n.X + j], 1, scale.B[j], gamma, state)
            accept.prob[j + 1] = accept.prob[j + 1] + ans[[2]]
            state[n.X + j] = ans[[1]]
        }

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
    #print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg, accept.prob / n)
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
Nobs = 20 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)

a4 = 2
b4 = .6
tau1 = 5
tau0 = 0.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
burn_in = 5000 / del_t
R = diag(1 / 20, 3) # observational error
inv_R = diag(20,3)
mu = 0
sigma2 = 10
mu_truth = c(-10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
non_zero = c(4, 5, 7, 8, 12, 24, 29) - 3
param_i = c(1, 2, 4, 9)
n.X = 3 * (N + 1)
n.theta = 33
n.sigma = 3
n.param = n.X + n.theta + n.sigma
q = rep(0.1, n.theta) #runif(n.theta)
q[non_zero] = 0.9
n <- 5e4

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(.6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
#save(X, file = '../burninX6_by_10')
load('../burninX6_by_10')
X = X[, 1:(N + 1)]
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)

init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)

init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1,mu_truth, sigma = diag(1 / 50, n.theta))

#load('l63_linch_T_20_2e3_cwise_1_spikes')
#ans = to_save[[1]]
#pm = ans[[1]][, 1:(n.sigma + n.theta)]
#init[(n.X + 1):(n.X + n.theta)] = colMeans(pm[1e3:2e3,1:n.theta])


#sigma_Y = mean(diag(var(t(Y))))
#tau0 = sqrt(sigma_Y / (10 * K))/10
#tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K)))*2

load('../l63_linch_T_20_5e5_1')
var1 = cov(to_save[[1]][[1]][, 1:33])
scale_vec = 2.7 * sqrt(diag(var1))

scale = 0.6 * scale_vec
scale[non_zero[c(1, 7)]] = 0.3 * scale_vec[non_zero[c(1, 7)]]
scale[non_zero[c(3)]] = 0.7 * scale_vec[non_zero[c(3)]]
scale[non_zero[c(2, 4)]] = 0.4 * scale_vec[non_zero[c(2, 4)]]
scale[non_zero[c(5)]] = 0.08 * scale_vec[non_zero[c(5)]]
scale[non_zero[c(6)]] = 0.07 * scale_vec[non_zero[c(6)]]
scale[c(3)] = 0.4 * scale_vec[c(3)]
scale[c(5)] = 0.7 * scale_vec[c(5)]
scale[c(10)] = 0.1 * scale_vec[c(10)]
scale[c(27, 33)] = 0.3 * scale_vec[c(27, 23)]
scale[c(11, 15, 16, 17, 19)] = 0.11 * scale_vec[c(11, 15, 16, 17, 19)]
scale[c(6, 7, 8, 13, 14, 24)] = 0.18 * scale_vec[c(6, 7, 8, 13, 14, 24)]
scale[c(12)] = 0.07 * scale_vec[c(12)]
scale[c(31, 32)] = 0.2 * scale_vec[c(31, 32)]
scale[c(22, 23, 25, 26)] = 0.3 * scale_vec[c(22, 23, 25, 26)]
scale[c(18)] = 0.05 * scale_vec[c(18)]
scale[c(28, 29)] = 0.18 * scale_vec[c(28, 29)]
scale[c(20)] = 0.12 * scale_vec[c(20)]

scale.X = 0.0001


ans = linchpin(n, init)

chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from previous run previous run", "\n priors spike slab ", " time period ",
                            tf, " lam_0 is 10"))

print(chain_info)

attr = list('to' = to, 'tf' = tf, 'Nobs' = Nobs, 'del_t' = del_t, 'a4' = a4, 'b4' = b4, 'tau0' = tau0, 'tau1' = tau1, 'tau_o' = tau_o, 'lam_o' = lam_o,
            'K' = K, 'N' = N, 'burn_in' = burn_in, 'R' = R, 'inv_R' = inv_R, 'mu_truth' = mu_truth, 'non_zero' = non_zero, 'param_i' = param_i,
            'n.X' = n.X, 'n.theta' = n.theta, 'n.sigma' = n.sigma, 'n.param' = n.param, 'q' = q, 'seq.Y' = seq.Y, 'n' = n, 'scale_vec' = scale_vec,
            'scale' = scale, 'scale.X' = scale.X)


to_save = list(ans, chain_info)
save(to_save,attr, file = "l63_linch_T_20_5e4_cwise_1_spikes_init_theta_try")
pm = ans[[1]][, 1:(n.sigma + n.theta)]

print(matrix(colMeans(pm), nrow = 3))

pm2 = ans[[1]][, (n.sigma + n.theta + 1):(n.sigma + 2 * n.theta)]
print(matrix(colMeans(pm2), nrow = 3))