set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

make_tilde <- function(X, t) {
    X_vec = c(1, X[1], X[2], X[3], X[4], X[1] ^ 2, X[2] ^ 2, X[3] ^ 2, X[4] ^ 2, X[1] * X[2], X[1] * X[3], X[1] * X[4],
              X[2] * X[3], X[2] * X[4], X[3] * X[4], t, t ^ 2)
    return(X_vec)
}
# drifet function for Lorenz-63
drift_fun <- function(X, t, B) {
    #print(make_tilde(X,t))
    tildeX = matrix(make_tilde(X, t), nrow = 17, ncol = 1)
    B_mat = matrix(B, nrow = N.l96)
    ans = B_mat %*% tildeX
    return(ans)
}

drift_fun_true <- function(X, theta) {
    ans = matrix(, nrow = N.l96, ncol = 1)
    for (i in 0:(N.l96 - 1)) {
        ans[i + 1, 1] = (X[(i + 1) %% N.l96 + 1] - X[(i - 2) %% N.l96 + 1]) * X[(i - 1) %% N.l96 + 1] - X[i + 1] + theta
    }
    return(ans)
}

ludfun <- function(state, gamma) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
    B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = N.l96)

    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE)))
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################
    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3


    f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

ludfun.X <- function(state, gamma, all) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
    all[1:n.X] = state
    X_n = matrix(all[1:n.X], nrow = N.l96, ncol = N + 1)
    B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = N.l96)

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE)))
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################
    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    #f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
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

        if (i %% (n / 5) == 0) {
            print(i)
            print(accept.prob[1]/i)
            print(matrix(accept.prob[2:(n.theta + 1)] / i, nrow = N.l96))
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

        X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
        theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1:n.theta] = theta

        Sigma = numeric(length = n.sigma)
        f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(theta))
        del_X = t(diff(t(X_n)))
        beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
        for (j in 1:n.sigma) {
            Sigma[j] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[j])
        }

        param_mat[i, (n.theta + 1):(n.theta + n.sigma)] = Sigma
        init = state
    }
    #print(accept.prob / n)
    save.X = state[1:n.X]
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg, save.X,accept.prob / n)
    return(final_output)
}


# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = N.l96, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun_true(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


# hyper-parameters
to = 0 # initial time
tf = 10 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt
a4 = 2
b4 = .5
tau1 = 10
tau0 = 0.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 #/ del_t
R = diag(.05, N.l96) # observational error
inv_R = diag(1 / (0.05), N.l96)
mu = 0
sigma2 = 10
mu_truth = c(rep(8, 4), as.numeric(diag(rep(-1, 4))), rep(0, 16), 0, 0, -1, 0, 0, 1, 0, 1, 0, -1, rep(0, 5), -1, 1, 0, 1, 0, -1, rep(0, 11)) #c(-10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
non_zero = c(1, 2, 3, 4, 5, 10, 15, 20, 39, 42, 44, 46, 52, 53, 55, 57)
param_i = 1:4
n.theta = 68
n.sigma = N.l96
q = rep(0.1, n.theta) #runif(n.theta)
q[non_zero] = 0.9
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
n.X = N.l96 * (N + 1)
n.param = n.X + n.theta + n.sigma

n <- 1e5



X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.5, N.l96)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
X = X[, 1:(N + 1)]
Y = X[, seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)

init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1, mu_truth + rnorm(n.theta, sd = 1.5), sigma = diag(1 / 50, n.theta)) #rmvnorm(1, c(10, 28, 8 / 3), sigma = diag(0.5, 3)) # random initial values for MCMC


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


sigma_Y = mean(diag(var(t(Y))))
tau0 = sqrt(sigma_Y / (10 * K)) * 2
tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K))) / 4

tau0 = sqrt(sigma_Y / (10 * K)) * 1.5
tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K))) / 2

load('l96_linch_spike_5e2')
var1 = cov(to_save[[1]][[1]][, 1:n.theta])
scale_vec = 1.2 * sqrt(diag(var1))

scale = rep(n.theta)

scale = scale_vec # tf = 5 Sigma = 0.1  change of above
scale[c(1)] = .4 * scale_vec[c(1)] ## Sigma = .5
scale[c(2)] = 3.2 * scale_vec[c(2)]
scale[c(3)] = 12 * scale_vec[c(3)]
scale[c(4)] = 5.7 * scale_vec[c(4)]
scale[c(5)] = 0.7 * scale_vec[c(5)]
scale[c(6)] = 3.5 * scale_vec[c(6)]
scale[c(14, 16)] = 3 * scale_vec[c(14, 16)]
scale[c(9, 13, 18, 62)] = 2.5 * scale_vec[c(9, 13, 18, 62)]
scale[c(10)] = 0.4 * scale_vec[c(10)]
scale[c(8, 12, 29, 32, 55)] = 0.5 * scale_vec[c(8, 12, 29, 32, 55)]
scale[c(67)] = 0.3 * scale_vec[c(67)]
scale[c(21, 24, 25, 26, 27, 28, 30, 31, 35, 65)] = 0.6 * scale_vec[c(21, 24, 25, 26, 27, 28, 30, 31, 35, 65)]
scale[c(17, 19)] = 2.5 * scale_vec[c(17, 19)]
scale[c(7, 20, 61, 64)] = 2 * scale_vec[c(7, 20, 61, 64)]
scale[c(23)] = 1.8 * scale_vec[c(23)]
scale[c(49, 52)] = 0.7 * scale_vec[c(49, 52)]
scale[c(40, 42, 43)] = 1.2 * scale_vec[c(40, 42, 43)]
scale[c(22, 37, 38, 41, 45, 46, 47, 48, 49, 50, 54, 56, 58, 59, 60, 61)] = 1.5 * scale_vec[c(22, 37, 38, 41, 45, 46, 47, 48, 49, 50, 54, 56, 58, 59, 60, 61)]
scale[c(11, 15, 39, 53, 66, 68)] = 0.8 * scale_vec[c(11, 15, 39, 53, 66, 68)]

scale.X = 0.0018

# These values are generate ffrom the same code with difference scale and scale.X, which can be found in attr loaded with this file
load('l96_5e5_cwise_spikes_interp_diffuse_init_theta_try')
init[(n.X + 1):(n.X + n.theta)] = colMeans(to_save[[1]][[1]][3e4:5e4, 1:(n.theta)])
#init[1:n.X] = to_save[[1]][[2]]

ans = linchpin(n, init)

chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from previous run ", "\n priors spike slab ", " time period ",
                                tf, " Sigma is .5"))

print(chain_info)

attr = list('to' = to, 'tf' = tf, 'Nobs' = Nobs, 'N.l96' = N.l96, 'del_t' = del_t, 'a4' = a4, 'b4' = b4, 'tau0' = tau0, 'tau1' = tau1,
            'K' = K, 'N' = N, 'burn_in' = burn_in, 'R' = R, 'inv_R' = inv_R, 'mu_truth' = mu_truth, 'non_zero' = non_zero, 'param_i' = param_i,
            'n.X' = n.X, 'n.theta' = n.theta, 'n.sigma' = n.sigma, 'n.param' = n.param, 'q' = q, 'seq.Y' = seq.Y, 'n' = n, 'scale_vec' = scale_vec,
            'scale' = scale, 'scale.X' = scale.X)

to_save = list(ans, chain_info)
save(to_save,attr, file = "l96_1e5_cwise_spikes_interp_diffuse_init_theta_try")
pm = ans[[1]][, 1:(n.sigma + n.theta)]

print(matrix(colMeans(pm), nrow = N.l96))

pm2 = ans[[1]][, (n.sigma + n.theta + 1):(n.sigma + 2 * n.theta)]
print(matrix(colMeans(pm2), nrow = N.l96))