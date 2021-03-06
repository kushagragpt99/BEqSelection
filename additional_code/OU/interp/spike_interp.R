set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

make_tilde <- function(X, t) {
    X_vec = c(X, X ^ 2, t, t ^ 2)
    return(X_vec)
}
# drifet function for Lorenz-63
drift_fun <- function(X, t, B) {
    #print(make_tilde(X,t))
    tildeX = matrix(make_tilde(X, t), nrow = 4, ncol = 1)
    B_mat = matrix(B, nrow = 1)
    ans = B_mat %*% tildeX
    return(ans)
}

drift_fun_true <- function(X, theta) {
    ans = -theta * X
    return(ans)
}

ludfun <- function(state, gamma) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = state[1:n.X]
    B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 1)

    # Extracting observed data
    X_t = X_n[seq(2, N + 1, N / K)]

    p1 = sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE))
    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = diff(X_n)
    beta_tmp = sum((del_X / del_t - f[-(N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * log(b4 + beta_tmp)

    return(p1 + p2 + p3)

}

ludfun.X <- function(state, gamma, all) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
    all[1:n.X] = state
    X_n = all[1:n.X]
    B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 1)

    # Extracting observed data
    X_t = X_n[seq(2, N + 1, N / K)]

    p1 = sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE))

    B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
    p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
    #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

    #f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
    f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(B_vec))
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
    del_X = diff(X_n)
    beta_tmp = sum((del_X / del_t - f[-(N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * log(b4 + beta_tmp)

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
            print(accept.prob / i)
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

        X_n = state[1:n.X]
        theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1:n.theta] = theta

        f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(theta))
        del_X = diff(X_n)
        beta_tmp = sum((del_X / del_t - f[-(N + 1)]) ^ 2) * del_t / 2
        Sigma = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp)

        param_mat[i, (n.theta + 1):(n.theta + n.sigma)] = Sigma
        init = state
    }
    Xfinal = state[1:n.X]
    print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg, accept.prob / n, Xfinal)
    return(final_output)
}


# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = numeric(length = N + 1)
    X[1] = X0
    for (i in 2:(N + 1))
        X[i] = X[i - 1] + (drift_fun_true(X[i - 1], theta)) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


# hyper-parameters
to = 0 # initial time
tf = 2 # final time
Nobs = 20 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
a4 = 2
b4 = 1
tau1 = 10
tau0 = 0.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1) 
burn_in = 5000 #/ del_t
R = 0.04 #diag(1 / 2, 3) # observational error
inv_R = 1 / (0.04) #solve(R)
mu = 0
sigma2 = 10
mu_truth = c(-2, 0, 0, 0) #c(-10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
non_zero = 1
param_i = 1
n.X = (N + 1)
n.theta = 4
n.sigma = 1
n.param = n.X + n.theta + n.sigma
q = rep(0.1, n.theta) #runif(n.theta)
q[non_zero] = 0.9
n <- 1e5

X_total = euler_maruyama(-0.5, del_t, burn_in+N, 2, 1) # generating sample from Lorenz-63
X = X_total[(burn_in):(N + burn_in)]
#load('../../burninX')
Y = X[seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] = rmvnorm(1, mu_truth + rnorm(n.theta, sd = 1.5), sigma = diag(1 / 50, n.theta))

X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
y.index = 1
#X.interp[,1] = X[,1]
for (i in seq(2, N + 1, N / K)) {
    if (i == 2) {
        X.interp[2] = Y[1]
    } else {
        X.interp[(i - N / K + 1):i] = seq(Y[y.index], Y[y.index + 1], (Y[y.index + 1] - Y[y.index]) * K / N)[-1]
        y.index = y.index + 1
    }

}

init[(1:n.X)] <- as.numeric(X.interp)


sigma_Y = var(Y)
tau0 = sqrt(sigma_Y / (10 * K)) * 25
tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K))) * 5


load('../ou_linch_tf_1e4_1')
var1 = cov(ans2[[1]][, 1:n.theta])
scale_vec = 1 * sqrt(diag(var1))
scale = scale_vec
scale[c(1)] = 6 * scale_vec[c(1)] #5.5
scale[c(2)] = 2.4 * scale_vec[c(2)] #2
scale[c(3)] = 5 * scale_vec[c(3)] #4.7
scale[c(4)] = 10 * scale_vec[c(4)] #6.6

scale.X = 0.015 #0.0016


ans = linchpin(n, init)
#plot.ts(ans[[1]][, param_i])
#plot.ts(ans[[1]][, non_zero])
chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from interp ", "\n priors spike slab ", " time period ",
                            tf, " Sigma is .1"))

print(chain_info)
attr = list('to' = to, 'tf' = tf, 'Nobs' = Nobs, 'del_t' = del_t, 'a4' = a4, 'b4' = b4,
            'K' = K, 'N' = N, 'burn_in' = burn_in, 'R' = R, 'inv_R' = inv_R, 'mu_truth' = mu_truth, 'non_zero' = non_zero, 'param_i' = param_i,
            'n.X' = n.X, 'n.theta' = n.theta, 'n.sigma' = n.sigma, 'n.param' = n.param, 'q' = q, 'seq.Y' = seq.Y, 'n' = n, 'scale_vec' = scale_vec,
            'scale' = scale, 'scale.X' = scale.X)

to_save = list(ans, chain_info)
save(to_save,attr, file = "ou_linch_spike_interp_T2_1e5_try")
pm = ans[[1]][, 1:(n.sigma + n.theta)]

print(colMeans(pm))

pm2 = ans[[1]][, (n.sigma + n.theta + 1):(n.sigma + 2 * n.theta)]
print(colMeans(pm2))

#pdf(file = 'spike_interp.pdf')
#plot.ts(pm[(n / 4):n,])
#dev.off()