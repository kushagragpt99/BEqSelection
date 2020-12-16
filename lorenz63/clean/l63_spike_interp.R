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
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
    #if (index == 0) {
    ##print('0')
    #all[1:n.X] = state
    #} else {
    ##print(index)
    #all[n.X+index] = state
    #}

    #X_n = matrix(all[1:n.X], nrow = 3, ncol = N + 1)
    #B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #B_mat = matrix(B_vec, nrow = 3)

    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 3)

    # all the elements of theta should be positive
    #if (min(theta) <= 0)
    #return(-Inf)

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
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
    #if (index == 0) {
    ##print('0')
    #all[1:n.X] = state
    #} else {
    ##print(index)
    #all[n.X+index] = state
    #}
    all[1:n.X] = state
    X_n = matrix(all[1:n.X], nrow = 3, ncol = N + 1)
    B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    B_mat = matrix(B_vec, nrow = 3)

    #X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    #B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #B_mat = matrix(B_vec, nrow = 3)

    # all the elements of theta should be positive
    #if (min(theta) <= 0)
    #return(-Inf)

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

linchpin <- function(n, init, scale_vec) {
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
            print(matrix(accept.prob[2:(n.theta + 1)] / i, nrow = 3))
            #print(c(i, accept.prob / i))
        }

        all = init
        chain = metrop(ludfun.X, initial = init[1:n.X], nbatch = 1, scale = scale.X, gamma = gamma, all = all)
        accept.prob[1] = accept.prob[1] + chain$accept
        state[1:n.X] = chain$batch

        #ans = MH.X(init[1:n.X], 1, scale.X, gamma, init[(n.X + 1):(n.X + n.theta)])
        #accept.prob[1] = accept.prob[1] + ans[[2]]
        #init[1:n.X] = ans[[1]]

        for (j in 1:n.theta) {
            ####### are you considering updated values of B for next steps??????
            #all = init
            #chain = metrop(ludfun, initial = init[n.X + j], nbatch = 1, scale = scale.B[j], gamma = gamma, all = all, index = j)
            #accept.prob[j + 1] = accept.prob[j + 1] + chain$accept
            #init[n.X + j] = chain$batch

            ans = MH.B(j, init[n.X + j], 1, scale.B[j], gamma, state)
            accept.prob[j + 1] = accept.prob[j + 1] + ans[[2]]
            state[n.X + j] = ans[[1]]
        }

        #chain = metrop(ludfun, init, 1, scale = scale, gamma = gamma)
        #state = chain$batch
        #accept.prob = accept.prob + chain$accept

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
alpha1 = 20 # Prior for \sigma is Gamma (alpha1, beta1)
alpha2 = 56 # Prior for \rho is Gamma (alpha2, beta2)
alpha3 = 6 # Prior for \beta is Gamma (alpha3, beta3)
beta1 = 0.5
beta2 = 0.5
beta3 = 0.5
a4 = 2
b4 = 0.6
#tau1 = 10/2
#tau0 = 0.5*2

tau1 = 10/2  # just trying
tau0 = 0.5 # 0.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
burn_in = 5000 / del_t
R = diag(0.05, 3) # observational error
inv_R = diag(1/(0.05), 3)
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
n <- 5e5

X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(.6, 3)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
#load('burninX')
X = X[, 1:(N + 1)]
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)

#init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)

init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1, mu_truth + rnorm(n.theta), sigma = diag(1 / 50, n.theta))
load('l63_linch_T_20_1e4_cwise_1_spikes_interp_diffuse_6_by_10_scale_try')
ans = to_save[[1]]
pm = ans[[1]][, 1:(n.sigma + n.theta)]
init[(n.X + 1):(n.X + n.theta)] = colMeans(pm[8e3:1e4,1:n.theta])


X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
y.index = 1
#X.interp[,1] = X[,1]
for (i in seq(2, N + 1, N / K)) {
    if (i == 2) {
        X.interp[, 2] = Y[, 1]
    } else {
        X.interp[1, (i - N / K + 1):i] = seq(Y[1, y.index], Y[1, y.index + 1], (Y[1, y.index + 1] - Y[1, y.index]) * K / N)[-1]
        X.interp[2, (i - N / K + 1):i] = seq(Y[2, y.index], Y[2, y.index + 1], (Y[2, y.index + 1] - Y[2, y.index]) * K / N)[-1]
        X.interp[3, (i - N / K + 1):i] = seq(Y[3, y.index], Y[3, y.index + 1], (Y[3, y.index + 1] - Y[3, y.index]) * K / N)[-1]
        y.index = y.index + 1
    }

}

init[(1:n.X)] <- as.numeric(X.interp)
init[1:n.X] = ans[[2]]

#sigma_Y = mean(diag(var(t(Y))))
#tau0 = sqrt(sigma_Y / (10 * K)) 
#tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K))) /2 # 2

load('l63_linch_T_20_5e5_1')
var1 = cov(to_save[[1]][[1]][, 1:33])
scale_vec = 1.8 * sqrt(diag(var1))
scale = scale_vec
scale[ c(1)] = 0.7 * scale_vec[c(1)]
scale[ c(2)] = 1.2 * scale_vec[c(2)]
scale[ c(3)] = 2.0 * scale_vec[c(3)]
scale[ c(4)] = 1.2 * scale_vec[c(4)]
scale[ c(5)] = 1.3 * scale_vec[c(5)]
scale[ c(6)] = 0.7 * scale_vec[c(6)]
scale[ c(7)] = 0.4 * scale_vec[c(7)]
scale[ c(8)] = 0.55 * scale_vec[c(8)]
scale[ c(9)] = 0.22 * scale_vec[c(9)]
scale[ c(10)] = 0.25 * scale_vec[c(10)]
scale[ c(11)] = 0.3 * scale_vec[c(11)]
scale[ c(12)] = 0.25 * scale_vec[c(12)]
scale[ c(13)] = 0.35 * scale_vec[c(13)]
scale[ c(14)] = 0.5 * scale_vec[c(14)]
scale[ c(15)] = 0.4 * scale_vec[c(15)]
scale[ c(16)] = 0.22 * scale_vec[c(16)]
scale[ c(17)] = 0.4 * scale_vec[c(17)]
scale[ c(18)] = 0.2 * scale_vec[c(18)]
scale[ c(19)] = 0.2 * scale_vec[c(19)]
scale[ c(20)] = 0.25 * scale_vec[c(20)]
scale[ c(21)] = 0.2 * scale_vec[c(21)]
scale[ c(22)] = 0.6 * scale_vec[c(22)]
scale[ c(23)] = 0.7 * scale_vec[c(23)]
scale[ c(24)] = 0.6 * scale_vec[c(24)]
scale[ c(25)] = 0.7 * scale_vec[c(25)]
scale[ c(26)] = 1.2 * scale_vec[c(26)]
scale[ c(27)] = 1.1 * scale_vec[c(27)]
scale[ c(28)] = 0.5 * scale_vec[c(28)]
scale[ c(29)] = 0.7 * scale_vec[c(29)]
scale[ c(30)] = 0.5 * scale_vec[c(30)]
scale[ c(31)] = 0.45 * scale_vec[c(31)]
scale[ c(32)] = 0.7 * scale_vec[c(32)]
scale[ c(33)] = 0.7 * scale_vec[c(33)]


scale.X = 0.0032 #0.002 ##################################### tried increasing this?


ans = linchpin(n, init)
#plot.ts(ans[[1]][, param_i])
#plot.ts(ans[[1]][, non_zero])
chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from previous run ", "\n priors spike slab ", " time period ",
                            tf, " lam_0 is 10"))

print(chain_info)

attr = list('to' = to, 'tf' = tf, 'Nobs' = Nobs, 'del_t' = del_t, 'a4' = a4, 'b4' = b4, 'tau0' = tau0, 'tau1' = tau1,
            'K' = K, 'N' = N, 'burn_in' = burn_in, 'R' = R, 'inv_R' = inv_R, 'mu_truth' = mu_truth, 'non_zero' = non_zero, 'param_i' = param_i,
            'n.X' = n.X, 'n.theta' = n.theta, 'n.sigma' = n.sigma, 'n.param' = n.param, 'q' = q, 'seq.Y' = seq.Y, 'n' = n, 'scale_vec' = scale_vec,
            'scale' = scale, 'scale.X' = scale.X)

to_save = list(ans, chain_info)
save(to_save,attr, file = "l63_linch_T_20_5e5_cwise_1_spikes_interp_diffuse_6_by_10_scale_try")
pm = ans[[1]][, 1:(n.sigma + n.theta)]

print(matrix(colMeans(pm), nrow = 3))

pm2 = ans[[1]][, (n.sigma + n.theta + 1):(n.sigma + 2 * n.theta)]
print(matrix(colMeans(pm2), nrow = 3))

#100 minutes