set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

# logarithm of the unnormalized posterior
ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #Sigma_vec = state[(3 * N + 7):(3 * N + 12)]
    Sigma = diag(state[(n.X + n.theta + 1):n.param], 3)

    # all the elements of theta should be positive
    if (min(theta) <= 0)
        return(-Inf)
    # \Sigma should be positive semi-definite
    if (min(diag(Sigma)) < 0)
        return(-Inf)

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

    #######################################################################
    p1 = sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))
    ######################################################################

    # p2 is the log of prior of X conditional on theta
    p2 = 0
    inv_Sig = solve(Sigma)
    for (k in 1:N) {
        del_X = matrix(X_n[, k + 1] - X_n[, k], nrow = 3, ncol = 1) # try using diff() function
        f_k = drift_fun(X_n[, k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + t(del_X / del_t - f_k) %*% inv_Sig %*% (del_X / del_t - f_k)
    }
    p2 = -0.5 * p2 * del_t

    ########################################################################
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    #del_X = t(diff(t(X_n)))
    #p2 = sum(dmvnorm(t(del_X - f[, - (N + 1)] * del_t), sigma = Sigma * del_t, log = TRUE))
    ########################################################################

    # store inv.lam_o globally
    p2 = p2 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o) - (N / 2) * determinant(Sigma * del_t, logarithm = TRUE)$modulus

    # p3 is the log of priors of theta
    p3 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    ## add inverse gamma
    p4 <- sum(dinvgamma(diag(Sigma), shape = 2, scale = 1 / 6, log = TRUE))
    return(p1 + p2 + p3 + p4)

}

MwG_update <- function(state, h) {
    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    Sigma = state[(n.X + n.theta + 1):n.param]
    Sigma_mat = diag(Sigma, 3)
    #X_t = X_n[, seq(2, N + 1, N / K)]
    #X_n = t(X_n)


    K_mwg = ((del_t / Sigma[1]) * (theta[1] - 1 / del_t) ^ 2 + (del_t / Sigma[2]) * (theta[2] - X_n[3, 1]) ^ 2
         + (del_t / Sigma[3]) * X_n[2, 1] ^ 2 + 1)
    C_mwg = ((theta[1] * del_t - 1) / Sigma[1] * (theta[1] * X_n[2, 1] - X_n[1, 2] / del_t) + (theta[2] - X_n[3, 1]) / Sigma[2] *
         (X_n[2, 2] - X_n[2, 1] * (1 - del_t)) + (X_n[2, 1] / Sigma[3]) * (X_n[3, 2] - X_n[3, 1] * (1 - theta[3] * del_t)))
    X_n[1, 1] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))

    K_mwg = ((theta[1] ^ 2 * del_t) / Sigma[1] + (del_t / Sigma[2]) * (1 - 1 / del_t) ^ 2 + del_t * X_n[1, 1] ^ 2 / Sigma[3] + 1)
    C_mwg = (theta[1] / Sigma[1] * (X_n[1, 2] - X_n[1, 1] * (1 - theta[1] * del_t)) + (1 - del_t) / Sigma[2] * (X_n[2, 2] / del_t
        + X_n[1, 1] * (X_n[3, 1] - theta[2])) + X_n[1, 1] / Sigma[3] * (X_n[3, 1] * (theta[3] * del_t - 1) + X_n[3, 2]))
    X_n[2, 1] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))

    K_mwg = del_t * X_n[1, 1] ^ 2 / Sigma[2] + del_t / Sigma[3] * (theta[3] - 1 / del_t) ^ 2 + 1
    C_mwg = ((del_t * X_n[1, 1]) / Sigma[2] * (theta[2] * X_n[1, 1] - X_n[2, 1] - (X_n[2, 2] - X_n[2, 1]) / del_t)
        + (1 - theta[3] * del_t) / Sigma[3] * (X_n[3, 2] / del_t - X_n[1, 1] * X_n[2, 1]))
    X_n[3, 1] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))

    t_i = 2
    for (i in 2:N) {

        K_mwg = (1 / (Sigma[1] * del_t ^ 2) + (theta[1] - 1 / del_t) ^ 2 / Sigma[1] + (theta[2] - X_n[3, i]) ^ 2 / Sigma[2]
            + X_n[2, i] ^ 2 / Sigma[3])
        C_mwg = (1 / (del_t * Sigma[1]) * (X_n[1, i - 1] / del_t + theta[1] * (X_n[2, i - 1] - X_n[1, i - 1])) + (theta[1] - 1 / del_t) / Sigma[1]
             * (theta[1] * X_n[2, i] - X_n[1, i + 1] / del_t) + (theta[2] - X_n[3, i]) / Sigma[2] * (X_n[2, i] + (X_n[2, i + 1] - X_n[2, i]) / del_t)
             + X_n[2, i] / Sigma[3] * (theta[3] * X_n[3, i] + (X_n[3, i + 1] - X_n[3, i]) / del_t))

        if (i == t_i) {
            #t_i = t_i + N / K
            y.index = (i - 2) * (K / N) + 1
            C_mwg = C_mwg * del_t + Y[1, y.index] / R[1, 1]
            K_mwg = K_mwg * del_t + 1 / R[1, 1]
            X_n[1, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))
        } else {
            X_n[1, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / (K_mwg * del_t)))
        }

        K_mwg = 1 / (Sigma[2] * del_t ^ 2) + theta[1] ^ 2 / Sigma[1] + (1 - 1 / del_t) ^ 2 / Sigma[2] + X_n[1, i] ^ 2 / Sigma[3]
        C_mwg = (1 / (del_t * Sigma[2]) * (X_n[2, i - 1] / del_t + X_n[1, i - 1] * (theta[2] - X_n[3, i - 1]) - X_n[2, i - 1])
            + theta[1] / Sigma[1] * (theta[1] * X_n[1, i] + (X_n[1, i + 1] - X_n[1, i]) / del_t) + (1 - 1 / del_t) / Sigma[2]
            * (X_n[1, i] * (theta[2] - X_n[3, i]) - X_n[2, i + 1] / del_t) + X_n[1, i] / Sigma[3] * (theta[3] * X_n[3, i]
            + (X_n[3, i + 1] - X_n[3, i]) / del_t))

        if (i == t_i) {
            #t_i = t_i + N / K
            y.index = (i - 2) * (K / N) + 1
            C_mwg = C_mwg * del_t + Y[2, y.index] / R[2, 2]
            K_mwg = K_mwg * del_t + 1 / R[2, 2]
            X_n[2, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))
        } else {
            X_n[2, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / (K_mwg * del_t)))
        }

        K_mwg = 1 / (Sigma[3] * del_t ^ 2) + X_n[1, i] ^ 2 / Sigma[2] + (theta[3] - 1 / del_t) ^ 2 / Sigma[3]
        C_mwg = (1 / (Sigma[3] * del_t) * (X_n[3, i - 1] / del_t + X_n[1, i - 1] * X_n[2, i - 1] - theta[3] * X_n[3, i - 1])
            + X_n[1, i] / Sigma[2] * (theta[2] * X_n[1, i] - X_n[2, i] - (X_n[2, i + 1] - X_n[2, i]) / del_t)
            + (theta[3] - 1 / del_t) / Sigma[3] * (X_n[1, i] * X_n[2, i] - X_n[3, i + 1] / del_t))

        if (i == t_i) {

            #t_i = floor(t_i + N / K)
            y.index = (i - 2) * (K / N) + 1
            t_i = 2 + floor(y.index * N / K)
            C_mwg = C_mwg * del_t + Y[3, y.index] / R[3, 3]
            K_mwg = K_mwg * del_t + 1 / R[3, 3]
            X_n[3, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / K_mwg))
        } else {
            X_n[3, i] = rnorm(1, C_mwg / K_mwg, sqrt(1 / (K_mwg * del_t)))
        }

    }

    X_n[, N + 1] = rmvnorm(1, X_n[, N] - t(drift_fun(X_n[, N], theta)), del_t * Sigma_mat)

    #X_n = t(X_n)
    # update for theta
    prop = rmvnorm(1, mean = theta, sigma = h * diag(1, 3))
    prop_vec = c(as.vector(X_n), prop, Sigma)
    old_vec = c(as.vector(X_n), theta, Sigma)
    #print(ludfun(old_vec))
    u = runif(1)
    if (log(u) < (ludfun(prop_vec) - ludfun(old_vec))) {
        theta = prop
        accept_prob_theta = accept_prob_theta + 1
    }

    f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    Sigma[1] = rgamma(1, shape = N / 2 + 2, rate = 6 + beta_tmp[1])
    Sigma[2] = rgamma(1, shape = N / 2 + 2, rate = 6 + beta_tmp[2])
    Sigma[3] = rgamma(1, shape = N / 2 + 2, rate = 6 + beta_tmp[3])

    update = c(as.vector(X_n), theta, Sigma)
    return(update)

}


MwG <- function(init, n, h) {
    X_avg = numeric(length = n.X)
    old_update = init
    param_mat = matrix(, nrow = n, ncol = 6)

    for (t in 1:n) {
        #print(t)
        new_update = MwG_update(old_update, h)
        param_mat[t,] = new_update[(n.X + 1):n.param]
        X_avg = X_avg * (t - 1) / t + new_update[1:n.X] / t
        old_update = new_update
    }
    final_output = list(param_mat, X_sum)
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

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
R = diag(2, 3) # observational error
inv_R = solve(R)
n.X = 3 * (N + 1)
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma

accept_prob_theta = 0


X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.param)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- c(10, 28, 8 / 3) # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = 6 # inital \Sigma should also be positive semi definite

ans = MwG(init, 1e2, 0.2)

