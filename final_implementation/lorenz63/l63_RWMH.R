set.seed(1)
library(mvtnorm)
#library(matrixcalc)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

# logarithm of the unnormalized posterior
#ludfun <- function(state) {
    ## State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    ## \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    #X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    #theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #Sigma = diag(state[(n.X + n.theta + 1):n.param], 3)

    #if (min(diag(Sigma)) < 0)
        #return(-Inf)

    ## Extracting observed data
    #X_t = X_n[, seq(2, N + 1, N / K)]

    ########################################################################
    #p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))) + dmvnorm(t(X_n[, 1] - tau_o), sigma = lam_o, log = TRUE)
    ##- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    #######################################################################
    #p2 = dmvnorm(theta, mean = prior_mean, sigma = prior_cov, log = TRUE)
    ##p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    #del_X = t(diff(t(X_n)))
    #p3 = sum(dmvnorm(t(del_X / del_t - f[, - (N + 1)]), sigma = Sigma, log = TRUE))
    ### add inverse gamma
    #p4 <- sum(dinvgamma(diag(Sigma), shape = a4, rate = b4, log = TRUE))
    #return(p1 + p2 + p3 + p4)

#}

ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #Sigma_vec = state[(3 * N + 7):(3 * N + 12)]
    Sigma = diag(state[(n.X + n.theta + 1):n.param], 3)

    # all the elements of theta should be positive
    #if (min(theta) <= 0)
    #return(-Inf)
    # \Sigma should be positive semi-definite
    if (min(diag(Sigma)) < 0)
        return(-Inf)

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]


    # pi is the log of likelihood
    # This doesn't need a loop
    #p1 = 0
    #for (k in 1:K) {
        #Y.t = t(t(Y[, k]))
        #X_t.t = t(t(X_t[, k]))
        #p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    #}
    #p1 = -0.5 * p1

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
    p3 = dmvnorm(theta, mean = prior_mean, sigma = prior_cov, log = TRUE)

    ## add inverse gamma
    p4 <- sum(dinvgamma(diag(Sigma), shape = a4, rate = b4, log = TRUE))
    return(p1 + p2 + p3 + p4)

}

MH <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale <- rep(.0001, n.param) #0.00005
    scale[(n.X + 1):(n.X + n.theta)] <- .01 #0.01
    scale[(n.X + n.theta + 1):(n.param)] <- .025 #0.025
    accept.prob = 0
    state = init

    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

        #if (i <= floor((4 * n) / 5)) {
        #scale[1:n.X] = .0022 - (0.0022 - 0.0005) * (5 * i) / (4 * n)
        #scale[(n.X + 1):(n.X + n.theta)] = .38 - (0.38 - .06) * (5 * i) / (4 * n)
        #}

        chain = metrop(ludfun, init, 1, scale = scale)
        state = chain$batch
        accept.prob = accept.prob + chain$accept
        X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
        theta = state[(n.X+1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        Sigma = state[(n.X + n.theta + 1):n.param]
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1:n.theta] = theta

        param_mat[i, (n.theta+1):(n.sigma + n.theta)] = Sigma
        init = state
    }
    Xfinal = state[1:n.X]
    print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg, Xfinal)
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
Nobs = 20 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
a4 = 10 #2
b4 = (a4 - 1) * 0.6 #.5
mu_truth = c(10, 28, 8 / 3)

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 #/ del_t
R = diag(.05, 3) # observational error
inv_R = diag(1 / (0.05), 3)
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
n.X = 3* (N + 1)
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma
n = 1e4

load('interp/l63_linch_T_20_1e5_cwise_1_spikes_interp_diffuse_6_by_10_scale_try')
pm = colMeans(to_save[[1]][[1]])
prior_mean = c(pm[4], pm[2], - pm[9])
prior_cov = diag(3, 3)


load('burninX6_by_10')
X = X[, 1:(N + 1)]
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.param)
init[1:n.X] = as.numeric(X)
init[(n.X + 1):(n.X + n.theta)] = mu_truth
init[(n.X+n.theta+1):n.param] = 0.6


#scale[c(6007, 6010, 6012)] <- 100
ans = MH(n, init)
pm = ans[[1]]
colMeans(pm)

#pdf(file = 'RWMH.pdf')
#plot.ts(pm[(n/4):n,])
#dev.off()
save(ans, file = "l96_1e6_MH")
