set.seed(1)
library(mvtnorm)
library(matrixcalc)
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
    Sigma = diag(state[(n.X + n.theta + 1):n.param])

    # all the elements of theta should be positive
    if (min(theta) <= 0)
        return(-Inf)
    # \Sigma should be positive semi-definite
    if (is.positive.semi.definite(Sigma) == FALSE)
        return(-Inf)

    # Euler - Muryami approximation expansions
    X_t = X_n[, seq(2, N + 1, N / K)]


    # pi is the log of likelihood
    # This doesn't need a loop
    p1 = 0
    for (k in 1:K) {
        Y.t = t(t(Y[, k]))
        X_t.t = t(t(X_t[, k]))
        p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    }
    p1 = -0.5 * p1

    #######################################################################
    # p1 = sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))
    ######################################################################

    # p2 is the log of prior of X conditional on theta
    p2 = 0
    inv_Sig = solve(Sigma*del_t)
    for (k in 1:N) {
        del_X = matrix(X_n[, k + 1] - X_n[, k], nrow = 3, ncol = 1) # try using diff() function
        f_k = drift_fun(X[, k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + t(del_X / del_t - f_k) %*% inv_Sig %*% (del_X / del_t - f_k)
    }
    p2 = -0.5 * p2

    ########################################################################
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    #del_X = t(diff(t(X_n)))
    #p2 = sum(dmvnorm(t(del_X / del_t - f[, - (N + 1)]), sigma = Sigma * del_t, log = TRUE))
    ########################################################################

    # store inv.lam_o globally
    p2 = p2 - 0.5 * t( t(t(X_n[, 1])) - tau_o ) %*% inv.lam_o %*% ( t(t(X_n[, 1])) - tau_o ) - (N / 2) * determinant(Sigma * del_t, logarithm = TRUE)$modulus

    # p3 is the log of priors of theta
    p3 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    ## add inverse gamma
    p4 <- sum(dinvgamma(diag(Sigma), shape = 3, scale = 1/2, log = TRUE))
    return(p1 + p2 + p3 + p4)

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


X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- c(10, 28, 8 / 3) # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = c(6,6,6)                                      # inital \Sigma should also be positive semi definite

scale <- rep(.2, n.param)
#scale[c(6007, 6010, 6012)] <- 100
chain = metrop(ludfun, init, nbatch = 1e4, scale = scale) # running MH

chain$accept
out <- chain$batch[, (n.X + 1):n.param]
plot.ts(out)
print(colMeans(out))

