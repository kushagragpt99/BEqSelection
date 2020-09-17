set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = matrix(, nrow = 40, ncol = 1)
    for (i in 0:39) {
        ans[i+1,1] = (X[(i+1)%%40 + 1] - X[(i-2)%%40 + 1]) * X[(i-1)%%40 + 1] - X[i + 1] + theta
    }
    return(ans)
}

make_Sigma <- function(sigma_ar, n) {
    Sigma = matrix(, nrow = n, ncol = n)
    k = 1
    start = 1
    for (i in 1:n) {
        len = n - i + 1
        if (i != 1)
            start = start + n - i + 2
        for (j in 1:len) {
            Sigma[i, (i - 1 + j)] = sigma_ar[start + j - 1]
            Sigma[(i - 1 + j), i] = sigma_ar[start + j - 1]
        }
    }
    return(Sigma)
}
#Sig = make_Sigma(sigma_init, 4)
#make_Sigma(1:10, 4)

# logarithm of the unnormalized posterior
ludfun <- function(state) {
    # State is the vector storing the vectors of length 40*N + 47. The first 40*(N+1) terms are Xs. The next term is the drift partameter \theta. 
    # The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = 40, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)]
    Sigma = diag(state[(n.X + n.theta + 1):n.param])

    # \Sigma should be positive semi-definite
    if (min(state[(n.X + n.theta + 1):n.param]) <= 0) {
        #print("mehh")
        return(-Inf)
    }

    X_t = X_n[, seq(2, N + 1, N / K)]

    # pi is the log of likelihood
    p1 = 0
    for (k in 1:K) {
        Y.t = t(t(Y[, k]))
        X_t.t = t(t(X_t[, k]))
        p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    }
    p1 = -0.5 * p1

    # p2 is the log of prior of X conditional on theta
    p2 = 0

    inv_Sig = solve(Sigma)
    for (k in 1:N) {
        del_X = matrix(X_n[, k + 1] - X_n[, k], nrow = 40, ncol = 1)
        f_k = drift_fun(X_n[, k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + t(del_X / del_t - f_k) %*% inv_Sig %*% (del_X / del_t - f_k)
    }
    p2 = -0.5 * p2
    p2 = p2 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% solve(lam_o) %*% (t(t(X_n[, 1])) - tau_o) - (N / 2) * log(det(Sigma * del_t))

    # p3 is the log priors of theta and \Sigma
    p3 = (theta - mu)^2 / (2 * sigma^2) + sum(dinvgamma(diag(Sigma), shape = 3, scale = 2, log = TRUE))

    return(p1 + p2 + p3)

}

# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = 40, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + det(del_t * Sigma) * rmvnorm(1, sigma = diag(1, 40))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

component_wise_MH <- function(init, iter, hx, htheta, hsigma) {

    X_n = matrix(state[1:n.X], nrow = 40, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)]
    Sigma = diag(state[(n.X + n.theta + 1):n.param])

    MH_chain = matrix(, nrow = iter, ncol = nparam)
    MH_chain[1,] = init
    accept.prob = rep(0, n.sigma + 2)

    for (i in 2:iter) {

        # all Xs treated as one component
        prop.X = c(rmvnorm(1, mean = MH_chain[i - 1, 1:n.X], sigma = diag(hx, n.X)), MH_chain[i - 1, (n.X + 1):n.param])
        ux <- runif(1)
        if (log(ux) < (ludfun(prop.X) - ludfun(MH_chain[i - 1,]) )) {
            MH_chain[i,] = prop.X
            accept.prob[1] = accept.prob[1] + 1
        }
        else
            MH_chain[i,] = MH_chain[i - 1,]

        # updating theta
        prop.theta = c(MH_chain[i, (n.X)], rmvnorm(1, mean = MH_chain[i, (n.X + 1):(n.X + n.theta)], sigma = diag(htheta, n.theta)),
                       MH_chain[i, (n.X + n.theta + 1):n.param])
        uh <- runif(1)
        if (log(ux) < (ludfun(prop.theta) - ludfun(MH_chain[i,]))) {
            MH_chain[i, ] = prop.theta
            accept.prob[2] = accept.prob[2] + 1
        }

        #updating sigma
        prop.sigma_ar = MH_chain[i, (n.X + n.theta + 1):n.param]
        for (k in 1:n.sigma) {
            prop.sigma_comp = rnorm(1, mean = MH_chain[i, (n.X + n.theta + k)], sd = sqrt(hsigma))
            temp = prop.sigma_ar[k]
            prop.sigma_ar[k] = prop.sigma_comp
            if (is.positive.semi.definite(make_Sigma(prop.sigma_ar)) == FALSE) {
                prop.sigma_ar[k] = temp
            }
            else {
                usigma = runif(1)
                prop.sigma = c(MH_chain[i, (n.X + n.theta + 1):n.param], prop.sigma_ar)
                if (log(usigma) < (ludfun(prop.sigma) - ludfun(MH_chain[i, ]))) {
                    MH_chain[i,] = prop.sigma
                    accept.prob[2+k] = accept.prob[2+k] + 1
                }
                else {
                    prop.sigma_ar[k] = temp
                }
            }
        }
    }
    ans = list(MH_chain, accept.prob)
    return(ans)
}

# hyper-parameters
to = 0 # initial time
tf = 10 # final time
Nobs = 8 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 40), nrow = 40, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(1, 40) # prior covariance matrix of X[0]
mu = 8
sigma = 2

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
R = diag(1, 40) # observational error
inv_R = solve(R)

n.X = 40 * N + 40
n.theta = 1
n.sigma = 40
n.param = n.X + n.theta + n.sigma

X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, 8, diag(2, 40)) # generating sample from Lorenz-63
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 40), sigma = R)) # observations from Lorenz-63
init = runif(n.param, 0, 5) # random initial values for MCMC

chain = metrop(ludfun, init, nbatch = 1e4, scale = 0.1) # running MH
out = chain$batch
print(chain$accept)
print(mean(out[, n.X + 1]))
