set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = -theta * X
    return(ans)
}

ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = state[1:n.X]
    theta = state[(n.X + 1)] # vector of \sigma, \rho and \beta    


    # all the elements of theta should be positive
    if (min(theta) <= 0)
        return(-Inf)

    # Extracting observed data
    X_t = X_n[seq(2, N + 1, N / K)]


    # pi is the log of likelihood
    # This doesn't need a loop
    #p1 = 0
    ##print(dim(Y))
    #for (k in 1:K) {
        #Y.t = Y[k]
        #X_t.t = X_t[k]
        #p1 = p1 + (Y.t - X_t.t) * inv_R * (Y.t - X_t.t)
    #}
    #p1 = -0.5 * p1
    #p1 = p1 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o)

    #######################################################################
    p1 = (sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE)))
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################

    #p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3
    p2 = (alpha - 1) * log(theta) - theta/beta

    f = sapply(X_n, drift_fun, theta)  #sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = diff(X_n)
    beta_tmp = sum((del_X / del_t - f[- (N + 1)]) ^ 2) * del_t / 2
    p3 =  - (a4 + N / 2) * log(b4 + beta_tmp)

    return(p1 + p2 + p3)

}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta+n.sigma)
    scale = rep(0.014, n.X + n.theta) # 0.0047 #0.005
    scale[(n.X + 1):(n.X + n.theta)] = 0.8
    accept.prob = 0

    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

        if (i == 2e4) {
            scale[1:n.X] = rep(0.008, n.X)
        }

        chain = metrop(ludfun, init, 1, scale = scale)
        state = chain$batch
        accept.prob = accept.prob + chain$accept
        X_n = state[1:n.X]
        theta = state[(n.X + 1)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1:n.theta] = theta


        f = sapply(X_n, drift_fun, theta) #sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
        del_X = diff(X_n)
        beta_tmp = sum((del_X / del_t - f[-(N + 1)]) ^ 2) * del_t / 2
        Sigma = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp)

        param_mat[i, n.theta+n.sigma] = Sigma
        init = state
    }

    print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg)
    return(final_output)
}


# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = numeric(length = N+1)
    X[1] = X0
    for (i in 2:(N + 1))
        X[i] = X[i - 1] + drift_fun(X[i - 1], theta) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


# hyper-parameters
to = 0 # initial time
tf = 2 # final time
Nobs = 50 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha = 4
beta = 0.5
a4 = 2
b4 = 1

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 0 #5000 / del_t
R = 0.04 #diag(2, 3) # observational error
inv_R = 1/(0.04) #solve(R)
n.X = (N + 1)
n.theta = 1
n.sigma = 1
n.param = n.X + n.theta + n.sigma
n = 1e5

X = euler_maruyama(-0.5, del_t, N , 2, 1) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
#load('../burninX')
Y = X[seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) #+ rnorm(n.X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- rnorm(1,2,sd = 0.2) #rmvnorm(1, c(10, 28, 8 / 3), sigma = diag(0.5, 3)) # random initial values for MCMC

X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
y.index = 1
#X.interp[,1] = X[,1]
for (i in seq(2, N + 1, N / K)) {
    if (i == 2) {
        X.interp[2] = Y[1]
    } else {
        X.interp[(i - N / K + 1):i] = seq(Y[y.index], Y[ y.index + 1], (Y[y.index + 1] - Y[ y.index]) * K / N)[-1]
        y.index = y.index + 1
    }

}

init[(1:n.X)] <- as.numeric(X.interp)


ans = linchpin(n, init)
pm = ans[[1]]
colMeans(pm)
#plot.ts(pm)
save(ans, file = "OU_linch_1e5_Nobs_init_tf")
# 11:30