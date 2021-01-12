set.seed(1)
library(mvtnorm)
#library(matrixcalc)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = matrix(, nrow = N.l96, ncol = 1)
    for (i in 0:(N.l96 - 1)) {
        ans[i + 1, 1] = (X[(i + 1) %% N.l96 + 1] - X[(i - 2) %% N.l96 + 1]) * X[(i - 1) %% N.l96 + 1] - X[i + 1] + theta
    }
    return(ans)
}

# logarithm of the unnormalized posterior
ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    Sigma = diag(state[(n.X + n.theta + 1):n.param], N.l96)

    if (min(diag(Sigma)) < 0)
        return(-Inf)

    # Extracting observed data
    X_t = X_n[, seq(2, N + 1, N / K)]

    #######################################################################
    p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))) + dmvnorm(t(X_n[, 1] - tau_o), sigma = lam_o, log = TRUE)
    #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################
    p2 = dnorm(theta, mean = alpha, sd = beta, log = TRUE)
    #p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = t(diff(t(X_n)))
    p3 = sum(dmvnorm(t(del_X / del_t - f[, - (N + 1)]), sigma = Sigma, log = TRUE))
    ## add inverse gamma
    p4 <- sum(dinvgamma(diag(Sigma), shape = a4, rate = b4, log = TRUE))
    return(p1 + p2 + p3 + p4)

}

MH <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale <- rep(.0015, n.param)
    scale[(n.X + 1):(n.X + n.theta)] <- .08
    scale[(n.X + n.theta + 1):(n.param)] <- 2.2
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
        X_n = matrix(state[1:n.X], nrow = N.l96, ncol = N + 1)
        theta = state[(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        Sigma = state[(n.X + n.theta + 1):n.param]
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1] = theta

        param_mat[i, 2:(n.sigma + n.theta)] = Sigma
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
    X = matrix(, nrow = N.l96, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

# hyper-parameters
to = 0 # initial time
tf = 10 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, N.l96), nrow = N.l96, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, N.l96) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha = 8 # changed later
beta = 2
a4 = 10 #2
b4 = (a4 - 1) * 0.5 #.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 #/ del_t
R = diag(.05, N.l96) # observational error
inv_R = diag(1 / (0.05), N.l96)
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
n.X = N.l96 * (N + 1)
n.theta = 1
n.sigma = N.l96
n.param = n.X + n.theta + n.sigma
n = 1e6


X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.5, N.l96)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
X = X[, 1:(N + 1)]
#load('../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, N.l96), sigma = R)) # observations from Lorenz-63
init = numeric(n.param)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- c(8) # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = .5 # inital \Sigma should also be positive semi definite

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
load('l96_1e5_cwise_spikes_interp_diffuse_init_theta_try')
alpha = mean(colMeans(to_save[[1]][[1]][5e3:1e4, 1:N.l96]))
init[(n.X + 1):(n.X + n.theta)] = alpha
init[(n.X + n.theta + 1):n.param] = colMeans(to_save[[1]][[1]][5e3:1e4, 69:72])

#scale[c(6007, 6010, 6012)] <- 100
ans = MH(n, init)
pm = ans[[1]]
colMeans(pm)

#pdf(file = 'RWMH_interp.pdf')
#plot.ts(pm[(n/4):n,])
#dev.off()

save(ans, file = "l96_interp_1e6_MH")
