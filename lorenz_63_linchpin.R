set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    


    # all the elements of theta should be positive
    if (min(theta) <= 0)
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
    p1 = p1 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o)

    #######################################################################
    #p1 = (sum(dmvnorm(t(Y - X_t), sigma = R, log = TRUE))
         #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
    ######################################################################

    p2 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3


    f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = t(diff(t(X_n)))
    beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
    p3 = 3 * lgamma(a4 + N/2) - (a4 + N/2) * sum(log(b4 + beta_tmp))

    return(p1 + p2 + p3)

}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = 6)
    scale = rep(0.001, n.X + n.theta)
    scale[(n.X + 1):(n.X + n.theta)] = 0.08
    accept.prob = 0

    for (i in 1:n) {
        chain = metrop(ludfun, init, 1, scale = scale)
        state = chain$batch
        accept.prob = accept.prob + chain$accept
        X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
        theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg  + state[1:n.X]
        param_mat[i,1:3] = theta

        Sigma = numeric(length = 3)
        f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
        del_X = t(diff(t(X_n)))
        beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
        Sigma[1] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
        Sigma[2] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
        Sigma[3] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])

        param_mat[i,4:6] = Sigma
        init = state
    }

    print(accept.prob / n)
    X_avg = X_avg / n
    final_output = list(param_mat, X_avg)
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
Nobs = 10 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(10, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha1 = 20 # Prior for \sigma is Gamma (alpha1, beta1)
alpha2 = 56 # Prior for \rho is Gamma (alpha2, beta2)
alpha3 = 6 # Prior for \beta is Gamma (alpha3, beta3)
beta1 = 0.5*2
beta2 = 0.5*2
beta3 = 0.5*2
a4 = 2*2
b4 = 6

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000/del_t
R = diag(2, 3) # observational error
inv_R = solve(R)
n.X = 3 * (N + 1)
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) + rnorm(n.X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1,c(10, 28, 8 / 3), sigma = diag(3,3)) # random initial values for MCMC

ans = linchpin(1e4, init)
pm = ans[[1]]
colMeans(pm)
plot.ts(pm)
#save(ans, file = "l63_linch_burnin_1e7")

#save(X, file = 'burninX')
#load('l63_linch_burnin')
#params = ans[[1]]
#load('burninX')
#Xe = ans[[2]]
#Xe = matrix(Xe, nrow = 3)
#sum((Xe - as.vector(X))^2)
#pdf("L63_sigma_linchpin_1e4.pdf", height = 6, width = 6)
#plot(density(params[,1]), ylab = expression(sigma), main = expression(paste('density plot ', sigma)))
#dev.off()

#pdf("L63_rho_linchpin_1e4.pdf", height = 6, width = 6)
#plot(density(params[, 2]), ylab = expression(rho), main = expression(paste('density plot ', rho)))
#dev.off()

#pdf("L63_beta_linchpin_1e4.pdf", height = 6, width = 6)
#plot(density(params[, 3]), ylab = expression(beta), main = expression(paste('density plot ', beta)))
#dev.off()

#pdf("L63_Sigma_z_linchpin_1e4.pdf", height = 6, width = 6)
#plot(density(params[, 6]), ylab = expression(Sigma[z]), main = expression(paste('density plot ', Sigma[z])))
#dev.off()

#pdf("L63_TS_linchpin_1e4.pdf", height = 6, width = 6)
#plot.ts(params, main = 'Time series plots')
#dev.off()

#pdf("L63_butterfly_xz_EM_1e4.pdf", height = 6, width = 6)
#plot(X[1,], X[3,], type = 'l', lwd = 2, xlab = 'X(t)', ylab = 'Z(t)', main = 'Euler-Maruyama')
#dev.off()

#pdf("L63_butterfly_xz_linchpin_1e4.pdf", height = 6, width = 6)
#plot(Xe[1,], Xe[3,], type = 'l', lwd = 2, xlab = 'X(t)', ylab = 'Z(t)', main = 'Linchpin')
#dev.off()