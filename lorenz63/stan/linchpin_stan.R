set.seed(1)
library(mvtnorm)
library(mcmc)
library(invgamma)
library(rstan)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = 6)
    scale = rep(0.001, n.X + n.theta)
    scale[(n.X + 1):(n.X + n.theta)] = 0.08
    accept.prob = 0

    for (i in 1:n) {
        #chain = metrop(ludfun, init, 1, scale = scale)
        #state = chain$batch
        #accept.prob = accept.prob + chain$accept
        state = p2[i,]
        X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
        theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
        X_avg = X_avg + state[1:n.X]
        param_mat[i, 1:3] = theta

        Sigma = numeric(length = 3)
        f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
        del_X = t(diff(t(X_n)))
        beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
        Sigma[1] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
        Sigma[2] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
        Sigma[3] = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])

        param_mat[i, 4:6] = Sigma
        #init = state
    }

    #print(accept.prob / n)
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
beta1 = 0.5 
beta2 = 0.5 
beta3 = 0.5 
a4 = 2 
b4 = 6

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 5000 / del_t
R = diag(2, 3) # observational error
inv_R = solve(R)
n.X = 3 * (N + 1)
n.theta = 3
n.sigma = 3
n.param = n.X + n.theta + n.sigma
seq_t = seq(2, N + 1, N / K)
n = 2e4
burn_in_n = n/2

#X_total = euler_maruyama(c(0,0,25), del_t, N + burn_in, c(10, 28, 8 / 3), diag(6, 3)) # generating sample from Lorenz-63
#X = X_total[, (burn_in):(N + burn_in)]
load('../burninX')
Y = X[, seq(2, N + 1, N / K)] + t(rmvnorm(K, mean = rep(0, 3), sigma = R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- X #as.numeric(X) + rnorm(n.X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1, c(10, 28, 8 / 3), sigma = diag(3, 3)) # random initial values for MCMC

initf <- function() {
    print('you shall not pass***************************************8')
    return(list(X = init[(1:n.X)], theta = init[(n.X + 1):(n.X + n.theta)]))
}
#initf <- function() {
    #print('you shall not pass***************************************8')
    #return(list(X_n = matrix(init[(1:n.X)],nrow = 3), theta = init[(n.X + 1):(n.X + n.theta)]))
#}
model = stan_model('linchpin.stan')

fit <- sampling(model, list(N = N, K = K, n_X = n.X, n_theta = n.theta, n_sigma = n.sigma, y = Y, seq_t = seq_t, inv_R = inv_R,
                inv_lam_0 = inv.lam_o, tau_0 = tau_o[, 1], del_t = del_t, a1 = alpha1, a2 = alpha2, a3 = alpha3, b1 = beta1,
                b2 = beta2, b3 = beta3, a4 = a4, b4 = b4), iter = n, warmup = burn_in_n, chains = 1, init = initf,
                control = list(max_treedepth = 4))

chain_info = capture.output(cat("no of samples from MC is ", n, " \n using warmup ", burn_in_n,
                 "max tree depth is ", 4, " \n starting from init ", "\n priors centered at truth",
                 " time period ", tf))

print(chain_info)

p2 = extract(fit, inc_warmup = FALSE, permuted = FALSE)
p2 = p2[,1,1:(n.X+n.theta)]
#p2 = p2[,1,1:(n.X+n.theta)]
#to_save = list(fit, chain_info)
#save(to_save, file = "nuts_linchpin_td_4") ######not 5, 3
#p2 = extract(to_save[[1]], inc_warmup = TRUE)
#colMeans(p1[, 1,])

#al = c(as.numeric(X) + rnorm(n.X, sd = 0.1), colMeans(p2$theta))
ans = linchpin(n/2, init)
colMeans(ans[[1]])
to_save = list(ans, chain_info)
save(to_save, file = "nuts_linchpin_td_4")

#p2.1 = matrix(0,nrow = 500, ncol = n.X+n.theta)
#for (i in 1:500) {
    #p2.1[i,] = c(as.numeric(al[i,,]),10,28,3)
#}

