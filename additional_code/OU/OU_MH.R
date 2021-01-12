set.seed(1)
library(mvtnorm)
#library(matrixcalc)
library(mcmc)
library(invgamma)

# drifet function for Lorenz-63
drift_fun <- function(X, theta) {
    ans = -theta*X
    return(ans)
}

# logarithm of the unnormalized posterior
ludfun <- function(state) {
    # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
    # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.

    X_n = state[1:n.X]
    theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
    #Sigma_vec = state[(3 * N + 7):(3 * N + 12)]
    Sigma = state[(n.X + n.theta + 1):n.param]

    # all the elements of theta should be positive
    if (min(theta) <= 0)
        return(-Inf)
    # \Sigma should be positive semi-definite
    if (min(diag(Sigma)) < 0)
        return(-Inf)

    # Extracting observed data
    X_t = X_n[seq(2, N + 1, N / K)]


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
    p1 = sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE))
    ######################################################################

    # p2 is the log of prior of X conditional on theta
    p2 = 0
    inv_Sig = 1/(Sigma)
    for (k in 1:N) {
        del_X = X_n[k + 1] - X_n[k] # try using diff() function
        f_k = drift_fun(X_n[k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + (del_X / del_t - f_k) * inv_Sig * (del_X / del_t - f_k)
    }
    p2 = -0.5 * p2 * del_t

    ########################################################################
    #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    #del_X = t(diff(t(X_n)))
    #p2 = sum(dmvnorm(t(del_X - f[, - (N + 1)] * del_t), sigma = Sigma * del_t, log = TRUE))
    ########################################################################

    # store inv.lam_o globally
    p2 = p2 - (N / 2) * log(del_t*Sigma) #determinant(Sigma * del_t, logarithm = TRUE)$modulus

    # p3 is the log of priors of theta
    p3 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 #+ (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    ## add inverse gamma
    p4 <- sum(dinvgamma(diag(Sigma), shape = 2, scale = 1 , log = TRUE))
    return(p1 + p2 + p3 + p4)

}

# Numerical method to sample from SDE
euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = numeric(length = N+1)
    X[1] = X0
    for (i in 2:(N + 1))
        X[i] = X[i - 1] + t(drift_fun(X[i - 1], theta)) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

# hyper-parameters
to = 0 # initial time
tf = 20 # final time
Nobs = 2 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1) # prior mean for X[0], i.e. initial state of Lorenz-63 oricess
lam_o = diag(1, 3) # prior covariance matrix of X[0]
inv.lam_o = solve(lam_o)
alpha1 = 4 # Prior for \sigma is Gamma (alpha1, beta1)
alpha2 = 56 # Prior for \rho is Gamma (alpha2, beta2)
alpha3 = 6 # Prior for \beta is Gamma (alpha3, beta3)
beta1 = 0.5
beta2 = 0.5
beta3 = 0.5

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
R = 0.04 #diag(2, 3) # observational error
inv_R = 1/(0.04) #solve(R)
n.X = (N + 1)
n.theta = 1
n.sigma = 1
n.param = n.X + n.theta + n.sigma


X = euler_maruyama(-0.5, del_t, N, 2, 1) # generating sample from Lorenz-63
Y = X[seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
init = numeric(n.param)
init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- 2 # random initial values for MCMC
init[(n.X + n.theta + 1):(n.param)] = 1 # inital \Sigma should also be positive semi definite

scale <- rep(.001, n.param)
scale[(n.X + 1):(n.X + n.theta)] <- .05
scale[(n.X + n.theta + 1):(n.param)] <- .2
#scale[c(6007, 6010, 6012)] <- 100
chain = metrop(ludfun, init, nbatch = 5e3, scale = scale) # running MH

print(chain$accept)
out <- chain$batch[, (n.X + 1):n.param]
plot.ts(out)
print(colMeans(out))

#save(chain, file = "OU_1e5_MH")
