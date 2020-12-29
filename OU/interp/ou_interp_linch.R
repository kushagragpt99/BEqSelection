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

    p1 = (sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE)))
    p2 = dnorm(theta, mean = alpha, sd = beta, log = TRUE)

    f = sapply(X_n, drift_fun, theta) #sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, theta)
    del_X = diff(X_n)
    beta_tmp = sum((del_X / del_t - f[-(N + 1)]) ^ 2) * del_t / 2
    p3 = -(a4 + N / 2) * log(b4 + beta_tmp)

    return(p1 + p2 + p3)

}

linchpin <- function(n, init) {
    X_avg = numeric(length = n.X)
    param_mat = matrix(, nrow = n, ncol = n.theta + n.sigma)
    scale = rep(0.012, n.X + n.theta) # 0.015
    scale[(n.X + 1):(n.X + n.theta)] = 1.25 # 1.5
    accept.prob = 0
    state = init
    for (i in 1:n) {
        if (i %% (n / 10) == 0) print(c(i, accept.prob / i))

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

        param_mat[i, n.theta + n.sigma] = Sigma
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
    X = numeric(length = N + 1)
    X[1] = X0
    for (i in 2:(N + 1))
        X[i] = X[i - 1] + drift_fun(X[i - 1], theta) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

load('ou_linch_spike_interp_T2_1e5_try')
# hyper-parameters
to = attr$to # initial time
tf = attr$tf # final time
Nobs = attr$Nobs # no of observations (Y) per time step
del_t = attr$del_t # discrete approximation of dt
alpha = 2 # changed later
beta = 2
a4 = 10 #attr$a4
b4 = (a4-1) #attr$b4

K = attr$K # no of real life observations, i.e. size of Y
N = attr$N # no of discretizations of the Lorenz-63, i.e. size of X
seq.Y = attr$seq.Y
N = attr$N
burn_in = attr$burn_in
R = attr$R #diag(2, 3) # observational error
inv_R = attr$inv_R #solve(R)
n.X = attr$n.X
n.theta = 1
n.sigma = attr$n.sigma
n.param = attr$n.param
n = 1e6

alpha = -1 * mean(to_save[[1]][[1]][(attr$n / 2):attr$n, 1])

X_total = euler_maruyama(-0.5, del_t, burn_in + N, 2, 1) # generating sample from Lorenz-63
X = X_total[(burn_in):(N + burn_in)]
#load('../burninX')
Y = X[seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
init = numeric(n.X + n.theta)
init[(1:n.X)] <- as.numeric(X) #+ rnorm(n.X) #runif(n.param, 0, 5)
init[(n.X + 1):(n.X + n.theta)] <- alpha #rmvnorm(1, c(10, 28, 8 / 3), sigma = diag(0.5, 3)) # random initial values for MCMC

X.interp = X #matrix(-50,nrow = 3, ncol = N + 1)
y.index = 1
#X.interp[,1] = X[,1]
for (i in seq(2, N + 1, N / K)) {
    if (i == 2) {
        X.interp[2] = Y[1]
    } else {
        X.interp[(i - N / K + 1):i] = seq(Y[y.index], Y[y.index + 1], (Y[y.index + 1] - Y[y.index]) * K / N)[-1]
        y.index = y.index + 1
    }

}

init[(1:n.X)] <- as.numeric(X.interp)

ans = linchpin(n, init)
pm = ans[[1]]
colMeans(pm)
#pdf(file = 'inf_interp_tune.pdf')
#plot.ts(pm)
#dev.off()
save(ans, file = "OU_linch_interp_1e6_T2_fin")
# 11:30