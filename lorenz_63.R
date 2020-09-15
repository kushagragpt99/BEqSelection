set.seed(1)
library(mvtnorm)
library(matrixcalc)
library(mcmc)

drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

ludfun <- function(state) {
    
    X_n = matrix(state[1:(3 * (N + 1))], nrow = 3, ncol = N + 1)
    Sigma_vec = state[(3 * N + 7):(3 * N + 12)]
    Sigma = matrix(c(Sigma_vec[1], Sigma_vec[2], Sigma_vec[3], Sigma_vec[2], Sigma_vec[4], Sigma_vec[5], Sigma_vec[3], Sigma_vec[5], Sigma_vec[6]), nrow = 3)
    theta = state[(3 * N + 4):(3 * N + 6)]
    if (min(theta) <= 0)
        return(-Inf)
    if (is.positive.semi.definite(Sigma) == FALSE)
        return(-Inf)

    X_t = X_n[, seq(1, N + 1, N / K)]
    inv_R = solve(R)

    p1 = 0
    for (k in 1:K) {
        Y.t = t(t(Y[, k]))
        X_t.t = t(t(X_t[,k]))
        p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
    }
    p1 = -0.5 * p1

    p2 = 0
    
    inv_Sig = solve(Sigma)
    for (k in 1:N) {
        del_X = matrix(X_n[, k + 1] - X_n[, k], nrow = 3, ncol = 1)
        f_k = drift_fun(X[, k], theta)
        #print(dim(del_X))
        #print(dim(f_k))
        p2 = p2 + t(del_X/del_t - f_k) %*% inv_Sig %*% (del_X/del_t - f_k)
    }
    p2 = -0.5 * p2
    p2 = p2 - 0.5 * t(t(t(X_n[,1])) - tau_o) %*% solve(lam_o) %*% (t(t(X_n[,1])) - tau_o) - (N / 2) * log(det(Sigma * del_t))

    
    p3 = (alpha1 - 1) * log(theta[1]) - theta[1] / beta1 + (alpha2 - 1) * log(theta[2]) - theta[2] / beta2 + (alpha3 - 1) * log(theta[3]) - theta[3] / beta3

    return(p1 + p2 + p3)

}


euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + det(del_t * Sigma) * rmvnorm(1, sigma = diag(1, 3))
    return(X)
}
# X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))

to = 0
tf = 20
Nobs = 10
del_t = 0.01
tau_o = matrix(rep(0, 3), nrow = 3, ncol = 1)
lam_o = diag(1, 3)
alpha1 = 20
alpha2 = 56
alpha3 = 6
beta1 = 0.5
beta2 = 0.5
beta3 = 0.5

K = (tf - to) * Nobs
N = (tf - to) / del_t
R = diag(2,3)

X = euler_maruyama(rmvnorm(1, tau_o, lam_o), del_t, N, c(10,28, 8/3), diag(6,3))
Y = X[, seq(1, N + 1, N / K)] + t(rmvnorm(K + 1, mean = rep(0, 3), sigma = R))
init = runif(3 * N + 12, 0, 5)
init[(3*N + 7):(3*N + 12)] = c(1,0,0,1,0,1)

chain = metrop(ludfun, init, nbatch = 1e3, scale = 0.2)

#R = diag(2, nrow = 3)
#inv_R = solve(R)
#fun <- function(arg1, arg2) {
    #X = t(t(arg2[1:3]))
    #Y = t(t(arg2[4:6]))
    #ans = arg1 + t(Y - X) %*% inv_R %*% (Y - X)
    #return(ans)
#}

#X1 <- matrix(1:3, nrow = 3)
#Y1 <- matrix(4:6, nrow = 3)
#X2 <- matrix(2:4, nrow = 3)
#Y2 <- matrix(5:7, nrow = 3)
#list1 = list(X1, X2)
#list2 = list(Y1, Y2)

#result.df <- expand.grid(list1, list2)
#result.list <- lapply(apply(result.df, 1, identity), unlist)
#listy = list(0)
#listy = append(listy, result.list)

#ans = Reduce(fun, listy)
#Z <- list(X, Y)
#Z[[3]] = X
#Z[[4]] = 0
#x <- matrix(1:10, ncol = 2)
#y <- x + 300
#z <- x + 50
#MATS <- list(x, y, z)
#Reduce(fun, t(MATS)%*%MATS)
#apply(MATS, 1, sum)

#sq = matrix(1:12, nrow = 3)