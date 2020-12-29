set.seed(1)
library(mcmcse)
library(mvtnorm)
library(invgamma)
library(scales)


## LORENZ 63
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}

# TRUTH

# spike slab table
load('L63/truth/l63_linch_T_20_5e4_cwise_1_spikes_init_theta_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(matrix(gamma_est, nrow = 3))

col = as.numeric(attr$mu_truth != 0)
pdf(file = 'L63/truth/L63_truth_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n",xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend('right',legend = c("1 in L63","0 in L63"), col = c(col[1]+1, col[3]+1), pch = c(16,16), bty = 'n')
dev.off()

# inference model density plots for parameters with gamma >= 0.5
load('L63/truth/l63_linch_2e6_fin_var')
params = to_save[[1]][[1]]
n = dim(params)[1]
params = params[(n/4+1):n,]  # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(10, 28, 8 / 3, .6, .6, .6)
names = c(expression(sigma), expression(rho), expression(beta), expression(Sigma[x]), expression(Sigma[y]), expression(Sigma[z]))
for (i in 1:length(true_params)) {
    #print(i)
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('L63/truth/L63_truth_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}
Xtruth = to_save[[1]][[2]]
# trajectory
#X_total = euler_maruyama(c(0, 0, 25), del_t = 0.01, attr$N + attr$burn_in, c(10, 28, 8 / 3), diag(.6, 3))
#X = X_total[,attr$burn_in:(attr$burn_in+attr$N)]
#X.mat = matrix(to_save[[1]][[2]], nrow = 3)
#pdf(file = 'L63/truth/L63_truth_traj.pdf')
#plot(X[1,], X[3,], type = 'l', xlab = 'X', ylab = 'Z', lwd = 1.5)
#lines(X.mat[1,], X.mat[3,], type = 'l', col = 'red', lwd = 1.5)
#dev.off()

# INTERP

 #spike slab table

load('L63/interp/l63_linch_T_20_1e5_cwise_1_spikes_interp_diffuse_6_by_10_scale_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(signif(matrix(gamma_est, nrow = 3)),2)
col = as.numeric(attr$mu_truth != 0)
pdf(file = 'L63/interp/L63_interp_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n", xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend('right', legend = c("1 in L63", "0 in L63"), col = c(col[1] + 1, col[3] + 1), pch = c(16, 16), bty = 'n')
dev.off()

# inference model density plots for parameters with gamma >= 0.5
load('L63/interp/l63_interp_5e6')
params = ans[[1]]
n = dim(params)[1]
params = params[(4*n / 5 + 1):n,] # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(10, 28, 8 / 3, .6, .6, .6)
names = c(expression(sigma), expression(rho), expression(beta), expression(Sigma[x]), expression(Sigma[y]), expression(Sigma[z]))
for (i in 1:length(true_params)) {
    #print(i)
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('L63/interp/L63_interp_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}
Xinterp = ans[[2]]

# trajectory 
X_total = euler_maruyama(c(0, 0, 25), del_t = 0.01, attr$N + attr$burn_in, c(10, 28, 8 / 3), diag(.6, 3))
X = X_total[, attr$burn_in:(attr$burn_in + attr$N)]
X.mat.true = matrix(Xtruth, nrow = 3)
X.mat.interp = matrix(Xinterp, nrow = 3)
pdf(file = 'L63/L63_traj.pdf')
plot(X[1,], X[3,], type = 'l', xlab = 'X', ylab = 'Z')
#lines(X.mat.true[1,], X.mat.true[3,], type = 'l', col = 'red')
lines(X.mat.interp[1,], X.mat.interp[3,], type = 'l', col = 'red')
#legend('bottomright', legend = c('truth', 'X init truth', 'X init interpolation'), lty = c(1,1,1), col = c('black', 'red', 'blue'),bty = 'n')
legend('bottomright', legend = c('truth','estimate'), lty = c(1, 1), col = c('black', 'red'), bty = 'n')
dev.off()

#X.mat = matrix(to_save[[1]][[2]], nrow = 3)
#pdf(file = 'L63/interp/L63_interp_traj.pdf')
#plot(X[1,], X[3,], type = 'l', xlab = 'X', ylab = 'Z', lwd = 1.5)
#lines(X.mat[1,], X.mat[3,], type = 'l', col = 'red', lwd = 1.5)
#dev.off()

# Chaoticity plots

mu_est = colMeans(params[, 1:3])
X_total = euler_maruyama(c(0, 0, 25), del_t = 0.01, attr$N + attr$burn_in, mu_est, diag(.6, 3))
X.mu_est = X_total[, attr$burn_in:(attr$burn_in + attr$N)]

pdf('L63/l63_mu_est.pdf')
plot(X[1,], X[3,], type = 'l', xlab = 'X', ylab = 'Z', ylim = range(X[3,], X.mu_est[3,]), xlim = range(X[1,], X.mu_est[1,]))
lines(X.mu_est[1,], X.mu_est[3,], type = 'l', col = 'red')
op <- par(cex = 1.2)
legend('bottomright', legend = c('truth', 'estimate'), lty = c(1,1), col = c('black', 'red'), bty = 'n')
dev.off()


## Butterfly
# drift function for Lorenz-63
drift_fun <- function(X, mu_truth) {
    ans = c(mu_truth[1] * (X[2] - X[1]), mu_truth[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - mu_truth[3] * X[3])
    return(t(t(ans)))
}

euler_maruyama <- function(X0, del_t, N, mu_truth, Sigma) {
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], mu_truth)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}

# hyper-parameters
to = 0 # initial time
tf = 20 # final time
Nobs = 20 # no of observations (Y) per time step
del_t = 0.01 # discrete approximation of dt
a4 = 10 #2
b4 = (a4 - 1) * 0.6 #.6
mu_truth = c(10, 28, 8 / 3)

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
seq.Y = seq(2, N + 1, N / K)
N = tail(seq.Y, 1)
burn_in = 5000 / del_t
n.X = 3 * (N + 1)
R = diag(.05, 3) # observational error
inv_R = diag(20, 3)

X_total = euler_maruyama(c(0, 0, 25), del_t, N + burn_in, mu_truth, diag(.06, 3)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]
X = X[, 1:(N + 1)]

Sigma = numeric(length = 3)
f = sapply(split(X, rep(1:ncol(X), each = nrow(X))), drift_fun, mu_truth)
del_X = t(diff(t(X)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,3)

X_1 = X + rnorm(n.X, sd = 1)
X_.5 = X + rnorm(n.X, sd = .5)
X_.1 = X + rnorm(n.X, sd = .1)
X_.05 = X + rnorm(n.X, sd = .05)
X_.01 = X + rnorm(n.X, sd = .01)

pdf(file = 'L63/butterfly_sd1.pdf', height = 5, width = 5)
plot(X[1,], X[3,], type = 'l', lwd = 2)
lines(X_1[1,], X_1[3,], type = 'l', col = alpha(2,0.8))
dev.off()

Sigma = numeric(length = 3)
f = sapply(split(X_1, rep(1:ncol(X_1), each = nrow(X_1))), drift_fun, mu_truth)
del_X = t(diff(t(X_1)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,3)

load('L63/l63_jitter_1e5_sd_1')
round(colMeans(to_save[[1]][[1]]),2)

pdf(file = 'L63/butterfly_sd_5_by_10.pdf', height = 5, width = 5)
plot(X[1,], X[3,], type = 'l', lwd = 2, xlab = 'X', ylab = 'Z')
lines(X_.5[1,], X_.5[3,], type = 'l', col=alpha(3, 0.8))
dev.off()

Sigma = numeric(length = 3)
f = sapply(split(X_.5, rep(1:ncol(X_.5), each = nrow(X_.5))), drift_fun, mu_truth)
del_X = t(diff(t(X_.5)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,2)

load('L63/l63_jitter_1e5_sd_5_by_10')
round(colMeans(to_save[[1]][[1]]),2)

pdf(file = 'L63/butterfly_sd_1_by_10.pdf', height = 5, width = 5)
plot(X[1,], X[3,], type = 'l', lwd = 2, xlab = 'X', ylab = 'Z')
lines(X_.1[1,], X_.1[3,], type = 'l', col=alpha(4, 0.8))
dev.off()

Sigma = numeric(length = 3)
f = sapply(split(X_.1, rep(1:ncol(X_.1), each = nrow(X_.1))), drift_fun, mu_truth)
del_X = t(diff(t(X_.1)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,2)

load('L63/l63_jitter_1e5_sd_1_by_10')
round(colMeans(to_save[[1]][[1]]),2)

pdf(file = 'L63/butterfly_sd_5_by_100.pdf', height = 5, width = 5)
plot(X[1,], X[3,], type = 'l', lwd = 2, xlab = 'X', ylab = 'Z')
lines(X_.05[1,], X_.05[3,], type = 'l', col=alpha(5, 0.8))
dev.off()

Sigma = numeric(length = 3)
f = sapply(split(X_.05, rep(1:ncol(X_.05), each = nrow(X_.05))), drift_fun, mu_truth)
del_X = t(diff(t(X_.05)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,2)

load('L63/l63_jitter_1e5_sd_5_by_100')
round(colMeans(to_save[[1]][[1]]),2)

pdf(file = 'L63/butterfly_sd_1_by_100.pdf', height = 5, width = 5)
plot(X[1,], X[3,], type = 'l', lwd = 2, xlab = 'X', ylab = 'Z')
lines(X_.01[1,], X_.01[3,], type = 'l', col = alpha(6, 0.8))
dev.off()

Sigma = numeric(length = 3)
f = sapply(split(X_.01, rep(1:ncol(X_.01), each = nrow(X_.01))), drift_fun, mu_truth)
del_X = t(diff(t(X_.01)))
beta_tmp = rowSums((del_X / del_t - f[, - (N + 1)]) ^ 2) * del_t / 2
Sigma[1] = (b4 + beta_tmp[1]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[1])
Sigma[2] = (b4 + beta_tmp[2]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[2])
Sigma[3] = (b4 + beta_tmp[3]) / (a4 + N / 2 - 1) # rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp[3])
round(Sigma,3)

load('L63/l63_jitter_1e5_sd_1_by_100')
round(colMeans(to_save[[1]][[1]]), 2)


# lag plot
drift_fun <- function(X, theta) {
    ans = c(theta[1] * (X[2] - X[1]), theta[2] * X[1] - X[2] - X[1] * X[3], X[1] * X[2] - theta[3] * X[3])
    return(t(t(ans)))
}

euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    set.seed(1)
    X = matrix(, nrow = 3, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}

to = 0 # initial time
tf = 10 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 0 #5000 #/ del_t

X = euler_maruyama(c(0,0,25), del_t, N, c(10, 28, 8 / 3), diag(.06, 3)) # generating sample from Lorenz-63

Xtmp = euler_maruyama(c(0, 0, 25), del_t, N, c(10, 28, 8 / 3), diag(1, 3)) # generating sample from Lorenz-63

name = c(expression(X[1]), expression(X[2]), expression(X[3]))
for (i in 1:3) {
    pdf(paste('L63/intro_lag', i, '.pdf', sep = '_'), height = 3, width = 5)
    plot(seq(0, tf, del_t), X[i,], type = 'l', ylim = range(X[i,], Xtmp[i,]), ylab = name[i], xlab = 'Time')
    lines(seq(0, tf, del_t), Xtmp[i,], type = 'l', col = 'red')
    dev.off()

}

### LORENZ 96

drift_fun <- function(X, theta) {
    ans = matrix(, nrow = attr$N.l96, ncol = 1)
    for (i in 0:(attr$N.l96 - 1)) {
        ans[i + 1, 1] = (X[(i + 1) %% attr$N.l96 + 1] - X[(i - 2) %% attr$N.l96 + 1]) * X[(i - 1) %% attr$N.l96 + 1] - X[i + 1] + theta
    }
    return(ans)
}

euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = matrix(, nrow = attr$N.l96, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}

## TRUTH


## spike slab table
load('L96/truth/l96_5e4_cwise_spikes_truth_diffuse_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(signif(matrix(gamma_est, nrow = 4)),2)

col = as.numeric(attr$mu_truth != 0)
pdf(file = 'L96/truth/L96_truth_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n", xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend('right', legend = c("1 in L96", "0 in L96"), col = c(col[1] + 1, col[6] + 1), pch = c(16, 16), bty = 'n')
dev.off()

# inference model density plots for parameters with gamma >= 0.5
load('L96/truth/l96_linch_2e6_fin_var')
params = ans[[1]]
n = dim(params)[1]
params = params[(n / 4 + 1):n,] # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(8, 0.5, 0.5, 0.5, 0.5)
names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
for (i in 1:length(true_params)) {
    #print(i)
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('L96/truth/L96_truth_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}
Xtruth = ans[[2]]
# Trajectory
#X_total = euler_maruyama(rep(0, attr$N.l96), attr$del_t, attr$N + attr$burn_in, 8, diag(.5, attr$N.l96)) # generating sample from Lorenz-63
#X = X_total[, (attr$burn_in):(attr$N + attr$burn_in)]
#X.mat = matrix(ans[[2]], nrow = 4)
#names = c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]))
#for (i in 1:4) {
    #pdf(file = paste('L96/truth/L96_truth_traj',i,'.pdf',sep = '_'))
    #plot.ts(X[i,], ylab = names[i], ylim = range(X[i,], X.mat[i,]))
    #lines(X.mat[i,], type = 'l', col = 'red')
    #dev.off()
#}



# INTERP

## spike slab table
load('L96/interp/l96_1e5_cwise_spikes_interp_diffuse_init_theta_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(signif(matrix(gamma_est, nrow = 4)), 2)

col = as.numeric(attr$mu_truth != 0)
pdf(file = 'L96/interp/L96_interp_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n", xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend('right', legend = c("1 in L96", "0 in L96"), col = c(col[1] + 1, col[6] + 1), pch = c(16, 16), bty = 'n')
dev.off()

load('L96/interp/l96_linch_5e6_interp')
params = ans[[1]]
n = dim(params)[1]
params = params[(4*n / 5 + 1):n,] # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(8, 0.5, 0.5, 0.5, 0.5)
names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
for (i in 1:length(true_params)) {
    #print(i)
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('L96/interp/L96_interp_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}
Xinterp = ans[[2]]

# Trajectory

#X.mat = matrix(ans[[2]], nrow = 4)
#names = c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]))
#for (i in 1:4) {
    #pdf(file = paste('L96/interp/L96_interp_traj', i, '.pdf', sep = '_'))
    #plot.ts(X[i,], ylab = names[i], ylim = range(X[i,], X.mat[i,]))
    #lines(X.mat[i,], type = 'l', col = 'red')
    #dev.off()
#}

X_total = euler_maruyama(rep(0, attr$N.l96), attr$del_t, attr$N + attr$burn_in, 8, diag(.5, attr$N.l96)) # generating sample from Lorenz-63
X = X_total[, (attr$burn_in):(attr$N + attr$burn_in)]
X.mat.true = matrix(Xtruth, nrow = 4)
X.mat.interp = matrix(Xinterp, nrow = 4)
names = c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]))
for (i in 1:4) {
    pdf(file = paste('L96/L96_traj', i, '.pdf', sep = '_'))
    plot.ts(X[i,], ylab = names[i], ylim = range(X[i,], X.mat.true[i,], X.mat.interp[i,]))
    #lines(X.mat.true[i,], type = 'l', col = 'red')
    lines(X.mat.interp[i,], type = 'l', col = 'red')
    op <- par(cex = 1.2)
    #legend('topleft', legend = c('truth', 'X init truth', 'X init interpolation'), lty = c(1,1,1), col = c('black', 'red', 'blue'), bty = 'n')
    legend('topleft', legend = c('truth', 'estimate'), lty = c(1, 1), col = c('black', 'red'), bty = 'n')

    dev.off()
}

#pdf(file = '10th.pdf')
#plot.ts(ans[[1]][seq(1, dim(ans[[1]])[1], 50),])
#dev.off()

# Intro plot

drift_fun <- function(X, theta) {
    ans = matrix(, nrow = N.l96, ncol = 1)
    for (i in 0:(N.l96 - 1)) {
        ans[i + 1, 1] = (X[(i + 1) %% N.l96 + 1] - X[(i - 2) %% N.l96 + 1]) * X[(i - 1) %% N.l96 + 1] - X[i + 1] + theta
    }
    return(ans)
}

euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    set.seed(1)
    X = matrix(, nrow = N.l96, ncol = N + 1)
    X[, 1] = X0
    for (i in 2:(N + 1))
        X[, i] = X[, i - 1] + t(drift_fun(X[, i - 1], theta)) * del_t + rmvnorm(1, sigma = del_t * Sigma)
    return(X)
}


to = 0 # initial time
tf = 2 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
seq.Xdots = seq(1,N+1,0.5*N/K)
seq.Y = seq.Xdots[seq(1, length(seq.Xdots), 0.6* N / K)]

burn_in = 5000 #/ del_t
R = diag(2, N.l96) # observational error

X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.05, N.l96)) # generating sample from Lorenz-63
X = X_total[, (burn_in):(N + burn_in)]

X_total = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(0.05, N.l96)) # generating sample from Lorenz-63
Xtmp = X_total[, (burn_in):(N + burn_in)]
Xdots = X[,seq.Xdots]
Y = Xtmp[, seq.Y] + rnorm(seq.Y, sd = sqrt(R[1,1]))
name = c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]))
for (i in 1:4) {
    pdf(paste('L96/intro',i,'.pdf',sep = '_'))
    plot(seq(0,2,.01), X[i,], type = 'l', ylim = range(X, Xdots, Y), xlab = 'Time', ylab = name[i])
    points(seq.Xdots/100, Xdots[i,], col = 'red',pch = 1, cex = 0.5)
    points(seq.Y/100, Y[i,], col = 'blue', pch = 16)
    #segments(x0 = seq.Y/100, y0 = Xtmp[i, seq.Y], x1 = seq.Y/100, y1 = Y[i,], col = 'red')
    segments(x0 = seq.Y / 100, y0 = X[i, seq.Y], x1 = seq.Y / 100, y1 = Y[i,], col = 'blue')
    dev.off()
}


# lag intro plot

to = 0 # initial time
tf = 10 # final time
Nobs = 20 # no of observations (Y) per time step
N.l96 = 4
del_t = 0.01 # discrete approximation of dt

K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
burn_in = 0 #5000 #/ del_t

X = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(.05, N.l96)) # generating sample from Lorenz-63

Xtmp = euler_maruyama(rep(0, N.l96), del_t, N + burn_in, 8, diag(1, N.l96)) # generating sample from Lorenz-63

name = c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]))
for (i in 1:4) {
    pdf(paste('L96/intro_lag', i, '.pdf', sep = '_'), height = 3, width = 5)
    plot(seq(0,10,0.01),X[i,], type = 'l', ylim = range(X[i,], Xtmp[i,]), ylab = name[i], xlab = 'Time')
    lines(seq(0, 10, 0.01), Xtmp[i,], type = 'l', col = 'red')
    dev.off()

}

# RWMH vs linchpin
load('L96/truth/l96_1e5_MH')
param_MH = ans[[1]]
load('L96/truth/l96_linch_1e5_fin_var')
param_linch = ans[[1]]
linch = acf(param_linch[, 1], lag.max = 50, plot = FALSE)
pdf(file = 'L96/l96_rwmh_linch.pdf')
acf(param_MH[, 1], lag.max = 50, main = ' ')
lines(linch$acf, type = 'l',col = 'red', lwd = 2)
legend('topright', legend = c('RWMH', 'linchpin'), col = c('black', 'red'), lty = c(1, 1), bty = 'n')
dev.off()
ess(param_MH[,1])
ess(param_linch[,1])

## OU

drift_fun <- function(X, theta) {
    ans = -theta * X
    return(ans)
}

euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    X = numeric(length = N + 1)
    X[1] = X0
    for (i in 2:(N + 1))
        X[i] = X[i - 1] + drift_fun(X[i - 1], theta) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    return(X)
}
# TRUTH

# spike slab table
load('OU/truth/ou_linch_spike_truth_T2_1e5_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(gamma_est)

col = as.numeric(attr$mu_truth != 0)
pdf(file = 'OU/truth/OU_truth_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n", xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend(x = 3.3, y = 0.63, legend = c("1 in L96", "0 in L96"), col = c(col[1] + 1, col[2] + 1), pch = c(16, 16), bty = 'n')
dev.off()

#inference model density plots for parameters with gamma >= 0.5

#load('OU/truth/OU_linch_1e6_T2')
#params = ans[[1]]
#n = dim(params)[1]
#params = params[(n / 4 + 1):n,] # remove n/4 burnin
##selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
#true_params = c(2,1)
#names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
#for (i in 1:length(true_params)) {
    ##print(i)
    #data = params[, i]
    #qtl = quantile(data, probs = c(0.025, 0.975))
    #pdf(file = paste('OU/truth/OU_truth_density_infin', i, '.pdf', sep = '_'), height = 5, width = 5)
    #op <- par(cex.main = 1.5)
    #plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    #abline(v = qtl[1], lty = 2)
    #abline(v = qtl[2], lty = 2)
    #abline(v = true_params[i], col = "red")
    #abline(v = mean(data), col = "black")
    #dev.off()
#}

load('OU/truth/OU_linch_1e6_T2_fin')
params = ans[[1]]
n = dim(params)[1]
params = params[(n / 4 + 1):n,] # remove n/4 burnin
#selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(2, 1)
names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
for (i in 1:length(true_params)) {
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('OU/truth/OU_truth_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}

Xtruth = ans[[2]]
# trajectory
#X_total = euler_maruyama(-0.5, attr$del_t, attr$burn_in + attr$N, 2, 1) # generating sample from Lorenz-63
#X = X_total[(attr$burn_in):(attr$N + attr$burn_in)]
#pdf(file = 'OU/truth/OU_truth_traj.pdf')
#plot.ts(X, ylab = 'X', ylim = range(X, ans[[2]]))
#lines(ans[[2]], type = 'l', col = 'red')
#dev.off()

## INTERP

# spike slab table
load('OU/interp/ou_linch_spike_interp_T2_1e5_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(gamma_est)

col = as.numeric(attr$mu_truth != 0)
pdf(file = 'OU/interp/OU_interp_gamma.pdf', height = 4, width = 5)
plot(1:length(gamma_est), gamma_est, type = "n", xlab = 'index', ylab = expression(gamma))
for (i in 1:length(gamma_est)) {
    points(i, gamma_est[i], col = col[i] + 1, pch = 16)
}
abline(h = 0.5)
legend(x=3.3,y= 0.63, legend = c("1 in L96", "0 in L96"), col = c(col[1] + 1, col[2] + 1), pch = c(16, 16), bty = 'n')
dev.off()

#load('OU/interp/OU_linch_interp_1e6_T2')
#params = ans[[1]]
#n = dim(params)[1]
#params = params[(n / 4 + 1):n,] # remove n/4 burnin
##selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
#true_params = c(2, 1)
#names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
#for (i in 1:length(true_params)) {
    #data = params[, i]
    #qtl = quantile(data, probs = c(0.025, 0.975))
    #pdf(file = paste('OU/interp/OU_interp_density_infin', i, '.pdf', sep = '_'), height = 5, width = 5)
    #op <- par(cex.main = 1.5)
    #plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    #abline(v = qtl[1], lty = 2)
    #abline(v = qtl[2], lty = 2)
    #abline(v = true_params[i], col = "red")
    #abline(v = mean(data), col = "black")
    #dev.off()
#}

load('OU/interp/OU_linch_interp_1e6_T2_fin')
params = ans[[1]]
n = dim(params)[1]
params = params[(n / 4 + 1):n,] # remove n/4 burnin
#selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(2, 1)
names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
for (i in 1:length(true_params)) {
    data = params[, i]
    qtl = quantile(data, probs = c(0.025, 0.975))
    pdf(file = paste('OU/interp/OU_interp_density_fin', i, '.pdf', sep = '_'), height = 5, width = 5)
    op <- par(cex.main = 1.5)
    plot(density(data), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = qtl[1], lty = 2)
    abline(v = qtl[2], lty = 2)
    abline(v = true_params[i], col = "red")
    abline(v = mean(data), col = "black")
    dev.off()
}

Xinterp = ans[[2]]
# Trajectory
X_total = euler_maruyama(-0.5, attr$del_t, attr$burn_in + attr$N, 2, 1) # generating sample from Lorenz-63
X = X_total[(attr$burn_in):(attr$N + attr$burn_in)]
pdf(file = 'OU/OU_traj.pdf')
plot.ts(X, ylab = 'X', ylim = range(X, Xtruth, Xinterp))
#lines(Xtruth, type = 'l', col = 'red')
lines(Xinterp, type = 'l', col = 'red')
#legend('bottomright', legend = c('truth', 'X init truth', 'X init interpolation'), col = c('black', 'red', 'black'), lty = c(1,1,1), bty = 'n')
legend('bottomright', legend = c('truth', 'estimate'), lty = c(1, 1), col = c('black', 'red'), bty = 'n')
dev.off()

#pdf(file = 'OU/interp/OU_interp_traj.pdf')
#plot.ts(X, ylab = 'X', ylim = range(X, ans[[2]]))
#lines(ans[[2]], type = 'l', col = 'red')
#dev.off()

## Intro
#drift_fun <- function(X, theta) {
    #ans = -theta * X
    #return(ans)
#}
#euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
    #set.seed(1)
    #X = numeric(length = N + 1)
    #X[1] = X0
    #for (i in 2:(N + 1))
        #X[i] = X[i - 1] + drift_fun(X[i - 1], theta) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
    #return(X)
#}

#to = 0 # initial time
#tf = 2 # final time
#Nobs = 20 # no of observations (Y) per time step
#N.l96 = 4
#del_t = 0.01 # discrete approximation of dt

#K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
#N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
#seq.Xdots = seq(1, N + 1, 0.35 * N / K)
#seq.Y = seq.Xdots[seq(1, length(seq.Xdots), 1 * N / K)]

#burn_in = 5000 #/ del_t
#R = 0.05 # observational error

#X_total = euler_maruyama(-0.5, del_t, burn_in + N, 2, .1) # generating sample from Lorenz-63
#XOU = X_total[(burn_in):(N + burn_in)]

#X_total = euler_maruyama(-0.5, del_t, burn_in + N, 2, .3) # generating sample from Lorenz-63
#XtmpOU = X_total[(burn_in):(N + burn_in)]
#Xdots = XtmpOU[seq.Xdots]
#Y = XtmpOU[seq.Y] + rnorm(seq.Y, sd = sqrt(R))

#plot(XOU, type = 'l', ylim = range(XOU, XtmpOU))
#lines(XtmpOU, type = 'l', col = 'red')

#plot(XOU, type = 'l', ylim = range(XOU, Xdots, Y))
#points(seq.Xdots, Xdots, col = 'red')
#points(seq.Y, Y, col = 'blue', pch = 16)
##segments(x0 = seq.Y, y0 = XtmpOU[seq.Y], x1 = seq.Y, y1 = Y, col = 'red')
#segments(x0 = seq.Y, y0 = XOU[seq.Y], x1 = seq.Y, y1 = Y, col = 'blue')