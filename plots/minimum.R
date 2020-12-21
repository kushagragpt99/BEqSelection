set.seed(1)
library(mcmcse)
library(mvtnorm)
source('SimTools/R/supp.R')
source('SimTools/R/Smcmcclass.R')

## LORENZ 63
# TRUTH

# spike slab table
load('L63/truth/l63_linch_T_20_5e4_cwise_1_spikes_init_theta_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(matrix(gamma_est, nrow = 3))

# inference model density plots for parameters with gamma >= 0.5
load('L63/truth/l63_linch_1e6')
params = to_save[[1]][[1]]
n = dim(params)[1]
params = params[(n/4+1):n,]  # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(10, 28, 8 / 3, .6, .6, .6)
names = c(expression(sigma), expression(rho), expression(beta), expression(Sigma[x]), expression(Sigma[y]), expression(Sigma[z]))
for (i in 1:length(true_params)) {
    print(i)
    data = Smcmc(params[, i])
    pdf(file = paste('L63/truth/L63_truth_density',i,'.pdf',sep = '_'))
    plot(data, Q = c(0.05, 0.95),xlim = range(data,true_params[i]), main = names[i], xlab = " ")
    abline(v = true_params[i], col = 'red')
    dev.off()
}

## LORENZ 96
# TRUTH

## spike slab table
load('L96/truth/l96_5e4_cwise_spikes_truth_diffuse_try')
gamma_est = colMeans(to_save[[1]][[1]][, (attr$n.sigma + attr$n.theta + 1):(attr$n.sigma + 2 * attr$n.theta)])
print(matrix(gamma_est, nrow = 4))

# inference model density plots for parameters with gamma >= 0.5
load('L96/truth/l96_linch_1e6')
params = ans[[1]]
n = dim(params)[1]
params = params[(n / 4 + 1):n,] # remove n/4 burnin
selected_params = attr$mu_truth[which(gamma_est >= 0.5)]
true_params = c(8, 0.5, 0.5, 0.5, 0.5)
names = c(expression(theta), expression(Sigma[1]), expression(Sigma[2]), expression(Sigma[3]), expression(Sigma[4]))
for (i in 1:length(true_params)) {
    print(i)
    data = Smcmc(params[, i])
    pdf(file = paste('L96/truth/L96_truth_density', i, '.pdf', sep = '_'))
    plot(data, Q = c(0.05, 0.95), xlim = range(data, true_params[i]), main = names[i], xlab = " ")
    abline(v = true_params[i], col = 'red')
    dev.off()
}

