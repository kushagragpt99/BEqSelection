    set.seed(1)
    library(mvtnorm)
    library(mcmc)
    library(invgamma)

    make_tilde <- function(X, t) {
        X_vec = c(X, X^2, t, t ^ 2)
        return(X_vec)
    }
    # drifet function for Lorenz-63
    drift_fun <- function(X, t, B) {
        #print(make_tilde(X,t))
        tildeX = matrix(make_tilde(X, t), nrow = 4, ncol = 1)
        B_mat = matrix(B, nrow = 1)
        #print(B)
        #print(dim(tildeX))
        ans = B_mat %*% tildeX
        return(ans)
    }

    drift_fun_true <- function(X, theta) {
        ans = -theta*X
        return(ans)
    }

    ludfun <- function(state, gamma) {
        # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
        # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
        #if (index == 0) {
        ##print('0')
        #all[1:n.X] = state
        #} else {
        ##print(index)
        #all[n.X+index] = state
        #}

        #X_n = matrix(all[1:n.X], nrow = 3, ncol = N + 1)
        #B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
        #B_mat = matrix(B_vec, nrow = 3)

        X_n = state[1:n.X]
        B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
        B_mat = matrix(B_vec, nrow = 1)

        # all the elements of theta should be positive
        #if (min(theta) <= 0)
        #return(-Inf)

        # Extracting observed data
        X_t = X_n[seq(2, N + 1, N / K)]


        # pi is the log of likelihood
        # This doesn't need a loop
        #p1 = 0
        ##print(dim(Y))
        #for (k in 1:K) {
            #Y.t = t(t(Y[, k]))
            #X_t.t = t(t(X_t[, k]))
            #p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
        #}
        #p1 = -0.5 * p1
        #p1 = p1 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o)

        #######################################################################
        p1 = sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE))
        #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
        ######################################################################
        B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
        p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
        #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

        f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(B_vec))
        #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
        del_X = diff(X_n)
        beta_tmp = sum((del_X / del_t - f[- (N + 1)]) ^ 2) * del_t / 2
        p3 = -(a4 + N / 2) * log(b4 + beta_tmp)

        return(p1 + p2 + p3)

    }

    ludfun.X <- function(state, gamma, all) {
        # State is the vector storing the vectors of length 3*N + 12. The first 3*(N+1) terms are Xs. The next three terms are the parameters \sigma, \rho & 
        # \beta. The remaining 6 terms are the \Sigma matrix. Definition of Sigma below shows how the symmetric matrix is constructed.
        #if (index == 0) {
        ##print('0')
        #all[1:n.X] = state
        #} else {
        ##print(index)
        #all[n.X+index] = state
        #}
        all[1:n.X] = state
        X_n = all[1:n.X]
        B_vec = all[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
        B_mat = matrix(B_vec, nrow = 1)

        #X_n = matrix(state[1:n.X], nrow = 3, ncol = N + 1)
        #B_vec = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta    
        #B_mat = matrix(B_vec, nrow = 3)

        # all the elements of theta should be positive
        #if (min(theta) <= 0)
        #return(-Inf)

        # Extracting observed data
        X_t = X_n[seq(2, N + 1, N / K)]


        # pi is the log of likelihood
        # This doesn't need a loop
        #p1 = 0
        ##print(dim(Y))
        #for (k in 1:K) {
            #Y.t = t(t(Y[, k]))
            #X_t.t = t(t(X_t[, k]))
            #p1 = p1 + t(Y.t - X_t.t) %*% inv_R %*% (Y.t - X_t.t)
        #}
        #p1 = -0.5 * p1
        #p1 = p1 - 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o)

        #######################################################################
        p1 = sum(dnorm(Y - X_t, sd = sqrt(R), log = TRUE))
        #- 0.5 * t(t(t(X_n[, 1])) - tau_o) %*% inv.lam_o %*% (t(t(X_n[, 1])) - tau_o))
        ######################################################################
        B_cov_gamma = gamma * (tau1 ^ 2) + (1 - gamma) * (tau0 ^ 2)
        p2 = dmvnorm(B_vec, sigma = diag(B_cov_gamma), log = TRUE)
        #p2 = (-1 / 2) * sum((B_vec - mu) ^ 2) / sigma2

        #f = mapply(drift_fun, X = split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), t = del_t * (0:N), MoreArgs = list(B_vec))
        f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(B_vec))
        #f = sapply(split(X_n, rep(1:ncol(X_n), each = nrow(X_n))), drift_fun, B_vec, list(1,2))
        del_X = diff(X_n)
        beta_tmp = sum((del_X / del_t - f[- (N + 1)]) ^ 2) * del_t / 2
        p3 = -(a4 + N / 2) * log(b4 + beta_tmp)

        return(p1 + p2 + p3)

    }

    sample_gamma <- function(B_vec) {
        gamma = numeric(length = n.theta)
        for (i in 1:n.theta) {
            prob = q[i] * dnorm(B_vec[i], sd = tau1) / (q[i] * dnorm(B_vec[i], sd = tau1) + (1 - q[i]) * dnorm(B_vec[i], sd = tau0))
            gamma[i] = rbinom(1, 1, prob)
        }
        return(gamma)
    }

    MH.X <- function(init, n, scale, gamma, B_vec) {
        chain = matrix(, nrow = n, ncol = n.X)
        accept.prob = 0
        for (i in 1:n) {
            prop = sapply(init, function(t) rnorm(1, t, scale))
            prop_ludf = c(prop, B_vec)
            init_ludf = c(init, B_vec)
            if (log(runif(1)) < (ludfun(prop_ludf, gamma) - ludfun(init_ludf, gamma))) {
                init = prop
                accept.prob = accept.prob + 1
            }
            chain[i,] = init
        }
        ans = list(chain, accept.prob / n)
        return(ans)
    }

    MH.B <- function(index, init, n, scale, gamma, state) {
        chain = numeric(length = n)
        accept.prob = 0
        prop_ludf = state
        init_ludf = state
        for (i in 1:n) {
            prop = rnorm(1, init, scale)
            prop_ludf[n.X + index] = prop
            init_ludf[n.X + index] = init
            if (log(runif(1)) < (ludfun(prop_ludf, gamma) - ludfun(init_ludf, gamma))) {
                init = prop
                accept.prob = accept.prob + 1
            }
            chain[i] = init
        }
        ans = list(chain, accept.prob / n)
        return(ans)
    }

    linchpin <- function(n, init, scale_vec) {
        X_avg = numeric(length = n.X)
        param_mat = matrix(, nrow = n, ncol = 2 * n.theta + n.sigma)
        scale = rep(0.0001 * 1, n.X + n.theta)

        scale[(n.X + 1):(n.X + n.theta)] = scale_vec # tf = 5 Sigma = 0.1  change of above
        scale[n.X + c(1)] = 5.7 * scale_vec[c(1)]
        scale[n.X + c(2)] = 0.2 * scale_vec[c(2)]
        scale[n.X + c(3)] = 0.06 * scale_vec[c(3)]
        scale[n.X + c(4)] = 0.1 * scale_vec[c(4)]

        scale.X = 0.0025

        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec # tf=5 Sigma = 1
        #scale[n.X + c(1)] = 5.8 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 0.03 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 0.08 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 0.17 * scale_vec[c(4)]

        #scale.X = 0.007 # 0.588612 0.436650 0.244172 0.189654 0.202868

        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec # tf=20
        #scale[n.X + c(1)] = 1.7 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 1.5 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 1 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 1 * scale_vec[c(4)]

        #scale.X = 0.0085


        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec   ### run this for 1
        #scale[n.X + c(1)] = 2.55 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 0.06 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 0.1 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 0.3 * scale_vec[c(4)]

        #scale.X = 0.01

        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec  # Sigma = 0.1  0.46056 0.47489 0.67292 0.30559 0.27525
        #scale[n.X + c(1)] = 4.4 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 0.033 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 0.08 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 0.2 * scale_vec[c(4)]

        #scale.X = 0.0038

        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec # Sigma = 0.1  change of above
        #scale[n.X + c(1)] = 5 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 0.2 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 0.2 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 0.2 * scale_vec[c(4)]

        #scale.X = 0.0042

        #scale[(n.X + 1):(n.X + n.theta)] = scale_vec  # Sigma = 10
        #scale[n.X + c(1)] = 0.25 * scale_vec[c(1)]
        #scale[n.X + c(2)] = 0.15 * scale_vec[c(2)]
        #scale[n.X + c(3)] = 0.4 * scale_vec[c(3)]
        #scale[n.X + c(4)] = 0.7 * scale_vec[c(4)]

        #scale.X = 0.027
        scale.B = scale[(n.X + 1):(n.X + n.theta)]

        accept.prob = numeric(1 + n.theta)
        state = init

        for (i in 1:n) {
            gamma = sample_gamma(init[(n.X + 1):(n.X + n.theta)])
            param_mat[i, (n.theta + n.sigma + 1):(2 * n.theta + n.sigma)] = gamma

            if (i %% (n / 5) == 0) {
                print(i)
                print(accept.prob / i)
                #print(c(i, accept.prob / i))
            }

            all = init
            chain = metrop(ludfun.X, initial = init[1:n.X], nbatch = 1, scale = scale.X, gamma = gamma, all = all)
            accept.prob[1] = accept.prob[1] + chain$accept
            state[1:n.X] = chain$batch

            #ans = MH.X(init[1:n.X], 1, scale.X, gamma, init[(n.X + 1):(n.X + n.theta)])
            #accept.prob[1] = accept.prob[1] + ans[[2]]
            #init[1:n.X] = ans[[1]]

            for (j in 1:n.theta) {
                #all = init
                #chain = metrop(ludfun, initial = init[n.X + j], nbatch = 1, scale = scale.B[j], gamma = gamma, all = all, index = j)
                #accept.prob[j + 1] = accept.prob[j + 1] + chain$accept
                #init[n.X + j] = chain$batch
               # state = init  ################### this was not there in previous code
                ans = MH.B(j, init[n.X + j], 1, scale.B[j], gamma, state)
                accept.prob[j + 1] = accept.prob[j + 1] + ans[[2]]
                state[n.X + j] = ans[[1]]
            }
            #state = init
            #chain = metrop(ludfun, init, 1, scale = scale, gamma = gamma)
            #state = chain$batch
            #accept.prob = accept.prob + chain$accept

            X_n = state[1:n.X]
            theta = state[(n.X + 1):(n.X + n.theta)] # vector of \sigma, \rho and \beta 
            X_avg = X_avg + state[1:n.X]
            param_mat[i, 1:n.theta] = theta

            f = mapply(drift_fun, X = X_n, t = del_t * (0:N), MoreArgs = list(theta))
            del_X = diff(X_n)
            beta_tmp = sum((del_X / del_t - f[ - (N + 1)]) ^ 2) * del_t / 2
            Sigma = rinvgamma(1, shape = N / 2 + a4, rate = b4 + beta_tmp)

            param_mat[i, (n.theta + 1):(n.theta + n.sigma)] = Sigma
            init = state
        }
        print(accept.prob / n)
        X_avg = X_avg / n
        final_output = list(param_mat, X_avg, accept.prob / n)
        return(final_output)
    }


    # Numerical method to sample from SDE
    euler_maruyama <- function(X0, del_t, N, theta, Sigma) {
        X = numeric(length = N+1)
        X[1] = X0
        for (i in 2:(N + 1))
            X[i] = X[i - 1] + (drift_fun_true(X[i - 1], theta)) * del_t + rnorm(1, sd = sqrt(del_t * Sigma))
        return(X)
    }
    # X = euler_maruyama(c(1,1,1), 0.1, 20, c(1,2,3), diag(2,3))


    # hyper-parameters
    to = 0 # initial time
    tf = 5 # final time
    Nobs = 50 # no of observations (Y) per time step
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
    b4 = 1
    tau1 = 10
    tau0 = 0.5

    K = (tf - to) * Nobs # no of real life observations, i.e. size of Y
    N = (tf - to) / del_t # no of discretizations of the Lorenz-63, i.e. size of X
    burn_in = 0 #5000 / del_t
    R = 0.04 #diag(1 / 2, 3) # observational error
    inv_R = 1/(0.04) #solve(R)
    mu = 0
    sigma2 = 10
    mu_truth = c(-2, 0, 0, 0) #c(-10, 28, 0, 10, -1, rep(0, 3), -8 / 3, rep(0, 11), 1, rep(0, 4), -1, rep(0, 7))
    non_zero = 1
    param_i = 1
    n.X =  (N + 1)
    n.theta = 4
    n.sigma = 1
    n.param = n.X + n.theta + n.sigma
    q = rep(0.5, n.theta) #runif(n.theta)
    n <- 1e3

    X = euler_maruyama(-0.5, del_t, N, 2, .1) # generating sample from Lorenz-63
    #X = X_total[, (burn_in):(N + burn_in)]
    #load('../../burninX')
    Y = X[seq(2, N + 1, N / K)] + rnorm(K, sd = sqrt(R)) # observations from Lorenz-63
    init = numeric(n.X + n.theta)
    init[(1:n.X)] <- as.numeric(X) #runif(n.param, 0, 5)
    init[(n.X+1):(n.X+n.theta)] = mu_truth

    #init[(n.X + 1):(n.X + n.theta)] <- rmvnorm(1, mu_truth, sigma = diag(1 / 50, n.theta))
    #non_zero = 1
    #param_i = 1
    #load("../../l63_linch_reg_bsv_0001_T_20_pv_10_init")
    #init[(n.X + 1):(n.X + n.theta)] <- head(tail(ans[[1]], 1)[1, - c(1, 2, 3)], -3)
    #init[n.X + 5] = -0.8

    #load('OU_linch_1e6_Nobs_init') # tf=20
    #load('OU_linch_5e5_Nobs_init_tf')
    #init[1:n.X] = ans[[2]]
    #init[(n.X + 1):(n.X + n.theta)] = rmvnorm(1, mean = mu_truth, sigma = diag(0.2, n.theta))

    sigma_Y = var(Y)
    tau0 = sqrt(sigma_Y / (10 * K))
    tau1 = sqrt(sigma_Y * max((n.theta ^ 2.1) / (100 * K), log(K)))*20

    #load('ou_linch_1e4_1')
    #var1 = cov(to_save[[1]][[1]][, 1:n.theta])
    #scale_vec = 4 * sqrt(diag(var1))
    #ans2 = linchpin(n, init, scale_vec)
    load('ou_linch_tf_1e4_1')
    var1 = cov(ans2[[1]][, 1:n.theta])
    scale_vec = 1 * sqrt(diag(var1))
    ans = linchpin(n, init, scale_vec)
    #plot.ts(ans[[1]][, param_i])
    #plot.ts(ans[[1]][, non_zero])
    chain_info = capture.output(cat("no of samples from MC is ", n, " \n starting from truth ", "\n priors spike slab ", " time period ",
                                tf, " Sigma is 1"))

    print(chain_info)
    to_save = list(ans, chain_info)
    #save(to_save, file = "ou_linch_T_2_1e5_cwise_1_spikes_truth_Sigma_10")
    pm = ans[[1]][, 1:(n.sigma + n.theta)]

    print(colMeans(pm))

    pm2 = ans[[1]][, (n.sigma + n.theta + 1):(n.sigma + 2 * n.theta)]
    print(colMeans(pm2))