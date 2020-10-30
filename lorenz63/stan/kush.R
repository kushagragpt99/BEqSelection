library(mcmcse)
library(pracma)
library(rstan)


changepoint_uni <- function(chain) {
    n = length(chain)
    z = (cumsum(chain) - (1:n) * mean(chain)) / sqrt(n)
    cov = as.numeric(sqrt(mcse.multi(chain, method = "bartlett")$cov))

    brownian = z / cov
    for (i in 1:n) {
        if (abs(brownian[i]) > abs(brownian[i + 1])) {
            breakmax = i + 1
            break
        }
    }

    brownianwt = numeric(n - 1)
    for (i in 1:n - 1)
        brownianwt[i] = brownian[i] / sqrt((i / n) * (1 - i / n))

    breakpt = which.max(brownian)
    breakptwt = which.max(brownianwt)

    print(paste('\nSquared CUSUM method:\nThe estimated break pt. is at the index:', breakpt,
              '\nSquared CUSUM-wt method:\nThe estimated break pt. is at the index:', breakptwt, sep = ""))

    plot.ts(chain, main = "Sampled Markov chain", ylab = "chain", xlab = "Time")
    abline(v = c(breakmax, breakpt, breakptwt), col = c("red", "blue", "green"))

    #plot.ts(brownian, main = "CUSUM scaled process", ylab = "", xlab = "Time")
    #abline(v = breakpt, col = "red")

    #plot.ts(brownianwt, main = "CUSUM scaled process", ylab = "", xlab = "Time")
    #abline(v = breakptwt, col = "blue")
}

changepoint_multi <- function(chain) {
    p = ncol(chain)
    n = nrow(chain)
    z = matrix(0, nrow = n, ncol = p)
    #mean = matrix(colMeans(chain),1,p)

    z = apply(chain, 2, function(x, n = nrow(chain))(cumsum(x) - (1:n) * mean(x)) / sqrt(n)) #maybe later do z(0) = 0.
    cov = mcse.multi(chain, method = "bartlett")$cov
    covinv = solve(cov)
    brownian = z %*% (sqrtm(cov)$Binv)

    breakpts = numeric(p)
    breakpts = apply(abs(brownian), 2, function(x) { for (i in 1:length(x))
            if (x[i] > x[i + 1]) return(i + 1) })
            estimate = max(breakpts)
    which_chain = which.max(breakpts)
    print(paste('LocalMax method:\nThe estimated break pt. is at the index:', estimate,
              '\nBreak pt. corresponds to the chain in the component:', which_chain, sep = ""))

    ### old method continued
    qn = numeric(N)
    qnwt = numeric(N - 1)

    for (i in 1:N) {
        qn[i] = t(as.matrix(z[i,])) %*% covinv %*% as.matrix(z[i,])
        if (i != N)
            qnwt[i] = qn[i] / ((i / N) * (1 - i / N))
        }
    breakcusum = which.max(abs(qn))
    breakwt = which.max(abs(qnwt))
    print(paste('\nSquared CUSUM method:\nThe estimated break pt. is at the index:', breakcusum,
              '\nSquared CUSUM-wt method:\nThe estimated break pt. is at the index:', breakwt, sep = ""))


    plot.ts(chain, main = "Sampled Markov chain", ylab = "chain", xlab = "Time")
    abline(v = c(estimate, breakcusum, breakwt), col = c("red", "blue", "green"))
    #legend("topright",legend=c("LocalMax", "Squared", "Squared-wt"),col=c("red", "blue","green"), lty=1:2, cex=0.8)

    plot.ts(brownian, main = "Non-Squared CUSUM", ylab = "chain", xlab = "Time")
    abline(v = c(estimate, breakcusum, breakwt), col = c("red", "blue", "green"))
    #legend("topright",legend=c("LocalMax", "Squared", "Squared-wt"),col=c("red", "blue","green"), lty=1:2, cex=0.8)

    plot(qn, main = "Squared CUSUM", ylab = "", xlab = "Time", type = "l")
    abline(v = breakcusum, col = "blue")
    plot(qnwt, main = "Squared CUSUM weighted", ylab = "", xlab = "Time", type = "l")
    abline(v = breakwt, col = "green")
}


changepoint <- function(chain) {
    if (is.vector(chain))
        return(changepoint_uni(chain))
    else
        return(changepoint_multi(chain))
    }

changepoint_modf <- function(chain) {
    for (i in 1:dim(chain)[2]) {
        changepoint_uni(chain[, i])
    }
}

load('output_stan/L63_HMC_chain_1_mtd5_bp_pv_10')
p1 = extract(fit, inc_warmup = TRUE, permuted = FALSE)
p1 = p1[, 1, 1:36]
changepoint_modf(p1)