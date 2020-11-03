library(tensorflow)
library(sgmcmc)
use_condaenv("r-tensorflow")
installTF()
dataset = list("x" = rnorm(1000))
params = list("theta" = 0)
logLik = function(params, dataset) {
    distn = tf$distributions$Normal(params$theta, 1)
    return(tf$reduce_sum(distn$log_prob(dataset$x)))
}
stepsize = list("theta" = 1e4)
sgld = sgldSetup(logLik, dataset, params, stepsize)
nIters = 10 ^ 3 # Initialize location estimate
locEstimate = 0 # Initialise TensorFlow session
sess = initSess(sgld)
for (i in 1:nIters) {
    sgmcmcStep(sgld, sess)
    locEstimate = locEstimate + 1 / nIters * getParams(sgld, sess)$theta
}