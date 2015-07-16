# parameter learning for a bayesian factor model with normal and 
# non-normal errors on simulated data with conjugate priors.
# (note: posterior distribution fully tractable with known
# marginal distributions. sampling is direct from posterior. gibbs
# and metropolis hastings may be used in more sophisiticated settings.
# see johannes (2014) for more detail.)

# author: samir choksy

rm(list=ls(all=TRUE))
library(MCMCpack)

# normal errors
# training data and posterior sample sizes
T = 100
N = 5000

# simulated data - true parameters
numCoeff = 3
beta = c(beta0 = 1, beta1 = -1, beta2 = 5)
sigma = 5
X = matrix(rnorm(n=T*numCoeff, mean=0, sd=1), ncol=numCoeff)
Y = X%*%beta + rnorm(n=T, mean=0, sd=sigma)

# prior hyperparameters for normal inverse gamma conjugate prior
a = c(0,0,0) #prior means
A = diag(numCoeff) #prior precision
b = 1 #inverse gamma prior shape
B = 1 #inverse gamma prior rate

# update hyperparameters
AT = solve((t(X)%*%X)+solve(A))
aT = AT%*%(t(X)%*%Y+solve(A)%*%a)
bT = b + T
BT = t(Y)%*%Y+t(a)%*%solve(A)%*%a-t(aT)%*%solve(AT)%*%aT+B

# sample from posterior
posteriorBeta = matrix(0, N, 3)
posteriorSig = rinvgamma(N, bT/2, BT/2)
for(i in 1:N){posteriorBeta[i,] = mvrnorm(1, aT, sqrt(posteriorSig[i]*AT))}

# prior hyperparameters for jeffrey's prior
AT = solve(t(X)%*%X)
beta.hat = AT%*%t(X)%*%Y
sig.sq = t(Y-X%*%beta.hat)%*%(Y-X%*%beta.hat)/(T-numCoeff)
posteriorBetaJeffreys = matrix(0, N, 3)
posteriorSigJeffreys = rinvgamma(N, (T-numCoeff)/2, ((T-numCoeff)*sig.sq)/2)
for(i in 1:N){posteriorBetaJeffreys[i,] = mvrnorm(1, beta.hat, sqrt(posteriorSig[i])*AT)}
