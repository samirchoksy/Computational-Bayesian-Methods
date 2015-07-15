# simple gibbs sampling for non-normal data
# robust likelihoods with laplace, student-t and mixed normal errors
# author: samir choksy

library(MCMCpack)
library(ghyp)
library(VGAM)

# the following gibbs sampling methods use normal inverse gamma parameter
# prior. the last example, with a mixture of normal errors, features a 
# beta prior for the mixing probabilities.

# location-scale model y_t = mu + sigma*epsilon_t, where epsilon is
# a double exponential/laplace distributed random error

# simulate data
set.seed(78957)
T = 1000 # training data
N = 5000  # posterior sample size

# true parameter values
mu = rnorm(1, mean=.2, sd=1)
sigma = rnorm(1, mean=2, sd=1)

y <- mu + sigma*rlaplace(T, location=0, scale=1) # also: rlaplace(T,mu,sigma) 

# prior hyperparameters
a = 0
A = 100
b = 10
B = 40

# posterior arrays
p_mu = rep(0,T)
p_sigma = rep(0,T)
lambda  = rexp(T, rate=2)

# update hyperparameters
mu_lambda = sum(y^2/lambda)/sum(1/lambda)
s_lambda = sum(y^2/lambda) - ((sum(y/lambda))^2)/(sum(1/lambda))

AT = 1/((sum(1/lambda))+(1/A))
aT = ((mu_lambda*(sum(1/lambda))) + (a/A))*AT
bT = b+T
BT = B + s_lambda + (mu_lambda-a)^2/((sum(1/lambda)^-1)+A)

# sample from full conditionals for 1:N
for(i in 1:N){
    # draw mu and sigma from inverse gamma and normal distribution
    # calculate posteriors from priors and double exponential likelihood

    p_sigma[i] = rinvgamma(1, bT/2, BT/2)
    p_mu[i] = rnorm(1, aT, sqrt(AT*p_sigma[i]))

    # draw lambda
    for(t in 1:T){
        chi = ((y[t]-p_mu[t])^2/(p_sigma[i]))
        if(chi<.001){chi = chi+.001} # float error cf. 'ghyp'
        lambda[t] = rgig(1, 1/2, chi, 1)
    }

    # update hyperparameters
    mu_lambda = sum(y/lambda)/sum(1/lambda)
    s_lambda = sum(y^2/lambda) - ((sum(y/lambda))^2)/(sum(1/lambda))

    AT = 1/((sum(1/lambda))+(1/A))
    aT = ((mu_lambda*(sum(1/lambda))) + (a/A))*AT
    BT = B + s_lambda + (mu_lambda-aT)^2/((sum(1/lambda)^(-1))+A)
}
