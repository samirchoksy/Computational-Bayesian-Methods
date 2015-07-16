# kalman filter simulation and parameter estimation
# author: samir choksy

# prior parameters
alpha = rnorm(n=1, mean=5, sd=.5)
beta = rnorm(n=1, mean=.4, sd=.01)
tau = rnorm(n=1, mean=.2, sd=.03)
sigma = rnorm(n=1, mean=.25, sd=.05)
mu = rnorm(n=1, mean=1, sd=.5)

# containers for rolling means and variances
Ex = rep(0,T-1)
Ey = rep(0,T-1)
Sdx = rep(0,T-1)
Sdy = rep(0,T-1)

# simulated data
T = 1000

# latent and observed states
x = rep(0,T)
x0 = rnorm(n=1, mean=alpha, sd=2)
x[1] = rnorm(n=1, mean=alpha+beta*x0, sd=tau)
y = rep(0,T)

# simulate
for(i in 2:T){
    x[i] <- rnorm(n=1, mean=alpha+beta*x[i-1], tau)
    y[i-1] <- rnorm(n=1, mean=mu+x[i-1], sd=sigma)

    # populate moving average mean/sd
    Ex[i-1] = mean(x[1:i])
    Sdx[i-1] = sd(x[1:i])
    Ey[i-1] = mean(y[1:i])
    Sdy[i-1] = sd(y[1:i])
}

# plot sample means and variances for latent and observed states
par(mfrow=c(2,2))
plot(Ex, type="l", col="blue", main="Latent State: E[x]")
abline(h=mean(x))
plot(Sdx, type="l", col="red", mean="sd[X]")
abline(h=sd(x))
plot(Ey, type="l", col="blue", main="Observed State: E[y]")
abline(h=mean(y))
plot(Sdy, type="l", col="red", main="sd[y]")
abline(h=sd(y)) 
