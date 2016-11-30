library(data.table)

## prior: Beta
k <- 10000
alpha <- 2
beta <- 6
pi <- 0.5
n <- 100
sampled <- data.table(prop = dbeta(runif(k), alpha, beta))
sampled$dens <- dbinom(round(sampled$prop*n), size = n, prob = pi)
maxDens <- max(sampled$dens)
sampled$accepted <- ifelse(rbeta(k,alpha,beta) < sampled$dens/maxDens, TRUE, FALSE)
plot(density(sampled$prop))
hist(sampled$prop[sampled$accepted], freq = FALSE, col = "grey", breaks = 100)
temp <- rbeta(k, alpha + round(sampled$prop*n), beta + n - round(sampled$prop*n))
lines(density(temp), col = "red")
