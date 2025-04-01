
### parameter selection for kernel ridge regression ###

## package 'DRR' is needed for running kernel ridge regression

## package 'CVST' is needed for cross-validation

# install.packages("DRR")
library(CVST)
library(DRR)
library(dplyr)

set.seed(999)

## simulate data
n <- 100000   # sample size
ep_x <- rnorm(n, 0, 1)   # first stage error
ep_y <- rnorm(n, 0, 0.1) # second stage error
u <- rnorm(n, 0, 1)     # confounder
z <- rbinom(n, 1, 0.5)  # binary IV
mu <- runif(n, 0, 1)    # random variable for nonconstant but rank-preserving genetic effect
xi <- rnorm(n, 0, 1)    # random variable for rank-varying genetic effect

## define first stage model

## impute parameters according to Table 3.1
# alpha0 <- 15, 14.5, 14
# alphaZ <- 1, 0.5, 0.25, 0
# alphaL <- 1, 0.85, 0.5, 0.45, 0.38, 0.25, 0
# alphaM <- 1, 0
# alphaU <- 2, 2.5

x <- alpha0 + alphaZ*z + alphaL*(z*exp(mu) + (1 - z)*mu) + alphaM*z*xi + alphaU*u + ep_x # first stage model

## define structural function
# h_x <- 0*x null;
# h_x <- 0.05*x^3 cubic;
# h_x <- 0.8*exp(x/4) exponential; 
# h_x <- ifelse(x>15, 0.5*(x - 15), 0) threshold

y <- h_x + u + ep_y # second stage model

## select data based on values of the instrument and rank according to values of the exposure
data <- data.frame(x, y, z)

data_0 <- data %>% filter(z == 0) %>% arrange(x)
data_1 <- data %>% filter(z == 1) %>% arrange(x)

## generate t
n_1 <- sum(z)
n_0 <- n - n_1

t_0 <- seq(1/n_0, 1, by = 1/n_0)
t_1 <- seq(1/n_1, 1, by = 1/n_1)

## construct data for kernel ridge regression
dat_y1 <- constructData(t_1, data_y1$y)
dat_y0 <- constructData(t_0, data_y0$y)
dat_x1 <- constructData(t_1, data_x1$x)
dat_x0 <- constructData(t_0, data_x0$x)

## construct kernel ridge regression
krr <- constructKRRLearner()

## construct parameter space for selection
sigmas=10^((1:9)/3)
lambdas= 10^(-8:0)
params <- constructParams(kernel="rbfdot", sigma=sigmas, lambda=lambdas) # choose the Gaussian kernel

## select parameters using 10-folds cross-validation
opt_y1 <- CV(dat_y1, krr, params, fold=10, verbose=FALSE)
opt_y0 <- CV(dat_y0, krr, params, fold=10, verbose=FALSE)
opt_x1 <- CV(dat_x1, krr, params, fold=10, verbose=FALSE)
opt_x0 <- CV(dat_x0, krr, params, fold=10, verbose=FALSE)

## when parameters differ, choose the minimum (i.e. the least peanalised parameter)
sigma = min(opt_y1[[1]]$sigma, opt_y0[[1]]$sigma, opt_x1[[1]]$sigma, opt_x0[[1]]$sigma)
lambda = min(opt_y1[[1]]$lambda, opt_y0[[1]]$lambda, opt_x1[[1]]$lambda, opt_x0[[1]]$lambda)

