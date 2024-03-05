#' Robin Elahi
#' Maximum likelihood Lab 1
#' 4 March 2024

library(tidyverse)
theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank()))

#### FUNCTIONS ####

## negative log-likelihood of y-values given intercept, slope, and error sd
## we pass the parameters as a vector in the order a, b, sigma
lm_nll <- function(par, x, y) {
  a <- par["a"]; b <- par["b"]; sigma <- exp(par["log_sigma"])
  -sum(dnorm(y, mean = a + b*x, sd = sigma, log = TRUE))
}

## MLE for linear regression
lm_fit <- function(x, y, SE = "fisher", init = c("a" = 0, "b" = 0, "log_sigma" = 0), n_boot = 1000) {
  
  if (!(SE %in% c("bootstrap", "fisher"))) stop("SE options are bootstrap or fisher")
  if (length(x) != length(y)) stop("x and y lengths differ")
  
  ## Maximum-likelihood estimate
  MLE <- optim(
    par = init, # initial values
    fn = lm_nll,  # our negative log-likelihood function
    x = x,   # x values
    y = y,   # y values
    hessian = (SE == "fisher") # only return Hessian if Fisher information used
  )
  
  ## Standard error using either fisher information or bootstrapping
  if (SE == "fisher") {
    var_cov <- solve(MLE$hessian)
    fit_se <- sqrt(diag(var_cov))
  } else {
    n_obs <- length(x) ## number of observations
    
    ## initialize an empty matrix to store parameter estimates for each sample
    par_samples <- matrix(nrow = n_boot, ncol = 3)
    
    ## recalculate MLE for each sample
    for (i in 1:n_boot) {
      samp <- sample(1:n_obs, size = n_obs, replace = TRUE) # sampled observations
      new_x <- x[samp] # subset x to sampled observations
      new_y <- y[samp] # subset y to sampled observations
      
      sample_fit <- optim(
        par = init, # initial values
        fn = lm_nll, # log likelihood function
        x = new_x, # sampled x-values
        y = new_y # sampled y-values
      )
      
      par_samples[i,] <- sample_fit$par # store sample parameters
    }
    fit_se <- apply(par_samples, MARGIN = 2, FUN = sd) # calculate column SDs
  }
  
  ## nicely format output
  data.frame(
    coef = c("intercept", "slope"), 
    estimate = MLE$par[1:2], 
    SE = fit_se[1:2]
  )
  
}

#### LAB EXERCISE 1 ####

a <- 1 ## intercept
b <- 2 ## slope
sigma <- 1 ## error standard deviation
n <- 50 ## sample size

## simulate data
sim_data <- tibble(
  x = runif(n, min = -3, max = 3), ## x values
  y = rnorm(n, mean = a + b*x, sd = sigma) ## regression equation
)

## plot data
sim_data %>% ggplot(aes(x, y)) + geom_point() + geom_smooth(method = "lm")

## run regression
lm_fit(sim_data$x, sim_data$y, SE = "fisher")

#### LAB EXERCISE 2 ####

