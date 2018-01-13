# 
# ST502 Project 2
# 
# This program uses a Markov Chain Monte Carlo sampling method to estimate the 
# parameters of a simple logistic regression model* and the sampling
# distributions of the parameters. At the end, we compare the results to the
# maximum likelihood fit. 
#
# * The logistic regression model aims to predict whether a field goal will be
# good or bad given the distance from which it was kicked. The data is from the
# 2008 NFL season.
# 
# Last edited: January 13, 2018 - Edited project description in header
# 





library(ggplot2)





# Expit function
  # 1 / (1 + exp(-beta0 - beta1*x))
  # parameter beta is a vector in the form of (beta0, beta1)
  # x is the vector of values of the predictor variable
expit <- function(beta, x) {
  # put x in matrix form (dimensions 2 x n)
  x <- cbind(rep(1, length(x)), x)
  # compute the expit function for the given beta and x
  return(1 / (1 + exp(-x %*% beta)))
}





# Prior density function
  # beta1, beta2 ~ Normal
  # parameter beta is a vector in the form of (beta0, beta1)
  # the joint density is the product of the marginal densities ?????????????????????????????
prior <- function(beta) {
  return(prod(dnorm(beta, mean=0, sd=10)))
}





# Likelihood function
  # y | beta1, beta2 ~ Bernoulli
  # parameter beta is a vector in the form of (beta0, beta1)
  # y is a vector of the observed values (0's for bad kicks and 1's for good kicks)
  # the joint density is the product of the marginal densities since observations are indep.
likelihood <- function(beta, x, y) {
  return(prod(dbinom(y, prob=expit(beta, x), size=1)))
}





# Function to generate values for the posterior distribution using MCMC sampling
  # x is a vector of predictor values
  # y is a vector of observed values
  # niter is the number of iterations
  # startvals is a vector in the form (beta0, beta1)
  # proposalsd is the sd of the jumping distribution
betasampler <- function(x, y, niter, startvals, proposalsd) {
  # Matrix to store sampled betas
  beta <- matrix(rep(0, niter), nrow=niter, ncol=2)
  # Use the starting value as the first beta
  beta[1,] <- startvals
  # Loop through the desired number of iterations
  for(i in 2:niter) {
    # Temp value for beta is previous beta
    currentbeta <- beta[i-1,]
    # Get the next random draw from the N(0, proposalsd^2) jumping distribution
    newbeta <- currentbeta + rnorm(2, 0, proposalsd)
    # Find r ratio
    r <- prior(beta=newbeta) * likelihood(beta=newbeta, x=x, y=y) /
      (prior(beta=currentbeta) * likelihood(beta=currentbeta, x=x, y=y))
    # Accept this value with prob. min(1, r)
    if(runif(1) < r) {
      beta[i,] <- newbeta
    }
    else {
      beta[i,] <- currentbeta
    } # end ifs
  } # end loop
  # Return the vector of theta values
  return(beta)
} # end function




# Function for plotting and formatting histograms using ggplot2
histogram <- function(x) {
  ggplot(data.frame(x), aes(x)) + 
    geom_histogram(bins=18, col="gray95") +
    labs(x=NULL, y="Frequency") +
    theme(plot.title=element_text(size=20),
      axis.title=element_text(size=16), 
      axis.text=element_text(size=12), 
      panel.background=element_rect(fill="gray95"), 
      strip.background=element_rect(fill="white"),
      plot.margin=margin(20,20,20,20)) + 
    scale_fill_manual(palette)
}










# Read in NFL field goal data from the csv file
fgdata <- read.csv("nfl2008_fga.csv", header=TRUE)





# Approximate the posterior distribution of beta0 and beta1 using the MCMC sampling function
randomDraws <- betasampler(x=fgdata$distance, y=fgdata$GOOD, 
                           niter=100000, startvals=c(0,0), proposalsd=0.05)
burnin <- 7500
posterior_b0 <- randomDraws[-c(1:burnin), 1]
posterior_b1 <- randomDraws[-c(1:burnin), 2]

# Histograms showing the distributions of beta0 and beta1
histogram(posterior_b0) + 
  labs(title=expression(paste("Posterior distribution of ", beta[0])))
histogram(posterior_b1) + 
  labs(title=expression(paste("Posterior distribution of ", beta[1])))
# Estimates
mean(posterior_b0)
mean(posterior_b1)




# Fit a logistic regression model with GOOD as the response and distance as the predictor
mod <- glm(GOOD ~ distance, data=fgdata, family=binomial(link="logit"))
summary(mod)





