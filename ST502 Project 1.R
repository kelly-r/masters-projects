#
# ST502 R Project 1
#
# Program comparing the performance of five types of intervals 
# for the parameter p in the binomial distribution using 
# simulated data for various values of p and n
#
# Last edited January 12, 2018 - Changed formatting of results
#

library(binom)

# Function to compute the properties of intervals
intstats <- function(x, truth) {
  N <- nrow(x)
  # Find instances where the truth was above the interval
  # and the proportion of the total
  above <- x[,2] < truth
  prop.above <- sum(above)/N
  # Find instances where the truth was below the interval
  # and the proportion of the total
  below <- x[,1] > truth
  prop.below <- sum(below)/N
  # Find the proportion of intervals containing the truth
  prop.capture <- 1 - prop.above - prop.below
  
  # Get the lengths
  length <- x[,3]
  # Calculate mean of the interval lengths and standard error
  mean <- mean(length)
  SE <- sqrt(var(length)/N)
  
  out <- c(prop.capture, prop.above, prop.below, mean, SE)
  names(out) <- c("Proportion captured", "Proportion above", 
                  "Proportion below", "Average length", 
                  "Standard error")
  return(out)
}

# Function to compute the five types of confidence intervals for 
# a given combination of n and p
confints <- function(n, p) {
  # Generate 5000 random samples from the Binomial distribution
  N <- 5000
  set.seed(1)
  simdata <- rbinom(N, n, p)
  # Find which values are equal to 0
  is0 <- (simdata==0)
  # Find which values are equal to n
  isn <- (simdata==n)
  # Find which values are not equal to 0 or n
  ok <- !(is0 + isn)
  
  # Calculate p.hat for each sample
  p.hats <- simdata/n
  
  # (a) Normality based interval
  
  z <- qnorm(1-0.05/2) # Quantile of N(0,1) distribution
  SEs <- sqrt(p.hats*(1-p.hats)/n) # Estimated standard errors
  int.a.lower <- p.hats - z*SEs # Lower bounds
  int.a.upper <- p.hats + z*SEs # Upper bounds
  # Special case where y = 0: Set interval to (0, 0)
  int.a.lower[is0] <- 0
  int.a.upper[is0] <- 0
  # Special case where y = n: Set interval to (1, 1)
  int.a.lower[isn] <- 1
  int.a.upper[isn] <- 1
  # Intervals
  int.a <- cbind(int.a.lower, 
                 int.a.upper, 
                 int.a.upper - int.a.lower)
  colnames(int.a) <- c("2.5%", "97.5%", "length")
  # Properties of the intervals
  a <- intstats(int.a, p)
  
  # (b) Exact interval
  
  int.b.lower <- rep(NA, N)
  int.b.upper <- rep(NA, N)
  # Special case where y = 0: Set interval to (0, 0)
  int.b.lower[is0] <- 0
  int.b.upper[is0] <- 0
  # Special case where y = n: Set interval to (1, 1)
  int.b.lower[isn] <- 1
  int.b.upper[isn] <- 1
  # Intervals
  ci <- binom.confint(simdata[ok], n=n, 
                      conf.level=0.95, methods="exact")
  int.b.lower[ok] <- ci$lower
  int.b.upper[ok] <- ci$upper
  int.b <- cbind(int.b.lower,
                 int.b.upper,
                 int.b.upper - int.b.lower)
  colnames(int.b) <- c("2.5%", "97.5%", "length")
  # Properties of the intervals
  b <- intstats(int.b, p)

  # (c) Score interval
  
  term1 <- (p.hats + z^2/(2*n))/(1 + z^2/n)
  term2 <- (z*sqrt((p.hats*(1-p.hats) + z^2/(4*n))/n))/(1 + z^2/n)
  int.c.lower <- term1 - term2
  int.c.upper <- term1 + term2
  int.c <- cbind(int.c.lower, int.c.upper, 2*term2) # Intervals
  colnames(int.c) <- c("2.5%", "97.5%", "length")
  # Properties of the intervals
  c <- intstats(int.c, p) 
  
  # (d) Parametric raw percentile interval
  
  # Generate B=500 samples from the Binom(n, p.hat) distribution
  B <- 500 
  set.seed(1)
  bootdata <- sapply(p.hats, 
                     FUN=function(x, B, n) rbinom(B, n, x), 
                     B=B, n=n)
  # Calculate p.hat for each bootstrap sample
  boot.p.hats <- bootdata/n
  # Construct the interval from the quantiles of the bootstrap 
  # distribution of p-hat
  int.d.lower <- apply(boot.p.hats, MARGIN=2, 
                       FUN=function(x) quantile(x, 0.025))
  int.d.upper <- apply(boot.p.hats, MARGIN=2, 
                       FUN=function(x) quantile(x, 0.975))
  # Special case where y = 0: Set interval to (0, 0)
  int.d.lower[is0] <- 0
  int.d.upper[is0] <- 0
  # Special case where y = n: Set interval to (1, 1)
  int.d.lower[isn] <- 1
  int.d.upper[isn] <- 1
  int.d <- cbind(int.d.lower, int.d.upper, 
                 int.d.upper - int.d.lower)
  colnames(int.d) <- c("2.5%", "97.5%", "length")
  rownames(int.d) <- NULL
  # Properties of the intervals
  d <- intstats(int.d, p)
  
  # (e) Parametric bootstrap interval
  
  # Calculate the standard error of p-hat for each bootstrap 
  # sample
  boot.SEs <- sqrt(boot.p.hats*(1-boot.p.hats)/n)
  # Make a matrix where the rows are just the p.hats used to 
  # generate the values repeated (for calculation of T)
  p.hats.matr <- matrix(rep(p.hats, each=500), nrow=500)
  # Calcualte the t-statistic for each sample using the p-hats 
  # and standard errors
  boot.Ts <- (boot.p.hats - p.hats.matr)/boot.SEs
  boot.Ts[is.infinite(boot.Ts)] <- NA
  # Find quantiles of the distribution
  t.bounds <- apply(boot.Ts, MARGIN=2, 
                    FUN=function(x) quantile(x, c(0.025, 0.975), 
                                             na.rm=TRUE))
  # Construct the interval from the quantiles of the bootsrap 
  # distr. of p-hat
  int.e.lower <- p.hats - t.bounds[2,]*SEs
  int.e.upper <- p.hats - t.bounds[1,]*SEs
  # Special case where y = 0: Set interval to (0, 0)
  int.e.lower[is0] <- 0
  int.e.upper[is0] <- 0
  # Special case where y = n: Set interval to (1, 1)
  int.e.lower[isn] <- 1
  int.e.upper[isn] <- 1
  int.e <- cbind(int.e.lower, int.e.upper, 
                 int.e.upper - int.e.lower)
  colnames(int.e) <- c("2.5%", "97.5%", "length")
  # Properties of the intervals
  e <- intstats(int.e, p)
  
  out <- data.frame(a, b, c, d, e)
  colnames(out) <- c("Normal", "Exact", "Score", "Perc.", "t")

  return(round(out,3))
}

# Make a grid containing all combinations of the parameters
par <- expand.grid(c(10, 100, 1000), c(0.05, 0.2, 0.5))
colnames(par) <- c("n", "p")
npar <- nrow(par)

# Store results in a list in order of the parameters
results <- list()
for (i in 1:npar) {
  results[[i]] <- confints(n=par[i,1], p=par[i,2])
}

# Set size of new matrices and then fill them
proportion.captured <- matrix(NA, npar, 5)
proportion.above <- matrix(NA, npar, 5)
proportion.below <- matrix(NA, npar, 5)
average.length <- matrix(NA, npar, 5)
standard.deviation <- matrix(NA, npar, 5)
for (i in 1:npar) {
  proportion.captured[i,] <- t(results[[i]][1,])
  proportion.above[i,] <- t(results[[i]][2,])
  proportion.below[i,] <- t(results[[i]][3,])
  average.length[i,] <- t(results[[i]][4,])
  standard.deviation[i,] <- t(results[[i]][5,])
}
# Set row and column names
rows <- rep(NA, npar)
for(i in 1:npar) {
  rows[i] <- paste0('n = ', par[i,1], ', p = ', par[i,2])
}
rownames(proportion.captured) <- 
  rownames(proportion.above) <- 
  rownames(proportion.below) <- 
  rownames(average.length) <- 
  rownames(standard.deviation) <- rows
colnames(proportion.captured) <- 
  colnames(proportion.above) <- 
  colnames(proportion.below) <- 
  colnames(average.length) <- 
  colnames(standard.deviation) <-
  c("Normal", "Exact", "Score", "Perc.", "t")

# View results
proportion.captured
proportion.above
proportion.below
average.length
standard.deviation
