#
# ST 502 Project 3
# December 1, 2017
#
# This program is a simulation study designed to investigate the effects
# of the difference in means, population variances, and sample sizes on the Type I
# error and power of the two-sample t-test when the variances are assumed equal vs.
# when they are assumed unequal.
#

library(ggplot2)
library(dplyr)

set.seed(888)

# Function to estimate the Type I error using randomly generated data
# Doesn't depend on the difference in means; the null is always that the means are equal
# par is a vector of n1, n2, var1, var2, and var.equal
TypeISim <- function(N, par) {
  # Name each of the parameters
  n1 <- par[1]
  n2 <- par[2]
  var1 <- par[3]
  var2 <- par[4]
  var.equal <- par[5]
  
  # Generate data from normal distributions with the given variances
  sample1 <- replicate(N, rnorm(n1, mean=0, sd=sqrt(var1))) # has n1 rows and N columns
  sample2 <- replicate(N, rnorm(n2, mean=0, sd=sqrt(var2))) # has n2 rows and N columns
  # (It doesn't matter what the means are, just that the difference in means is 0)
  
  # Get the p-values for the two-sample t-tests for the difference in means
  p <- rep(NA, N)
  for (i in 1:N) {
    p[i] <- t.test(sample1[,i], sample2[,i], var.equal=var.equal)$p.value
  }
  
  # Return the proportion of tests that incorrectly rejected the null hypothesis
  return(sum(p < 0.05)/N)
}

# Function to estimate the Power using randomly generated data
# par is a vector of diff, n1, n2, var1, var2, and var.equal
PowerSim <- function(N, par, var.equal) {
  # Name each of the parameter values
  diff <- par[1]
  n1 <- par[2]
  n2 <- par[3]
  var1 <- par[4]
  var2 <- par[5]
  var.equal <- par[6]
  
  # Generate data from normal distributions with the given variances
  sample1 <- replicate(N, rnorm(n1, mean=0, sd=sqrt(var1))) # has n1 rows and N columns
  sample2 <- replicate(N, rnorm(n2, mean=diff, sd=sqrt(var2))) # has n2 rows and N columns
  # (It doesn't matter what the means are, just the difference)
  
  # Get the p-values for the two-sample t-tests for the difference in means
  p <- rep(NA, N)
  for (i in 1:N) {
    p[i] <- t.test(sample1[,i], sample2[,i], var.equal=var.equal)$p.value
  }
  
  # Return the proportion of tests that correctly rejected the null hypothesis
  return(sum(p < 0.05)/N)
}

# Create data frame with an observation for each combination of the following parameter values
# Differences in means: diff = -5, -1, 0, 1, 5 (only needed for power calculations)
# Sample 1 sizes: n1 = 10, 25, 60
# Sample 2 sizes: n2 = 10, 25, 60
# Group 1 variances: var1 = 1, 3, 9
# Group 2 variance: var2 = 1
# Variances equal: var.equal = TRUE, FALSE

# -- Type I Error parameters
par1 <- expand.grid(c(10, 25, 60),
                    c(10, 25, 60),
                    c(1, 3, 9),
                    c(1),
                    c(TRUE, FALSE))
colnames(par1) <- c("n1", 
                    "n2", 
                    "var1", 
                    "var2",
                    "var.equal")

# -- Power parameters
par2 <- expand.grid(c(-5, -1, 0, 1, 5),
                    c(10, 25, 60),
                    c(10, 25, 60),
                    c(1, 3, 9),
                    c(1),
                    c(TRUE, FALSE))
colnames(par2) <- c("diff", 
                    "n1", 
                    "n2", 
                    "var1", 
                    "var2",
                    "var.equal")


# Estimate the Type I error and power for each combination of the paramter values
N <- 10000
Type.I <- apply(par1, FUN=TypeISim, MARGIN=1, N=N)
Power <- apply(par2, FUN=PowerSim, MARGIN=1, N=N)

# Put the results in a data frame with the parameter values
results1 <- data.frame(par1, Type.I)
results2 <- data.frame(par2, Power)

# Create factor versions of variables and rename levels so we can plot them
results1 <- results1 %>% 
  mutate(n1 = as.factor(n1),
         n2 = as.factor(n2),
         var1 = as.factor(var1),
         var.equal = ifelse(var.equal == TRUE, "Equal variance", "Unequal variance"))
results2 <- results2 %>% 
  mutate(n1 = as.factor(n1),
         n2 = as.factor(n2),
         var1 = as.factor(var1),
         var.equal = ifelse(var.equal == TRUE, "Equal variance", "Unequal variance"))

# Plot the results

# -- Type I error
pdf("TypeI.pdf", width=7, height=5.25)
ggplot(results1, aes(x=var1, y=Type.I, fill=interaction(var.equal, var1))) + 
  geom_col(position="dodge") +
  scale_fill_manual(values=c("#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4")) +
  facet_grid(n1 ~ n2) +
  geom_hline(yintercept=0.05, color="gray50", lty=2) +
  labs(x="Sample 1 Variance", y="Type I Error") +
  guides(fill=guide_legend(title=NULL))
dev.off()

# -- Power
pdf("Power.pdf", width=7, height=5.25)
ggplot(results2, aes(x=diff, y=Power)) + 
  geom_line(aes(color=var1, lty=var.equal)) +
  facet_grid(n1 ~ n2) +
  labs(x=expression(mu[1] - mu[2]), y="Power") +
  guides(color=guide_legend(title=expression(sigma[1]^2)),
         lty=guide_legend(title=NULL))
dev.off()

