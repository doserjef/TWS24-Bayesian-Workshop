# 5b-lmm-nimble.R: this script displays how to fit a Bayesian LMM using 
#                  NIMBLE. For our example, we will use a hypothetical data 
#                  set on little owls, where we are seeking to understand how
#                  wing length varies across populations.
rm(list = ls())
# Load the nimble package to fit the Bayesian model with MCMC
library(nimble)
# For convergence diagnostics
library(MCMCvis)
library(coda)
# Needed for calculating LOO
library(loo)
# For some plots
library(ggplot2)

# Load the data and make some EDA plots -----------------------------------
load('data/little-owl-data.rda')
str(little.owl.dat)

boxplot(wing ~ pop, little.owl.dat)
# Does the variation in wing length appear to be larger across populations
# or within a population? How do you know?

# Fit a fixed effects model (a linear model) ------------------------------
n <- nrow(little.owl.dat)
nimble.consts <- list(n = n)
# Notice how the pop variable is passed into NIMBLE. We generate "dummy coded"
# variables for each population to fit the model with a means parameterization.
nimble.data <- list(wing = little.owl.dat$wing,
                    pop1 = ifelse(little.owl.dat$pop == 1, 1, 0),
                    pop2 = ifelse(little.owl.dat$pop == 2, 1, 0),
                    pop3 = ifelse(little.owl.dat$pop == 3, 1, 0),
                    pop4 = ifelse(little.owl.dat$pop == 4, 1, 0),
                    pop5 = ifelse(little.owl.dat$pop == 5, 1, 0))
# Initial values
nimble.inits <- list(beta.1 = 0, beta.2 = 0, beta.3 = 0, beta.4 = 0,
                     beta.5 = 0, sigma = 1)
# Write the model
nimble.fixed <- nimbleCode({
  # Priors ----------------------------
  # Pop 1
  beta.1 ~ dnorm(0, var = 1000)
  # Pop 2
  beta.2 ~ dnorm(0, var = 1000)
  # Pop 3
  beta.3 ~ dnorm(0, var = 1000)
  # Pop 4
  beta.4 ~ dnorm(0, var = 1000)
  # Pop 5
  beta.5 ~ dnorm(0, var = 1000)
  # Residual standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.1 * pop1[i] + beta.2 * pop2[i] + beta.3 * pop3[i] +
             beta.4 * pop4[i] + beta.5 * pop5[i]
    wing[i] ~ dnorm(mu[i], sd = sigma)
  }
})
nimble.model <- nimbleModel(code = nimble.fixed,
                            constants = nimble.consts, data = nimble.data,
                            inits = nimble.inits)
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.1', 'beta.2',
                                                        'beta.3', 'beta.4',
                                                        'beta.5', 'sigma'),
                             enableWAIC = TRUE)
nimble.mcmc <- buildMCMC(nimble.conf)
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.fixed <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                     thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

# The model converged, right?

# Fit a random effects model ----------------------------------------------
n <- nrow(little.owl.dat)
# Number of populations
n.pop <- length(unique(little.owl.dat$pop))
# Populations as a numeric vector
pop.nimble <- as.numeric(little.owl.dat$pop)
nimble.consts <- list(n = n, n.pop = n.pop, pop = pop.nimble)
nimble.data <- list(wing = little.owl.dat$wing)
# Initial values
nimble.inits <- list(mu.pop = 0, sigma = 1, sigma.pop = 1)
# Write the model
nimble.random <- nimbleCode({
  # Priors ----------------------------
  # Overall mean across populations
  mu.pop ~ dnorm(0, var = 10000)
  # SD across populations
  sigma.pop ~ dinvgamma(.001, .001)
  # Residual standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Population level random effects
  for (i in 1:n.pop) {
    beta.pop[i] ~ dnorm(mu.pop, sd = sigma.pop)
  }
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.pop[pop[i]]
    wing[i] ~ dnorm(mu[i], sd = sigma)
  }
})
nimble.model <- nimbleModel(code = nimble.random,
constants = nimble.consts, data = nimble.data,
inits = nimble.inits)
nimble.conf <- configureMCMC(nimble.model, monitors = c('mu.pop', 'sigma.pop',
                                                        'sigma', 'beta.pop'),
                             enableWAIC = TRUE)
nimble.mcmc <- buildMCMC(nimble.conf)
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.random <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                      thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

MCMCsummary(out.random$samples)

# Is the standard deviation of the population random effect larger than
# the residual standard deviation? What does this mean?

# Extract the random effect estimates for each population. Note that these are
# all differences from the overall population mean.
beta.pop.samples <- MCMCchains(out.random$samples, exact = FALSE, param = 'beta.pop')
# Get the mean values (of course could just look in the above summary as well.
apply(beta.pop.samples, 2, mean)

# Next extract the estimated means for each population. Note that we need to add
# the extracted random effects above to the overall mean
# Extract the overall mean wing length samples
mu.pop.samples <- MCMCchains(out.random$samples, param = 'mu.pop')
pop.wing.length.samples <- beta.pop.samples + mu.pop.samples[, 1]

# Generate a boxplot showing the estimates --------------------------------
plot.df <- data.frame(median = apply(pop.wing.length.samples, 2, median),
                      lowest = apply(pop.wing.length.samples, 2, quantile, 0.025),
                      low = apply(pop.wing.length.samples, 2, quantile, 0.25),
                      high = apply(pop.wing.length.samples, 2, quantile, 0.75),
                      highest = apply(pop.wing.length.samples, 2, quantile, 0.975),
                      population = levels(little.owl.dat$pop))

ggplot(data = plot.df, aes(x = population)) +
  geom_boxplot(aes(ymin = lowest, lower = low, middle = median,
                   upper = high, ymax = highest), stat = 'identity',
               fill = 'lightskyblue1') +
  theme_bw(base_size = 16) +
  labs(x = 'Population', y = 'Little owl wing length')