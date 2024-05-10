# 4b-poisson-glm-nimble.R: this script displays how to fit a Bayesian
#                          generalized linear model using NIMBLE. For our 
#                          example, we will look at how hooded warbler
#                          relative abundance varies across Pennsylvania
# Data citation:
#      These are real data that come from the USGS Breeding Bird Survey, and
#      we load them through the spAbundance package. After loading spAbundance,
#      take a look at ?bbsData for the citation information.
rm(list = ls())
library(nimble)
library(ggplot2)
# Load spAbundance to get the hooded warbler data
library(spAbundance)
# For some maps
library(sf)
library(stars)
# For convergence diagnostics and calculating LOO
library(MCMCvis)
library(coda)
library(loo)

# Load data and do some reformatting --------------------------------------
data(bbsData)
# This is a list, formatted for fitting Bayesian models in spAbundance.
str(bbsData)
# Extract HOooded WArbler (HOWA) and combine with covariates into a data frame
bbs.dat <- data.frame(howa = bbsData$y['HOWA', ],
                      bbsData$covs)
str(bbs.dat)
# Look at ?bbsData for details on what these variables are. Briefly, they include
# a few different bioclimatic variables, a couple of land cover variables, and
# a few variables we think may influence the probability of observing an
# individual bird.
# We will compare three models that differ in the types of covariates they include
# to predict hooded warbler abundance:
# (1) Climate variables
# (2) Land cover variables
# (3) Climate and land cover variables
# Note that we include a linear effect of day of year as well as
# a linear effect of time of day to try and account for variation in when HOWA
# sings and is available to be detected.

# Fit climate only model --------------------------------------------------
n <- nrow(bbs.dat)
# Specify constants
nimble.consts <- list(n = n)
# Specify the data
nimble.data <- list(howa = bbs.dat$howa, bio2 = c(scale(bbs.dat$bio2)), 
                    bio8 = c(scale(bbs.dat$bio8)), bio18 = c(scale(bbs.dat$bio18)), 
                    forest = c(scale(bbs.dat$forest)), devel = c(scale(bbs.dat$devel)),
                    day = c(scale(bbs.dat$day)), tod = c(scale(bbs.dat$tod)), 
                    obs = c(scale(bbs.dat$obs)))
# Specify initial values
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, beta.3 = 0, 
                     beta.4 = 0, beta.5 = 0)
# Write the NIMBLE model
nimble.pois.glm <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # bio2
  beta.1 ~ dnorm(0, var = 1000)
  # bio8
  beta.2 ~ dnorm(0, var = 1000)
  # bio18
  beta.3 ~ dnorm(0, var = 1000)
  # day 
  beta.4 ~ dnorm(0, var = 1000)
  # tod 
  beta.5 ~ dnorm(0, var = 1000)
  # Likelihood ------------------------
  for (i in 1:n) {
    log(lambda[i]) <- beta.0 + beta.1 * bio2[i] + beta.2 * bio8[i] + beta.3 * bio18[i] + 
                      beta.4 * day[i] + beta.5 * tod[i]
    howa[i] ~ dpois(lambda[i])
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dpois(lambda[i])
  }
})
# Specify the model
nimble.model <- nimbleModel(code = nimble.pois.glm,
                            constants = nimble.consts, data = nimble.data,
                            inits = nimble.inits)
# The monitors argument is what we use to specify the parameters that we want
# to keep track of (i.e., save their MCMC samples).
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1',
                                                        'beta.2', 'beta.3',
                                                        'beta.4', 'beta.5',
                                                        'lambda', 'fitted',
                                                        'logProb_howa'),
                             enableWAIC = TRUE)
# Build the model
nimble.mcmc <- buildMCMC(nimble.conf)
# Compile the model
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
# Run the model
n.iter <- 30000
n.burn <- 20000
n.thin <- 10
n.chain <- 3
out.climate <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                       thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
# Has the model converged
MCMCsummary(out.climate$samples)
# Just look at the beta parameters
MCMCsummary(out.climate$samples, exact = FALSE, param = 'beta')

# Fit land cover only model -----------------------------------------------
# Specify initial values
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, beta.3 = 0, 
                     beta.4 = 0)
# Write the NIMBLE model
nimble.pois.glm <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # forest 
  beta.1 ~ dnorm(0, var = 1000)
  # devel 
  beta.2 ~ dnorm(0, var = 1000)
  # day 
  beta.3 ~ dnorm(0, var = 1000)
  # tod 
  beta.4 ~ dnorm(0, var = 1000)
  # Likelihood ------------------------
  for (i in 1:n) {
    log(lambda[i]) <- beta.0 + beta.1 * forest[i] + beta.2 * devel[i] + 
                      beta.3 * day[i] + beta.4 * tod[i]
    howa[i] ~ dpois(lambda[i])
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dpois(lambda[i])
  }
})
# Specify the model
nimble.model <- nimbleModel(code = nimble.pois.glm,
                            constants = nimble.consts, data = nimble.data,
                            inits = nimble.inits)
# The monitors argument is what we use to specify the parameters that we want
# to keep track of (i.e., save their MCMC samples).
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1',
                                                        'beta.2', 'beta.3',
                                                        'beta.4', 
                                                        'lambda', 'fitted',
                                                        'logProb_howa'),
                             enableWAIC = TRUE)
# Build the model
nimble.mcmc <- buildMCMC(nimble.conf)
# Compile the model
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
# Run the model
n.iter <- 30000
n.burn <- 20000
n.thin <- 10
n.chain <- 3
out.land <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                    thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE,
                    WAIC = TRUE)
# Has the model converged
MCMCsummary(out.land$samples)
# Just look at the beta parameters
MCMCsummary(out.land$samples, exact = FALSE, param = 'beta')

# Fit climate + land cover model ------------------------------------------
# Specify initial values
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, beta.3 = 0, 
                     beta.4 = 0, beta.5 = 0, beta.6 = 0, beta.7 = 0)
# Write the NIMBLE model
nimble.pois.glm <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # forest 
  beta.1 ~ dnorm(0, var = 1000)
  # devel 
  beta.2 ~ dnorm(0, var = 1000)
  # day 
  beta.3 ~ dnorm(0, var = 1000)
  # tod 
  beta.4 ~ dnorm(0, var = 1000)
  # bio2
  beta.5 ~ dnorm(0, var = 1000)
  # bio8
  beta.6 ~ dnorm(0, var = 1000)
  # bio18
  beta.7 ~ dnorm(0, var = 1000)
  # Likelihood ------------------------
  for (i in 1:n) {
    log(lambda[i]) <- beta.0 + beta.1 * forest[i] + beta.2 * devel[i] + 
                      beta.3 * day[i] + beta.4 * tod[i] + 
                      beta.5 * bio2[i] + beta.6 * bio8[i] + beta.7 * bio18[i]
    howa[i] ~ dpois(lambda[i])
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dpois(lambda[i])
  }
})
# Specify the model
nimble.model <- nimbleModel(code = nimble.pois.glm,
                            constants = nimble.consts, data = nimble.data,
                            inits = nimble.inits)
# The monitors argument is what we use to specify the parameters that we want
# to keep track of (i.e., save their MCMC samples).
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1',
                                                        'beta.2', 'beta.3',
                                                        'beta.4', 'beta.5', 
                                                        'beta.6', 'beta.7',
                                                        'lambda', 'fitted',
                                                        'logProb_howa'),
                             enableWAIC = TRUE)
# Build the model
nimble.mcmc <- buildMCMC(nimble.conf)
# Compile the model
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
# Run the model
n.iter <- 30000
n.burn <- 20000
n.thin <- 10
n.chain <- 3
out.full <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                    thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE,
                    WAIC = TRUE)
# Has the model converged
MCMCsummary(out.full$samples)
# Just look at the beta parameters
MCMCsummary(out.full$samples, exact = FALSE, param = 'beta')

# Compare the three candidate models --------------------------------------
# WAIC
out.climate$WAIC
out.land$WAIC
out.full$WAIC
# Calculate LOO using the loo package
log.prob.climate <- MCMCchains(out.climate$samples, params = 'logProb', exact = FALSE)
log.prob.full <- MCMCchains(out.full$samples, params = 'logProb', exact = FALSE)
log.prob.land <- MCMCchains(out.land$samples, params = 'logProb', exact = FALSE)
loo.climate <- loo::loo(log.prob.climate)
loo.full <- loo::loo(log.prob.full)
loo.land <- loo::loo(log.prob.land)