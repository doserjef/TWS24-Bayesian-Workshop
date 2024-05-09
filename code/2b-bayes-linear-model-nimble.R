# 2a-bayes-linear-model-nimble.R: this script displays how to fit a Bayesian linear
#                                 using NIMBLE. For our example, we will look at 
#                                 how bird species richness varies in relation to 
#                                 elevation and landowner type (which may result 
#                                 in different management practices)
# Author: Jeffrey W. Doser
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

# Load data and some EDA --------------------------------------------------
# Change working directory as needed
load("data/bird-elev-data.rda")
str(dat)
ggplot(data = dat, aes(x = elev, y = bird.rich)) +
  geom_point(size = 3, pch = 21, fill = 'lightskyblue') +
  theme_bw(base_size = 18) +
  labs(x = 'Elevation', y = 'Bird Richness')

# Fit a simple linear model -----------------------------------------------
n <- nrow(dat)
nimble.consts <- list(n = n)
# 2. Data -----------------------------
# Notice how the land owner variable is passed into NIMBLE. We need to use 
# a "dummy coding" approach and choose whether we want to use a means 
# parameterization or an effects parameterization.
nimble.data <- list(richness = dat$bird.rich, 
		    elev = dat$elev,
                    lo.state = ifelse(dat$land.owner == 'state', 1, 0), 
                    lo.private = ifelse(dat$land.owner == 'private', 1, 0))
# 3. Initial values -------------------
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, beta.3 = 0, sigma = 1)
# 4. Write the NIMBLE model -----------
nimble.linear <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # Elevation regression coefficient
  beta.1 ~ dnorm(0, var = 1000) 
  # Effect of state landowner (relative to federal)
  beta.2 ~ dnorm(0, var = 1000)
  # Effect of private landowner (relative to federal)
  beta.3 ~ dnorm(0, var = 1000)
  # Standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.0 + 
             beta.1 * elev[i] + 
             beta.2 * lo.state[i] + 
             beta.3 * lo.private[i]
    richness[i] ~ dnorm(mu[i], sd = sigma)
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dnorm(mu[i], sd = sigma)
  }
})
# 5. Specify the model ----------------
nimble.model <- nimbleModel(code = nimble.linear, 
			    constants = nimble.consts, data = nimble.data, 
			    inits = nimble.inits)
# 6. Configure the model --------------
# The monitors argument is what we use to specify the parameters that we want
# to keep track of (i.e., save their MCMC samples). 
# NOTE: the logProb_richness is needed to calculate loo later.
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1', 
							'beta.2', 'beta.3',
							'sigma', 'mu', 'fitted', 
							'logProb_richness'), 
                             enableWAIC = TRUE)
# 7. Build the model ------------------
nimble.mcmc <- buildMCMC(nimble.conf)
# 8. Compile the model ----------------
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
# 9. Run the model --------------------
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.linear <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
		      thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
                      WAIC = TRUE)
# Summarize the model
# NOTE: when you set WAIC = TRUE, there are two components to the resulting 
#       NIMBLE model fit: (1) samples and (2) WAIC. 
str(out.linear)
MCMCsummary(out.linear$samples)

# Has the model converged? Use both visual assessments as well as
# numeric criteria to support your decision
# Are the Rhat values close to 1? How large are the effective sample sizes?
# Plot some traceplots of the resulting posteriors
MCMCtrace(out.linear$samples)

# Check the residuals as usual
fitted.samples <- MCMCchains(out.linear$samples, params = 'fitted', exact = FALSE)
fitted.means <- apply(fitted.samples, 2, mean)
res.vals <- dat$bird.rich - fitted.means

# Is there a pattern in the residuals?
plot(fitted.means, res.vals,
     main = "Residuals vs. fitted values", las = 1,
     xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)

# What does this plot suggest to you? Can we continue analyzing this model?

# It could also be fun to plot the normal probability plot (qqplot).
qqnorm(res.vals, pch = 19)
qqline(res.vals, col = 'red')

# Fit a model with a linear and quadratic term for elevation --------------
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, beta.3 = 0, 
		     beta.4 = 0, sigma = 1)
nimble.full <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # Elevation linear regression coefficient
  beta.1 ~ dnorm(0, var = 1000) 
  # Elevation quadratic regression coefficient
  beta.2 ~ dnorm(0, var = 1000)
  # Effect of state landowner (relative to federal)
  beta.3 ~ dnorm(0, var = 1000)
  # Effect of private landowner (relative to federal)
  beta.4 ~ dnorm(0, var = 1000)
  # Standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.0 + 
             beta.1 * elev[i] + 
	     beta.2 * elev[i]^2 + 
             beta.3 * lo.state[i] + 
             beta.3 * lo.private[i]
    richness[i] ~ dnorm(mu[i], sd = sigma)
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dnorm(mu[i], sd = sigma)
  }
})
nimble.model <- nimbleModel(code = nimble.full,  
			    constants = nimble.consts, data = nimble.data, 
			    inits = nimble.inits)
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1', 
							'beta.2', 'beta.3',
							'beta.4',
							'sigma', 'mu', 'fitted', 
							'logProb_richness'), 
                             enableWAIC = TRUE)
nimble.mcmc <- buildMCMC(nimble.conf)
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.full <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                    thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
                    WAIC = TRUE)

# Has the model converged? Perform the usual tests.
MCMCsummary(out.full$samples)
# Check residuals
fitted.samples.full <- MCMCchains(out.full$samples, params = 'fitted', exact = FALSE)
fitted.means.full <- apply(fitted.samples.full, 2, mean)
res.vals.full <- dat$bird.rich - fitted.means.full

# Is there a pattern in the residuals?
plot(fitted.means.mean, res.vals.full,
     main = "Residuals vs. fitted values", las = 1,
     xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)

# Intercept only model ----------------------------------------------------
nimble.inits <- list(beta.0 = 0, sigma = 1)
nimble.int <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # Standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.0
    richness[i] ~ dnorm(mu[i], sd = sigma)
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dnorm(mu[i], sd = sigma)
  }
})
nimble.model <- nimbleModel(code = nimble.int,  
			    constants = nimble.consts, data = nimble.data, 
			    inits = nimble.inits)
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0',
							'sigma', 'mu', 'fitted', 
							'logProb_richness'), 
                             enableWAIC = TRUE)
nimble.mcmc <- buildMCMC(nimble.conf)
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.int <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                   thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
                   WAIC = TRUE)

# Model with just elevation -----------------------------------------------
nimble.inits <- list(beta.0 = 0, beta.1 = 0, beta.2 = 0, sigma = 1)
nimble.elev <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # Elevation (linear)
  beta.1 ~ dnorm(0, var = 1000)
  # Elevation (quadratic)
  beta.2 ~ dnorm(0, var = 1000)
  # Standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.0 + beta.1 * elev[i] + beta.2 * elev[i]^2
    richness[i] ~ dnorm(mu[i], sd = sigma)
    # Fitted values (necessary for the posterior predictive checks)
    fitted[i] ~ dnorm(mu[i], sd = sigma)
  }
})
nimble.model <- nimbleModel(code = nimble.elev,  
			    constants = nimble.consts, data = nimble.data, 
			    inits = nimble.inits)
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1', 'beta.2',
							'sigma', 'mu', 'fitted', 
							'logProb_richness'), 
                             enableWAIC = TRUE)
nimble.mcmc <- buildMCMC(nimble.conf)
nimble.c.model <- compileNimble(nimble.model)
nimble.c.mcmc <- compileNimble(nimble.mcmc, project = nimble.model)
n.iter <- 30000
n.burn <- 10000
n.thin <- 10
n.chain <- 3
out.elev <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                   thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE, 
                   WAIC = TRUE)

# Posterior predictive checks ---------------------------------------------
# PPCs are a bit more complicated with NIMBLE compared to brms.
fitted.linear <- MCMCchains(out.linear$samples, params = 'fitted', exact = FALSE)
fitted.full <- MCMCchains(out.full$samples, params = 'fitted', exact = FALSE)
fitted.int <- MCMCchains(out.int$samples, params = 'fitted', exact = FALSE)
fitted.elev <- MCMCchains(out.elev$samples, params = 'fitted', exact = FALSE)

# Number of MCMC samples to draw the fitted curve for
n.plot.samples <- 100
# Index of the random n.plot.samples rows to take from the fitted samples 
indx <- sample(1:nrow(fitted.linear), n.plot.samples, replace = FALSE)
# Linear ------------------------------
max.y.val <- max(apply(fitted.linear[indx, ], 2, function(a) density(a)$y), 
               density(dat$bird.rich)$y)
max.x.val <- max(apply(fitted.linear[indx, ], 2, function(a) density(a)$x), 
               density(dat$bird.rich)$x)
min.x.val <- min(apply(fitted.linear[indx, ], 2, function(a) density(a)$x), 
               density(dat$bird.rich)$x)
for (i in 1:n.plot.samples) {
  if (i == 1) {
  plot(density(fitted.linear[indx[i], ]), xlab = 'Bird richness', ylab = 'Density', 
       col = alpha('gray', 0.8), ylim = c(0, max.y.val), main = 'Posterior Predictive Check', 
       xlim = c(min.x.val, max.x.val))    
  } else {
  lines(density(fitted.linear[indx[i], ]), xlab = 'Bird richness', ylab = 'Density', 
       col = alpha('gray', 0.8))    

  }
}
# Add the actual data
lines(density(dat$bird.rich), xlab = 'Bird richness', ylab = 'Density', 
     col = 'black')    
# Full --------------------------------
max.y.val <- max(apply(fitted.full[indx, ], 2, function(a) density(a)$y), 
               density(dat$bird.rich)$y)
max.x.val <- max(apply(fitted.linear[indx, ], 2, function(a) density(a)$x), 
               density(dat$bird.rich)$x)
min.x.val <- min(apply(fitted.linear[indx, ], 2, function(a) density(a)$x), 
               density(dat$bird.rich)$x)
for (i in 1:n.plot.samples) {
  if (i == 1) {
  plot(density(fitted.full[indx[i], ]), xlab = 'Bird richness', ylab = 'Density', 
       col = alpha('gray', 0.8), ylim = c(0, max.y.val), main = 'Posterior Predictive Check', 
       xlim = c(min.x.val, max.x.val))    
  } else {
  lines(density(fitted.full[indx[i], ]), xlab = 'Bird richness', ylab = 'Density', 
       col = alpha('gray', 0.8))    
  }
}
# Add the actual data
lines(density(dat$bird.rich), xlab = 'Bird richness', ylab = 'Density', 
     col = 'black')    
# Can make the same figure from the other two models accordingly.

# Compare the four candidate models using WAIC and LOO --------------------
out.linear$WAIC
out.full$WAIC
out.int$WAIC
out.elev$WAIC
# Calculate LOO using the loo package
log.prob.linear <- MCMCchains(out.linear$samples, params = 'logProb', exact = FALSE)
log.prob.full <- MCMCchains(out.full$samples, params = 'logProb', exact = FALSE)
log.prob.int <- MCMCchains(out.int$samples, params = 'logProb', exact = FALSE)
log.prob.elev <- MCMCchains(out.elev$samples, params = 'logProb', exact = FALSE)
loo.linear <- loo::loo(log.prob.linear)
loo.full <- loo::loo(log.prob.full)
loo.int <- loo::loo(log.prob.int)
loo.elev <- loo::loo(log.prob.elev)

# Plot the resulting model fit --------------------------------------------
# This is again a bit more complicated than using brms as we need to do calculate
# the richness values by hand.
# Old scratch from brms
elev.pred.vals <- seq(from = min(dat$elev),
		      to = max(dat$elev), length.out = 500)
# Predict richness in federal lands across the elevation range
pred.data <- data.frame(elev = elev.pred.vals,
			lo.state = 0, 
			lo.private = 0)
beta.samples.full <- MCMCchains(out.full$samples, params = 'beta', exact = FALSE)
str(beta.samples.full)
pred.plot.samples <- matrix(NA, nrow = nrow(beta.samples.full), ncol = length(elev.pred.vals))
for (i in 1:ncol(pred.plot.samples)) {
  pred.plot.samples[, i] <- beta.samples.full[, 'beta.0'] + 
                            beta.samples.full[, 'beta.1'] * pred.data$elev[i] + 
			    beta.samples.full[, 'beta.2'] * pred.data$elev[i]^2 + 
			    beta.samples.full[, 'beta.3'] * pred.data$lo.state[i] + 
			    beta.samples.full[, 'beta.4'] * pred.data$lo.private[i]
}

plot.dat <- data.frame(elev = elev.pred.vals,
		       est = apply(pred.plot.samples, 2, mean),
		       lower = apply(pred.plot.samples, 2, quantile, 0.025),
		       upper = apply(pred.plot.samples, 2, quantile, 0.975))
ggplot(data = plot.dat, aes(x = elev, y = est)) +
  geom_line(linewidth = 1.2, col = 'gray') +
  geom_line(aes(x = elev, y = lower), linewidth = 0.9, lty = 2, col = 'gray') +
  geom_line(aes(x = elev, y = upper), linewidth = 0.9, lty = 2, col = 'gray') +
  theme_bw(base_size = 24) +
  labs(x = "Elevation", y = "Bird Richness")

