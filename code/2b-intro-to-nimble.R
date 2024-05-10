# 2b-intro-to-nimble.R: this script introduces the NIMBLE package 
#                       for fitting Bayesian models using MCMC. Here we 
#                       showcase how NIMBLE can be used to fit a basic linear
#                       model to understand the relationship between body length
#                       and wingspan in nine blue-eyed hooktail dragonflys. This 
#                       can be compared directly to brms. Take a look at the 
#                       NIMBLE documentation for more details on how to write 
#                       models in NIMBLE: https://r-nimble.org/documentation-2.
# Citations: 
#   Data come from Kéry, M., & Royle, J. A. (2020). 
#                  Applied hierarchical modeling in ecology: Analysis of distribution, 
#                  abundance and species richness in R and BUGS: Volume 2: 
#                  Dynamic and advanced models. Academic Press.
rm(list = ls())
# Load the nimble package to fit the Bayesian model with MCMC 
library(nimble)
# For convergence diagnostics
library(MCMCvis)
library(coda)
# For some plots
library(ggplot2)

# Data prep ---------------------------------------------------------------
# Population that each individual dragonfly comes from
pop <- factor(c(rep('Navarra', 3), rep('Aragon', 3), rep('Catalonia', 3)),
	      levels = c('Navarra', 'Aragon', "Catalonia"))
# Wingspan for each individual
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9)
# Body length for each individual
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2)
# Sex of each individual
sex <- factor(c('M', 'F', 'M', 'F', 'M', 'F', 'M', 'F', 'M'), levels = c('M', 'F'))
# Ectoparasite load (number of mites on each individual)
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6)
# Color intensity (proportion of body that is yellow as opposed to black)
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57)
# Categorical variable indicating how many legs (out of 4) were damaged.
damage <- c(0, 2, 0, 0, 4, 2, 1, 0, 1)

# Put everything together in a data frame
hooktail.df <- data.frame(pop, wing, body, sex, mites, color, damage)

# Number of dragonflies
n <- nrow(hooktail.df)

# Fit the model using NIMBLE ----------------------------------------------
# 1. Specify constants ----------------
# Constant values to give to NIMBLE
nimble.consts <- list(n = n)
# 2. Data -----------------------------
nimble.data <- list(wing = wing, body = body)
# 3. Initial values -------------------
nimble.inits <- list(beta.0 = 0, beta.1 = 0, sigma = 1)
# 4. Write the NIMBLE model -----------
# NOTE: it is fairly common to write this in a separate .R script and then 
#       load the script here. The language used inside the curly brackets
#       uses BUGS syntax (the original syntax that motivated all current
#       Bayesian programming languages).
lm.nimble <- nimbleCode({
  # Priors ----------------------------
  # Intercept
  beta.0 ~ dnorm(0, var = 1000)
  # Body length regression coefficient
  beta.1 ~ dnorm(0, var = 1000) 
  # Standard deviation
  sigma ~ dinvgamma(.001, .001)
  # Likelihood ------------------------
  for (i in 1:n) {
    mu[i] <- beta.0 + beta.1 * body[i]
    # Notice how you can specify different arguments in dnorm depending on 
    # the paramters you are using.
    wing[i] ~ dnorm(mu[i], sd = sigma)
    # Fitted values
    fitted[i] ~ dnorm(mu[i], sd = sigma)
  }
})
# 5. Specify the model ----------------
nimble.model <- nimbleModel(code = lm.nimble, name = 'lmNimble', 
                            constants = nimble.consts, data = nimble.data, 
                            inits = nimble.inits)
# 6. Configure the model --------------
# The monitors argument is what we use to specify the parameters that we want
# to keep track of (i.e., save their MCMC samples). 
nimble.conf <- configureMCMC(nimble.model, monitors = c('beta.0', 'beta.1', 
                                                        'sigma', 'mu', 'fitted'))
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
out.nimble <- runMCMC(nimble.c.mcmc, niter = n.iter, nburnin = n.burn,
                      thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)

# Summarize the model (using MCMCsummary from MCMCvis package)
MCMCsummary(out.nimble)

# Summary of the model fit.
# Useful things to look at in the summary output:
#     1. mean: the posterior mean estimate of the parameter
#     2. sd: standard deviation of the posterior distribution
#     3. 2.5% and 97.5% CI: Lower (2.5%) and upper (97.5%) quantiles of the
#                           95% credible interval.
#     4. Rhat: Gelman-Rubin diagnostic for convergence diagnostics. Should
#              be below 1.1 to indicate your model has converged.
#     5. n.eff: effective sample size estimate, used for
#               assessing model convergence and autocorrelation.

# Convergence Diagnostics -------------------------------------------------
# Are the Rhat values close to 1? How large are the effective sample sizes?
# Plot some traceplots of the resulting posteriors
MCMCtrace(out.nimble)
# Does these plots look like "grass"? Are the three chains falling on top
# of each other?

# Assess autocorrelation to assess chain mixing
acfplot(out.nimble)
# Do the autocorrelation function (acf) plots all eventually hover around 0?

# Explore the model fit ---------------------------------------------------
# Extract the fitted mean values from the model and look at its structure.
mu.samples <- MCMCchains(out.nimble, params = 'mu', exact = FALSE)
str(mu.samples)
# Calculate the mean of the fitted samples
pred.vals <- apply(mu.samples, 2, mean)
str(pred.vals)
# Draw a line with the model fit
plot(body, wing, xlab = 'Body length', ylab = 'Wingspan', las = 1, 
     pch = 21, cex = 1.2, bg = 'lightskyblue1')
lines(x = body, y = pred.vals, col = 'salmon2', lwd = 4)

# Plot the observed vs the predicted --
fitted.samples <- MCMCchains(out.nimble, params = 'fitted', exact = FALSE)
fitted.means <- apply(fitted.samples, 2, mean)
plot(wing, fitted.means, pch = 19, xlab = "Observed", ylab = "Predicted")
# Add a 1-1 line for reference
abline(0, 1)
# Do the predicted values line up with the observed values?
# Are there any systematic differences? (e.g., do we always underpredict
# at high values?)

# Plot the residuals ------------------
res.vals <- wing - fitted.means

# Is there a pattern in the residuals?
plot(fitted.means, res.vals,
     main = "Residuals vs. pred.vals values", las = 1,
     xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)
# Recall what we do when we assess fit of a linear model. When looking
# at a plot of the residual vs. the fitted values, we could have: 
# No pattern -> the assumptions of the model are met.
# Pattern -> assumptions aren't met, which means the results of the model
#            may not be accurate.

