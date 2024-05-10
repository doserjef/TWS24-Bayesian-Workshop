# 2a-intro-to-brms.R: this script introduces the brms package that we 
#                     will use for fitting a variety of Bayesian models.
#                     Here we will specifically focus on a basic example
#                     where we want to understand the linear relationship
#                     between body length and wingspan in nine blue-eyed 
#                     hooktail dragonflys.
# Citations: 
#   Data come from Kéry, M., & Royle, J. A. (2020). 
#                  Applied hierarchical modeling in ecology: Analysis of distribution, 
#                  abundance and species richness in R and BUGS: Volume 2: 
#                  Dynamic and advanced models. Academic Press.

# Uncomment the below code to clear out your workspace. I normally put this
# at the top of all of my R scripts to make sure I keep a tidy workspace.
rm(list = ls())
# Load the brms package to fit a Bayesian model.
library(brms)
# For some plots
library(ggplot2)
# Code tip: it's a good idea to load all packages at the top of your R script.
#           This is useful for you to remind yourself what packages you need,
#           as well as for other people if you end up sharing your code.


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
str(hooktail.df)

# Fit the model using brms ------------------------------------------------
# In brms, the "brm" function is what we use to fit different models. By default, 
# the function will specify initial values and prior distributions that are 
# reasonably vague.
out <- brm(formula = wing ~ body,
           data = hooktail.df,
           family = gaussian(),
           warmup = 500,
           iter = 1000,
           chains = 3,
           init = 'random',
           cores = 1,
           seed = 123)

# Summary of the model fit.
summary(out)
# Useful things to look at in the summary output:
#     1. Estimate: the posterior mean estimate of the parameter
#     2. Est. Error: standard deviation of the posterior distribution
#     3. l-95%CI, u-95% CI: Lower (2.5%) and upper (97.5%) quantiles of the
#                           95% credible interval.
#     4. Rhat: Gelman-Rubin diagnostic for convergence diagnostics. Should
#              be below 1.1 to indicate your model has converged.
#     5. BULK_ESS, Tail_ESS: effective sample size estimates, used for
#                            assessing model convergence and autocorrelation.

# Convergence Diagnostics -------------------------------------------------
# Are the Rhat values close to 1? How large are the effective sample sizes?
# Plot some traceplots of the resulting posteriors
mcmc_plot(out, type = 'trace')
# Do these plots look like "grass"? Are the three chains falling on top
# of each other?

# Assess autocorrelation to assess chain mixing
mcmc_plot(out, type = 'acf')
# Do the autocorrelation function (acf) plots all eventually hover around 0?

# Some basic manipulations with the posterior samples ---------------------
# Extract the MCMC samples as a data frame. See ?as_draws for more details
out.samples <- as_draws_df(out)
# Each column is a parameter, rows are the MCMC samples (draws) for each chain
str(out.samples)
# What is the probability the effect of body is positive? 
mean(out.samples$b_body > 0)
paste0("The probability of a positive relationship of body length and wingspan is ", 
       round(mean(out.samples$b_body > 0), 2))

# Explore the model fit ---------------------------------------------------
# Extract the fitted values from the model and look at its structure.
fitted.vals <- fitted(out)
str(fitted.vals)
# Draw a line with the model fit
plot(body, wing, xlab = 'Body length', ylab = 'Wingspan', las = 1, 
     pch = 21, cex = 1.2, bg = 'lightskyblue1')
lines(x = body, y = fitted.vals[, "Estimate"], col = 'salmon2', lwd = 4)

# Plot the observed vs the predicted --
plot(wing, fitted.vals[, "Estimate"], pch = 19, xlab = "Observed", ylab = "Predicted")
# Add a 1-1 line for reference
abline(0, 1)
# Do the predicted values line up with the observed values?

# Plot the residuals ------------------
res.vals <- residuals(out)
# OR
# res.vals <- y - fitted.vals[, "Estimate"]

# Is there a pattern in the residuals?
plot(fitted.vals[, "Estimate"], res.vals[, "Estimate"],
     main = "Residuals vs. fitted.vals values", las = 1,
     xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)
# Recall what we do when we assess fit of a linear model. When looking
# at a plot of the residual vs. the fitted values, we could have: 
# No pattern -> the assumptions of the model are met.
# Pattern -> assumptions aren't met, which means the results of the model
#            may not be accurate.
# How do you feel about this model?

# Comparing models with different prior distributions ---------------------
# Look at the default priors used when we don't specify them when 
# running brm()
prior_summary(out)

# Let's set our own priors based on those we talked about in the lecture
# Use a Normal distribution with mean 0 and variance 1000 for the intercept
# and slope. Use an inverse-Gamma for the standard deviation
# NOTE: a key point of confusion when using different Bayesian software 
#       is that some software uses the variance parameter in the Normal 
#       distribution, some use the standard deviation (like rnorm), and some
#       even use the "precision" (the inverse of the variance). brms uses 
#       the standard deviation.
# Prior for the intercept
prior.intercept <- prior(normal(1, sqrt(1000)), class = 'Intercept')
# Prior for the regression coefficient on body length
prior.body <- prior(normal(1, sqrt(1000)), class = b, coef = 'body')
# Prior for the residual standard deviation
prior.resid.sd <- prior(inv_gamma(.001, .001), class = 'sigma')
# Combine all priors together
priors.all <- c(prior.intercept, prior.body, prior.resid.sd)
# Run the model, and explicitly specify the priors in the "prior" argument
out.2 <- brm(formula = wing ~ body,
             data = hooktail.df,
             family = gaussian(),
             prior = priors.all,
             warmup = 500,
             iter = 1000,
             chains = 3,
             init = 'random',
             cores = 1)

# Are there large differences between the two models that use different priors? 

# Do you think that increasing the sample size of our data would lead to more 
# or less differences when fitting models with different priors?

# Compare the estimates to lm() -------------------------------------------
out.lm <- lm(formula = wing ~ body, data = hooktail.df)
summary(out.lm)
# P-value is around 0.28. How to interpret that? Is that more intuitive than 
# our direct probability statement using our Bayesian model for the effect of 
# body length on wingspan? 

# Compare the fitted values for the frequentist and Bayes model
fitted.vals.lm <- fitted(out.lm)
plot(fitted.vals[, 'Estimate'], fitted.vals.lm, pch = 19, 
     xlab = 'Bayesian', ylab = 'Frequentist')
abline(0, 1)
cor(fitted.vals[, 'Estimate'], fitted.vals.lm)
