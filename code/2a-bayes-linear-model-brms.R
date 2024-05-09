# 2a-bayes-linear-model-brms.R: this script displays how to fit a Bayesian linear
#                               model using brms. For our example, we will look at 
#                               how bird species richness varies in relation to 
#                               elevation and landowner type (which may result in 
#                               different management practices).
rm(list = ls())
library(brms)
library(ggplot2)

# Load data and some EDA --------------------------------------------------
# Change working directory as needed
load("data/bird-elev-data.rda")
str(dat)
# Bird richness relationship with elevation
ggplot(data = dat, aes(x = elev, y = bird.rich)) +
  geom_point(size = 3, pch = 21, fill = 'lightskyblue') +
  theme_bw(base_size = 18) +
  labs(x = 'Elevation', y = 'Bird Richness')

# Bird richness relationship with landowner type
boxplot(bird.rich ~ land.owner, data = dat)

# Fit a simple linear model -----------------------------------------------
# Fit a model with bird richness as the response variable and a linear 
# effect of elevation and land owner. Note that the elevation variable is 
# standardized. 
out.linear <- brm(formula = bird.rich ~ elev + land.owner,
                  data = dat,
                  family = gaussian(),
                  warmup = 500,
                  iter = 1000,
                  chains = 3,
                  init = 'random',
                  cores = 1)

# Has the model converged? Use both visual assessments as well as 
# numeric criteria to support your decision.
summary(out.linear)
mcmc_plot(out.linear, type = 'trace')

# However, before we make any conclusions we should always check the residuals
res.vals <- residuals(out.linear)
pred.vals <- fitted(out.linear)

plot(pred.vals[, 1], res.vals[, 1], main = "Residuals vs. pred.vals values",
     las = 1, xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)

# What does this plot suggest to you? Can we continue analyzing this model?

# It could also be fun to plot the normal probability plot (qqplot). 
qqnorm(res.vals[, 1], pch = 19)
qqline(res.vals[, 1], col = 'red')

# Fit a model with a linear and quadratic term for elevation --------------
# Fit the same model as before, but now add in a quadratic effect of elevation
out.full <- brm(formula = bird.rich ~ elev + I(elev^2) + land.owner,
	        data = dat,
	        family = gaussian(),
	        warmup = 500,
	        iter = 1000,
	        chains = 3,
	        init = 'random',
	        cores = 1)

# Has the model converged? Perform the usual tests.
summary(out.full)

# Now let's check the residuals of this plot
res.vals.full <- residuals(out.full)
pred.vals.full <- fitted(out.full)

plot(pred.vals.full[, 1], res.vals.full[, 1], main = "Residuals vs. pred.vals values",
     las = 1, xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)

# What do these residuals look like? Better than the previous? Or not?
# What about the normal probability plot?

# Fit an intercept only model ---------------------------------------------
# Next for comparison, fit a model with just an intercept
out.int <- brm(formula = bird.rich ~ 1,
	        data = dat,
	        family = gaussian(),
	        warmup = 500,
	        iter = 1000,
	        chains = 3,
	        init = 'random',
	        cores = 1)
# Has the model converged? What do the residuals look like?

# Fit a model with just elevation -----------------------------------------
# Now, fit a model without the landowner effect
out.elev <- brm(formula = bird.rich ~ elev + I(elev^2),
	        data = dat,
	        family = gaussian(),
	        warmup = 500,
	        iter = 1000,
	        chains = 3,
	        init = 'random',
	        cores = 1,
	        seed = 123)

# Posterior predictive checks ---------------------------------------------
# Plots of the true data density compared to model-generated data densities
# for three of our candidate models
pp_check(out.int, ndraws = 100)
pp_check(out.linear, ndraws = 100)
pp_check(out.full, ndraws = 100)

# Histogram showing the standard deviation of the data points for the true 
# data set compared to the distribution of the standard deviations for a 
# variety of data types.
pp_check(out.int, ndraws = 100, type = 'stat', stat = 'sd')
pp_check(out.linear, ndraws = 100, type = 'stat', stat = 'sd')
pp_check(out.full, ndraws = 100, type = 'stat', stat = 'sd')

# Plot of the true vs fitted data values 
pp_check(out.int, ndraws = 100, type = 'scatter_avg')
pp_check(out.linear, ndraws = 100, type = 'scatter_avg')
pp_check(out.full, ndraws = 100, type = 'scatter_avg')

# Compare the four candidate models using WAIC and LOO --------------------
out.linear <- add_criterion(out.linear, c('waic', 'loo'))
out.full <- add_criterion(out.full, c('waic', 'loo'))
out.int <- add_criterion(out.int, c('waic', 'loo'))
out.elev <- add_criterion(out.elev, c('waic', 'loo'))
# Compare all three models with WAIC
loo_compare(out.int, out.linear, out.elev, out.full, criterion = c('waic'))
# Compare all three models with LOO
loo_compare(out.int, out.linear, out.elev, out.full, criterion = c('loo'))

# Plot the resulting model fit --------------------------------------------
# Display conditional effects plots using brms
elev.cond.effects <- conditional_effects(out.full, effects = 'elev')
# You can look at the default plot that brms shows
elev.cond.effects
# Or you can make your own and customize it as you desire with ggplot
str(elev.cond.effects)

ggplot(data = elev.cond.effects$elev, aes(x = elev, y = estimate__)) + 
  geom_line(linewidth = 1.2, col = 'gray') +
  geom_line(aes(x = elev, y = lower__), linewidth = 0.9, lty = 2, col = 'gray') +
  geom_line(aes(x = elev, y = upper__), linewidth = 0.9, lty = 2, col = 'gray') +
  theme_bw(base_size = 24) +
  labs(x = "Elevation", y = "Bird Richness")
