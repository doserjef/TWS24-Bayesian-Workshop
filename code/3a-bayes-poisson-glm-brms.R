# 3a-bayes-poisson-glm-brms.R: this script displays how to fit a Bayesian 
#                              generalized linear model using brms. For our
#                              example, we will look at how hooded warbler
#                              relative abundance varies across Pennsylvania.
# Data citation: 
#      These are real data that come from the USGS Breeding Bird Survey, and 
#      we load them through the spAbundance package. After loading spAbundance,
#      take a look at ?bbsData for the citation information.
rm(list = ls())
library(brms)
library(ggplot2)
# Load spAbundance to get the hooded warbler data 
library(spAbundance)
# For some maps
library(sf)
library(stars)

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
# NOTE: when specifying the family, we explicitly set the link = 'log'. 
# See ?brmsfamily for details on the different link functions you can use.
# NOTE: below we set cores = 3, which allows us to use parallelization to run
#       the three chains of the MCMC algorithm in parallel.
# NOTE: we are using the scale() function to standardize the climate variables.
#       Standardizing variables on vastly different scales can be pretty important
#       for ensuring model convergence.
out.clim <- brm(formula = howa ~ scale(bio2) + scale(bio8) + scale(bio18) + 
                                 scale(day) + scale(tod),
                  data = bbs.dat,
                  family = poisson(link = 'log'),
                  warmup = 500,
                  iter = 1000,
                  chains = 3,
                  init = 'random',
                  cores = 3)

# Has the model converged? 
summary(out.clim)

# What priors are used?
prior_summary(out.clim)

# Assess model fit with some posterior predictive check figures.
pp_check(out.clim, ndraws = 100)

# Oof, that doesn't look to good. We will ignore this for now, but what do you 
# think we could do to try and improve upon this? What other distributions work well
# with count data?

# Fit land cover only model -----------------------------------------------
out.land <- brm(formula = howa ~ scale(forest) + scale(devel) +
                                 scale(day) + scale(tod),
                data = bbs.dat,
                family = poisson(link = 'log'),
                warmup = 500,
                iter = 1000,
                chains = 3,
                init = 'random',
                cores = 3)
# Has the model converged?
summary(out.land)
# What priors are used? 

# Fit climate + land cover model ------------------------------------------
out.full <- brm(formula = howa ~ scale(forest) + scale(devel) + scale(bio2) + 
                                 scale(bio8) + scale(bio18) + 
                                 scale(day) + scale(tod),
                data = bbs.dat,
                family = poisson(link = 'log'),
                warmup = 500,
                iter = 1000,
                chains = 3,
                init = 'random',
                cores = 3)
# Has the model converged?
summary(out.full)
# What priors are used? 
prior_summary(out.full)
# Assess model fit with some posterior predictive check figures.

# Compare the three candidate models --------------------------------------
out.clim <- add_criterion(out.clim, c('waic', 'loo'))
out.land <- add_criterion(out.land, c('waic', 'loo'))
out.full <- add_criterion(out.full, c('waic', 'loo'))
# Compare all three models with WAIC
loo_compare(out.full, out.clim, out.land, criterion = c('waic'))
# Compare all three models with LOO
loo_compare(out.full, out.clim, out.land, criterion = c('loo'))

# Do a bit of interpretation for the top performing model -----------------
# Which is the top supported model? 
# out.land
# What is the probability of a positive effect of forest cover?
# What is the probability the effect of forest cover is greater than 
# the effect of time of day?

# Visualize relationships in the top performing model ---------------------
top.model.effects <- conditional_effects(out.land)
top.model.effects
# Can use ggplot2 to make a prettier plot.

# Predict relative abundance across Pennsylvania --------------------------
# bbsPredData contains all covariate values across a 12 x 12km grid of PA. 
# Let's predict relative abundance of hooded warbler across the state.
data(bbsPredData)
str(bbsPredData)
# Create data frame needed for brms
bbs.pred.df <- data.frame(forest = bbsPredData$forest, 
			  devel = bbsPredData$devel,
			  day = mean(bbs.dat$day), 
			  tod = mean(bbs.dat$tod))
# Predict relative abundance at the 
out.pred <- predict(out.land, bbs.pred.df)

plot.df <- data.frame(Easting = bbsPredData$x,
                      Northing = bbsPredData$y,
                      mu.0.med = out.pred[, 'Estimate'],
                      mu.0.ci.width = out.pred[, 'Q97.5'] - out.pred[, 'Q2.5'])
# proj4string for the coordinate reference system
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
coords.stars <- st_as_stars(plot.df, crs = my.crs)
coords.sf <- st_as_sf(as.data.frame(bbsData$coords), coords = c('X', 'Y'),
                      crs = my.crs)
# Plot of median estimate
ggplot() +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.med)) +
  geom_sf(data = coords.sf, col = 'grey') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 12) +
  labs(fill = '', x = 'Longitude', y = 'Latitude',
       title = 'Hooded Warbler Relative Abundance')
