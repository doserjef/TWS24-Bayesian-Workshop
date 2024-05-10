# 4d-binomial-glm-nimble.R: in this script, you will fit a Bayesian 
#                           generalized linear model using NIMBLE. For our 
#                           example here, we will look at the distribution 
#                           of eastern hemlock across Vewrmont and how it 
#                           relates to a variety of climate variables.
# Data citation: 
#   These data come from the most recent survey of the USFS Forest Inventory 
#   and Analysis data set:
#     Bechtold, W. A. and Patterson, P. L. (2005). The enhanced Forest Inventory and Analysis 
#     Program–National sampling design and estimation procedures. Number 80. USDA
#     Forest Service, Southern Research Station. 
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

# Load the data and make a few EDA plots ----------------------------------
load('data/eastern-hemlock-data.rda')
# y is the presence or absence of eastern hemlock at each of the inventory plots.
str(hemlock.df)

# CRS of the data
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

hemlock.sf <- st_as_sf(hemlock.df, coords = c('easting', 'northing'), crs = my.crs)
# Plot presence/absence of the species across the state.
hemlock.sf$val <- factor(ifelse(hemlock.sf$y == 1, 'Present', 'Absent'))
ggplot() +
  geom_sf(data = hemlock.sf, aes(col = val), size = 2) +
  scale_color_viridis_d() +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Longitude", y = "Latitude", col = "")
# Note that these aren't actually the true locations of the forest plots. There is 
# some noise that is added to the true plot locations to preserve plot integrity 
# and land owner privacy.

# Now let's take a look at the covariates
covs <- hemlock.df[, -c(1:3)]
# tmax = max monthly temperature
# tmin = min monthly temperature
# prcp = total annual precipitation
# pet = total growing season potential evapotranspiration 
#       (equal to the total water that would be lost from a 
#        forest ecosystem given local climatic conditions with an unlimited supply of water)
# aet = total growing season actual evapotranspiration 
#       (equal to the total water actually lost from a 
#        forest ecosystem given local climatic conditions and 
#        water availability—if water is not limiting, aet = pet)
# water_deficit = total growing season climatic water deficit (equal to the 
#                 difference between pet and aet—positive values indicate 
#                 that there an insufficient amount of water to meet atmospheric water demand of photosynthesis)
# vpd = mean growing season vapor pressure deficit (equal to the difference between the partial water pressure in the atmosphere and the saturation water pressure—it is a measure of the atmospheric water demand—higher vpd means greater water loss during photosynthesis due to dryer atmospheric conditions)
# elev = elevation
cov.names <- colnames(covs)
p.cov <- length(cov.names)
# Standardize covariates to make the plotting scales look nice
covs.std <- apply(covs, 2, scale)
plot.df <- data.frame(val = c(covs.std),
                      easting = rep(hemlock.df$easting, times = p.cov),
                      northing = rep(hemlock.df$northing, times = p.cov),
                      covariate = factor(rep(cov.names, each = nrow(hemlock.df))))
plot.sf <- st_as_sf(plot.df, coords = c('easting', 'northing'), crs = my.crs)
ggplot() +
  geom_sf(data = plot.sf, aes(col = val), size = 2) +
  scale_color_viridis_c(option = 'plasma') +
  theme_bw(base_size = 10) +
  facet_wrap(~covariate, nrow = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Longitude", y = "Latitude", col = "")

# Based on the patterns observed in the climate variables and the variables, 
# which climate variables may be related to the output? 

# How correlated are the covariates?
pairs(covs, upper.panel = NULL, pch = 19, cex = 0.5)
cor(covs)
# Very!

# Your task is to fit and compare three models that predict the distribution 
# of eastern hemlock across the state: 
# Model 1: A model with just a linear effect of elevation
# Model 2:  A model with linear effects of minimum temperature, actual evapotranspiration, 
#     and elevation
# Model 3: A model with linear and quadratic effects of minimum temperature, 
#     actual evapotranspiration, and elevation

# After fitting the models, you will need to assess convergence, perform 
# a posterior predictive check, do some model comparison, and determine which 
# model is the best of the three candidate models. You should be able to adapt
# code from the previous examples 

# Fit the three candidate models ------------------------------------------

# Convergence assessment --------------------------------------------------

# Posterior predictive checks ---------------------------------------------

# Model comparison --------------------------------------------------------

# Conditional effects plots -----------------------------------------------


