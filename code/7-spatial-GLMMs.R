# 7-spatial-GLMMs.R: this script showcases how to fit a spatial 
#                    GLMM to predict Hooded Warbler relative abundance
#                    across Pennsylvania. 
# Data citation:
#      These are real data that come from the USGS Breeding Bird Survey, and
#      we load them through the spAbundance package. After loading spAbundance,
#      take a look at ?bbsData for the citation information.
rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
library(stars)

# NOTE: more detailed information on fitting spatial GLMMs in the 
#       spAbundance R package is available at: 
#       https://www.jeffdoser.com/files/spabundance-web/articles/glmm

# Load the data -----------------------------------------------------------
data(bbsData)
# This is a list, formatted for fitting Bayesian models in spAbundance.
str(bbsData)
sp.names <- dimnames(bbsData$y)[[1]]
data.HOWA <- bbsData
data.HOWA$y <- data.HOWA$y[which(sp.names == "HOWA"), ]
# Observed number of HOWA at each site
data.HOWA$y
str(data.HOWA)
# This data list is in the format required for fitting spatial GLMMs in 
# spAbundance. It is a list with the following components: 
#    1. y: the resopnse variable
#    2. covs: a data frame or matrix of the covariates you may want to 
#             include in the model.
#    3. coords: a matrix of spatial coordinates for the sites in the data set.
str(data.list)

# Fit the model -----------------------------------------------------------
# In spAbundance, the total number of MCMC samples is split up into a series
# of MCMC batches, each with a specified number of samples. The total
# number of samples is the number of batches times the length of each batch.
n.batch <- 1000
batch.length <- 25
n.burn <- 15000
n.thin <- 10
n.chains <- 3
data.list$y <- sqrt(data.list$y)

out <- spAbund(formula = ~ scale(bio2) + scale(bio8) + scale(bio18) + scale(forest) +
                           scale(devel) + scale(day) + I(scale(day)^2) + scale(tod) +
                           (1 | obs),
               data = data.HOWA, n.neighbors = 15, 
               cov.model = 'exponential', NNGP = TRUE, family = 'Poisson',
               n.batch = n.batch, batch.length = batch.length,
               n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
               n.chains = n.chains, n.report = 100, n.omp.threads = 1)

summary(out)
# Regression coefficients
plot(out, 'beta', density = FALSE)
# Random effect variance for observer variability
plot(out, 'sigma.sq.mu', density = FALSE)
# Spatial parameters
plot(out, 'theta', density = FALSE)

# What variables have the strongest influence on the relative abundance of 
# Hooded Warblers across the state? 

# Predict relative abundance across the state -----------------------------
data(bbsPredData)
# Center and scale covariates by values used to fit model
bio2.pred <- (bbsPredData$bio2 - mean(data.HOWA$covs$bio2)) /
              sd(data.HOWA$covs$bio2)
bio8.pred <- (bbsPredData$bio8 - mean(data.HOWA$covs$bio8)) /
              sd(data.HOWA$covs$bio8)
bio18.pred <- (bbsPredData$bio18 - mean(data.HOWA$covs$bio18)) /
               sd(data.HOWA$covs$bio18)
forest.pred <- (bbsPredData$forest - mean(data.HOWA$covs$forest)) /
                sd(data.HOWA$covs$forest)
devel.pred <- (bbsPredData$devel - mean(data.HOWA$covs$devel)) /
               sd(data.HOWA$covs$devel)
day.pred <- 0
tod.pred <- 0
X.0 <- cbind(1, bio2.pred, bio8.pred, bio18.pred, forest.pred,
             devel.pred, day.pred, day.pred^2, tod.pred)
colnames(X.0) <- c('(Intercept)', 'scale(bio2)', 'scale(bio8)', 'scale(bio18)',
                   'scale(forest)', 'scale(devel)', 'scale(day)',
                   'I(scale(day)^2)', 'scale(tod)')
coords.0 <- as.matrix(bbsPredData[, c('x', 'y')])
out.sp.pred <- predict(out, X.0, coords.0, ignore.RE = TRUE, verbose = TRUE)
mu.0.quants <- apply(out.sp.pred$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot.df <- data.frame(Easting = bbsPredData$x,
                      Northing = bbsPredData$y,
                      mu.0.med = mu.0.quants[2, ],
                      mu.0.ci.width = mu.0.quants[3, ] - mu.0.quants[1, ])
# proj4string for the coordinate reference system
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
coords.stars <- st_as_stars(plot.df, crs = my.crs)
coords.sf <- st_as_sf(as.data.frame(data.HOWA$coords), coords = c('X', 'Y'),
                      crs = my.crs)
# Plot of median estimate
ggplot() +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.med)) +
  geom_sf(data = coords.sf, col = 'grey') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 12) +
  labs(fill = '', x = 'Longitude', y = 'Latitude',
       title = 'Hooded Warbler Relative Abundance')