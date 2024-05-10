# 9-N-mixture-models.R: this script showcases how to fit an N-mixture model 
#                       in spAbundance. The example data set uses point count
#                       data from the Hubbard Brook Experimental Forest. 
# The data are obtained from the spAbundance R package. See ?hbefCount2015 for 
# more detailed information. 
# Note there is no associated lecture notes for this exercise.
# For more resources on fitting N-mixture models in spAbundance, see the 
# tutorial at the following link, where I give much more 
# details on how to specify all of the arguments to fit 
# a Bayesian N-mixture model: https://www.jeffdoser.com/files/spabundance-web/articles/nmixturemodels.

rm(list = ls())
library(spAbundance)
library(ggplot2)
library(stars)

# Extract the data from one species ---------------------------------------
# hbefCount2015 contains count data for 12 foliage gleaning bird species
data.hbef <- hbefCount2015
# Let's work with Blackburnian Warbler
data.hbef$y <- data.hbef$y[which(dimnames(hbefCount2015$y)[[1]] == 'BLBW'), , ]
# Take a look at the data structure for a single-species N-mixture model.
str(data.hbef)
# spAbundance requires the data be provided in a specific format. The 
# data are stored in a list with four components: 
#   1. y: the repeated count data. The rows correspond to sites, columns 
#         correspond to the repat visits
#   2. abund.covs: a data frame or matrix consisting of the covariates that you 
#                  want to include in the model for abundance.
#   3. det.covs: a list of variables that you want to include in the model for 
#                modeling detection probability. Each component of the list is a 
#                different variable, which can either be a vector for site-level
#                variables or a site x visit matrix for observational level 
#                variables (which is what we have here).
#   4. coords: a matrix where the rows correspond to the sites and the two 
#              columns correspond to the easting and northing of the given site.

# Fit the model -----------------------------------------------------------
# In spAbundance, the total number of MCMC samples is split up into a series
# of MCMC batches, each with a specified number of samples. The total
# number of samples is the number of batches times the length of each batch.
# Number of MCMC batches
n.batch <- 4000
# Length of each batch
batch.length <- 25
# Total number of samples per chain is n.batch * batch.length
n.batch * batch.length
n.burn <- 50000
n.thin <- 20
n.chains <- 3
# Approx run time: 3 min
out <- NMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
            det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod),
            data = data.hbef,
            n.batch = n.batch,
            batch.length = batch.length,
            family = 'Poisson',
            accept.rate = 0.43,
            n.omp.threads = 1,
            verbose = TRUE,
            n.report = 100,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains)

# Quick summary of the model results
summary(out)
# Assess convergence 
plot(out, param = 'beta')

# What is the probability there is a negative linear effect of day of the year
# on detection probability?
apply(out$alpha.samples, 2, function(a) mean(a < 0))
# The above line of code gives the probability of each of the detection 
# parameters being less than 0.

# Posterior predictive check ----------------------------------------------
# Perform a numeric posterior predictive check to calculate a Bayesian 
# p-value. 
ppc.out <- ppcAbund(out, fit.stat = 'freeman-tukey', group = 1)
# Report the Bayesian p-value
summary(ppc.out)
# The Bayesian p-value is large, suggesting our model may not be a great
# fit for the data. Let's explore this further with a graphical check. This
# is a bit clunky at the moment, but I am hoping at some point to automate
# this within spAbundance like what is done in brms.
# Graphically plot the fitted data vs. the true data
# Get the fitted values
fitted.out <- fitted(out)
# Only use a subset of MCMC samples to make the resulting objects a bit 
# smaller
# Number of MCMC samples to draw the fitted curve for
n.plot.samples <- 100
# Index of the random n.plot.samples rows to take from the fitted samples 
indx <- sample(1:nrow(fitted.out$y.rep.samples), n.plot.samples, replace = FALSE)
y.rep.samples <- fitted.out$y.rep.samples[indx, , ]
# Number of data points
n.data <- sum(!is.na(data.hbef$y))
fitted.mat <- matrix(NA, n.plot.samples, n.data)
for (i in 1:n.plot.samples) {
  fitted.mat[i, ] <- c(y.rep.samples[i, , ])[which(!is.na(c(y.rep.samples[i, , ])))]
}
y.long <- data.hbef$y[which(!is.na(c(data.hbef$y)))]
# Linear ------------------------------
max.y.val <- max(apply(fitted.mat, 2, function(a) density(a)$y), 
		 density(y.long)$y)
max.x.val <- max(apply(fitted.mat, 2, function(a) density(a)$x), 
		 density(y.long)$x)
min.x.val <- min(apply(fitted.mat, 2, function(a) density(a)$x), 
		 density(y.long)$x)
for (i in 1:n.plot.samples) {
  if (i == 1) {
  plot(density(fitted.mat[i, ]), xlab = '# BLBW observed', ylab = 'Density',
       col = alpha('gray', 0.8), ylim = c(0, max.y.val), main = 'Posterior Predictive Check',
       xlim = c(min.x.val, max.x.val))
  } else {
  lines(density(fitted.mat[i, ]), xlab = '# BLBW observed', ylab = 'Density',
       col = alpha('gray', 0.8))
  }
}
# Add the actual data
lines(density(y.long), xlab = '# BLBW observed', ylab = 'Density',
     col = 'black')
# The visual PPC looks pretty good, despite the large p-value above, 
# showcasing the sensitivity of Bayesian p-values to single values.

# Summary figures ---------------------------------------------------------
# Generate a conditional effects plot
elev.pred.vals <- seq(min(data.hbef$abund.covs$elev), 
                      max(data.hbef$abund.covs$elev), length.out = 200)
# Scale predicted values by mean and sd used to fit the model
elev.pred.vals.s <- (elev.pred.vals - mean(data.hbef$abund.covs$elev)) /
                    sd(data.hbef$abund.covs$elev)
X.0 <- as.matrix(data.frame(intercept = 1, 
                            elev = elev.pred.vals.s, 
                            elev.2 = elev.pred.vals.s^2))
out.pred <- predict(out, X.0)
str(out.pred)
# Get the lower bound, median, and 95% credible interval
mu.0.quants <- apply(out.pred$mu.0.samples, 2, quantile, 
                     prob = c(0.025, 0.5, 0.975))
mu.plot.dat <- data.frame(mu.med = mu.0.quants[2, ],
                          mu.low = mu.0.quants[1, ],
                          mu.high = mu.0.quants[3, ],
                          elev = elev.pred.vals)
ggplot(mu.plot.dat, aes(x = elev, y = mu.med)) +
  geom_ribbon(aes(ymin = mu.low, ymax = mu.high), fill = 'grey70', alpha = 0.5) +
  geom_line() +
  theme_bw() +
  labs(x = 'Elevation', y = 'Blackburnian warbler abundance')