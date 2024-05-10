# 5c-glmm-brms.R: this script displays how to fit a Bayesian GLMM
#                 using brms. For our example, we will look at 
#                 modeling the species richness of an amphibian 
#                 community across some spatial region of interest.
rm(list = ls())
library(ggplot2)
library(brms)

# Overview ----------------------------------------------------------------
# We are interested in modeling species richness of an amphibian community
# across some spatial region of interest. We sample amphibian richness at
# 15 wetlands in each of 10 parks (a total of 150 measurements of richness).
# For each wetland, we obtain two variables that we are interested and think
# will influence amphibian richness: (1) the area of the wetland, and
# (2) the hydroperiod of the wetland. Area is measured as a continuous
# variable, and hydroperiod (the length of time a wetland retains water) is a
# categorical variable measured in three levels: Permanent, Semi-permanent,
# and Temporary. Our goal is to determine how richness is affected by
# area and hydroperiod. Because we have multiple measurements (i.e., wetlands)
# within a given area (i.e., park), we don't think our data meet the
# independence assumption, and so we want to include a random effect of
# "park" in our model.

# Read in the data --------------------------------------------------------
load("data/amphibian-data.rda")
str(amphibian.dat)
# Exploratory data analysis (EDA) -----------------------------------------

# Create an EDA plot (or two) that explores the relationship between
# area of the wetland, hydroperiod, and amphibian species richness. What
# does this tell you about their relationships?

ggplot(data = amphibian.dat, aes(x = area, y = rich, fill = hydroperiod)) +
  geom_point(size = 3, pch = 21) +
  scale_fill_viridis_d() +
  theme_bw(base_size = 14) +
  labs(x = 'Area', y = 'Richness', fill = 'Hydroperiod')

# Positive relationship between area and richness, richness is highest in
# permanent wetlands, followed by semi-permanent wetlands, then temporary wetlands.

# Fit the model -----------------------------------------------------------
# Fit a generalized linear mixed model with a Poisson distribution, random
# effect of park, and fixed effect of hydroperiod and area. Use brm as usual
# and don't standardize the area variable. What happens when you try to fit this model?


# Refit the model using the scale(area) instead of area. What happens now? 
# Does the model converge? Any problems?

# Now fit a GLM with a Poisson distribution without including the random effect
# for park (just include fixed effects of hydroperiod and standardized area).
# Do you notice any differences between the estimates from this model compared to the
# one with a random effect of park?

# Compute the WAIC of the Poisson GLM fit. Compare it to the Poisson GLMM fit. What
# do you see?

# Provide a broad interpretation of the results of the Poisson GLMM based on the summary
# output and the results of the model comparison we did.
