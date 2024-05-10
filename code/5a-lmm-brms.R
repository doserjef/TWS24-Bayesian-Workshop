# 5a-lmm-brms.R: this script displays how to fit a Bayesian LMM using 
#                brms. For our example, we will use a hypothetical data 
#                set on little owls, where we are seeking to understand how
#                wing length varies across populations.
rm(list = ls())
library(brms)
library(ggplot2)

# Load the data and make some EDA plots -----------------------------------
load('data/little-owl-data.rda')
str(little.owl.dat)

boxplot(wing ~ pop, little.owl.dat)
# Does the variation in wing length appear to be larger across populations 
# or within a population? How do you know? 

# Fit a fixed effects model (a linear model) ------------------------------
# Nothing new here. But do notice that here we specify the model in a 
# "means parameterization" by adding the "-1"
out.fixed <- brm(formula = wing ~ factor(pop) - 1,
                 data = little.owl.dat,
                 family = gaussian(),
                 warmup = 500,
                 iter = 1000,
                 chains = 3,
                 init = 'random',
                 cores = 1)
# The model converged, right?


# Fit a random effects model ----------------------------------------------
# We will use the brm() function as usual. To specify random intercepts, our
# notation is (1 | factor), where factor is our categorical variable. This 
# is just like the package lme4 if you have ever fit a frequentist random 
# effects or mixed effects model. If you
# get a warning about a divergent transition, try and increase your burn in
# (warm up) period. 
out.random <- brm(formula = wing ~ (1 | pop),
                  data = little.owl.dat,
                  family = gaussian(),
                  warmup = 500,
                  iter = 1000,
                  chains = 3,
                  init = 'random',
                  cores = 1)
# I received a warning that my Bulk ESS values were too low. If this is the 
# case we should run the model for longer.
out.random.2 <- brm(formula = wing ~ (1 | pop),
                    data = little.owl.dat,
                    family = gaussian(),
                    warmup = 1000,
                    iter = 2000,
                    chains = 3,
                    init = 'random',
                    cores = 1)

# Take a look at the credible intervals from the summary output. Notice that
# the credible interval for the random effect is on the standard deviation.
# That's because this is also what the prior is specified on.
prior_summary(out.random.2)
# The Intercept estimate in the population-level effects section corresponds
# to the mean wing length across all populations.

# Is the standard deviation of the population random effect larger than 
# the residual standard deviation? What does this mean?

# Extract the random effect estimates for each population. Do this using the
# ranef() function. Note that these population effects are all relative to
# the overall intercept.
ranef(out.random.2)

# Next extract the estimated means for each population. To get the estimated
# means for each population, we first need to extract the overall mean from the
# fixed effects. Do this using fixef(out), then add the intercept to the
# population effects you got previously. What is the interpretation of these
# means?
# Extract the intercept samples
# Note the summary = FALSE gives us the full posterior distribution instead
# of just the summary.
intercept.samples <- fixef(out.random.2, summary = FALSE)
# Extract the random effects for each population
random.effect.samples <- ranef(out.random.2, summary = FALSE)
str(random.effect.samples)
# Simplify the structure a bit
random.effect.samples <- random.effect.samples$pop[, , 1]
# Create the full wing lengths for each population
pop.wing.length.samples <- random.effect.samples + intercept.samples[, 1]

# Generate a boxplot showing the estimates --------------------------------
plot.df <- data.frame(median = apply(pop.wing.length.samples, 2, median), 
                      lowest = apply(pop.wing.length.samples, 2, quantile, 0.025), 
                      low = apply(pop.wing.length.samples, 2, quantile, 0.25),
                      high = apply(pop.wing.length.samples, 2, quantile, 0.75), 
                      highest = apply(pop.wing.length.samples, 2, quantile, 0.975), 
                      population = levels(little.owl.dat$pop))

ggplot(data = plot.df, aes(x = population)) + 
  geom_boxplot(aes(ymin = lowest, lower = low, middle = median, 
                   upper = high, ymax = highest), stat = 'identity', 
               fill = 'lightskyblue1') + 
  theme_bw(base_size = 16) + 
  labs(x = 'Population', y = 'Little owl wing length')