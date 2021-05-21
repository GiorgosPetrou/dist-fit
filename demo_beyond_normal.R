# Load Libraries ------
library(dplyr)
library(ggplot2)
source("./multi_fit.R")
set.seed(789)

# Load Data ------------
filled_cavity_df <- read.csv("./data.csv", header = TRUE)

# Step 1: Histogram  ----

# Estimating binwidth based on Freedman-Diaconis
bw <- 2 * IQR(filled_cavity_df$Corrected.Measured.Mean) /
  length(filled_cavity_df$Corrected.Measured.Mean)^(1 / 3)

# Plot histogram with density line
ggplot(data = filled_cavity_df, aes(x = Corrected.Measured.Mean)) +
  geom_histogram(aes(y=..density..),binwidth = bw, fill = "grey", colour = "black") +
  geom_density() + xlab(expression(paste("Mean U-value [W/(",m^2,"K)]",sep = ""))) + 
  ylab("Density") + theme_bw()

# Step 2: Cullen and Frey ----
descdist(filled_cavity_df$Corrected.Measured.Mean, boot = 1000)

# Step 3: Distribution Fitting ------

# Define the distributions to be fitted

fits <- multi_fit(data = filled_cavity_df$Corrected.Measured.Mean, dists = c("norm", "lnorm", "gamma", "weibull"))

# Step 4a: AIC ----

gof_summary(fits, params = TRUE)


# Step 4b: GOF Plots ----
plot_GOF(fits = fits, dists = c("gamma", "norm"), legend_items = c("gamma", "normal"), x_label = "Wall U-value")

# Best Fit Parameters ----
print(fits$gamma)

