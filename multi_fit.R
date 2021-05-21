# Load Libraries -----
library(fitdistrplus)
library(viridis)

# Fitting Function -----

# Fit multiple distributions
# Input:
#   Data: Vector of data to be fitted
#   dists: Vector of distribution names to be fitted
# Output:
#   List of fitted distributions
multi_fit <- function(data, dists = c("norm", "lnorm", "gamma", "weibull")){
  # Create a list of fitted distributions
  fits <- list()
  for (i in dists){
    fits[[i]] <- fitdist(data, i)
  }
  return(fits)
}

# AIC Functions -----

# Estimate corrected AIC (AICc)
# Input:
#   fit: A fitdist object
# Output:
#   AICc
estimate_aicc <- 
  function(fit){
    p <- length(fit$estimate)
    n <- fit$n
    LL <- fit$loglik
    aic <- -2*LL + 2*p
    sec_order_term <- (2*p+1)/(n - p -1)
    aicc <- aic + sec_order_term
    return(aicc)
  }

# Estimate AIC depending on sample and parameter size
# and order distributions according to their AIC
# Input:
#  fits: multi_fit object - a list of fitdist objects
# Output:
#  fits_aic_ordered: Dataframe of AIC ordered low to high
estimate_aic <- function(fits){
  no_fits <- length(fits)
  fits_aic <- data.frame("Distributions" = rep(NA,no_fits), 
                         "AIC" = rep(NA,no_fits))
  for (i in 1:no_fits){
    fits_aic[i,"Distributions"] <- names(fits)[i]
    k <- length(fits[[i]]$estimate) #  no. of fitted parameters
    n <- fits[[i]]$n # no. of data points
    if (k > (n/40)){
      fits_aic[i, "AIC"] <- estimate_aicc(fits[[i]]) # Estimating AICc
      #print("Estimating Corrected AIC")
    } else{
      fits_aic[i, "AIC"] <- fits[[i]]$aic  # Extracting AIC
      #print("Estimating AIC")
    }
  }
  # Order dataframe based on AIC
  fits_aic_ordered <- fits_aic[order(fits_aic$AIC),]
  
  return(fits_aic_ordered)
}

# Extract fitted distribution parameters
# Input:
#   fits: multi_fit object - a list of fitdist objects
# Output:
#   fits_params: Dataframe of fitted parameters for each distribution
extract_fitted_params <- function(fits){
  no_fits <- length(fits)
  fits_params <- data.frame("Distributions" = rep(NA,no_fits), 
                            "P1" = rep(NA,no_fits),
                            "P2" = rep(NA,no_fits))
  for (i in 1:no_fits){
    fits_params[i,"Distributions"] <- names(fits)[i]
    params_names <- names(fits[[i]]$estimate)
    params_val <- signif(as.vector(fits[[i]]$estimate),2)
    params_names_val <- paste(params_names, "=", params_val)
    
    fits_params[i,"P1"] <- params_names_val[1]
    fits_params[i,"P2"] <- params_names_val[2]
  }
  return(fits_params)
}

# Function: Estimate the difference in AIC
# Input:
#   fits: multi_fit object - a list of fitdist objects
# Output: 
#   daic: Vector of AIC differences
estimate_daic <- function(fits){
  ordered_aic <- estimate_aic(fits)
  daic <- ordered_aic$AIC - min(ordered_aic$AIC)
  return(daic)
}

# Function: Estimate the Akaike weights
# Input:
#  fits: multi_fit object - a list of fitdist objects
# Output: 
#   waic: Vector of AIC weights
estimate_waic <- function(fits){
  waic_num <- exp(-0.5*estimate_daic(fits))
  waic_den <- sum(waic_num)
  waic <- waic_num/waic_den
  return(waic)
}

# Summary of the goodness of fits
# Input:
#  fits: multi_fit object - a list of fitdist objects
# Output: 
#   gof_df: Dataframe of AIC metrics
gof_summary <- function(fits, params = TRUE){
  gof_df <- estimate_aic(fits)
  gof_df["dAIC"] <- estimate_daic(fits)
  gof_df["wAIC"] <- estimate_waic(fits)
  if (params == TRUE){
    params_df <- extract_fitted_params(fits)
    gof_df <- left_join(gof_df, params_df , by = "Distributions")
  }
  return(gof_df)
}


# Plotting Functions ---------

# Function to generate GOF plot and save to designated path for 
# Input:
#   path: path to save the plot
#   fits: multi_fit object - a list of fitdist objects
#   dists: vector of distributions to be plotted. Must be one of fits items names
#   legend_items: vector of names used for the legend items. Must be the same length as dists
#   x_label: Label use for the density and cdf plots.
# Output:
#   plot: pdf plot of GOF
plot_GOF <- function(filename, fits, dists, legend_items = dists, x_label = "data"){
  if(missing(filename)){
    filename <- "GOF_plot.pdf"
  }
  if(any(!dists %in% names(fits))){
    stop(paste("Dists items do not match fits items."))
  }
  pdf(filename)
  par(mfrow = c(2, 2))
  plot.legend <- legend_items
  denscomp(fits[dists], legendtext = plot.legend, xlab = x_label)
  qqcomp(fits[dists], legendtext = plot.legend)
  cdfcomp(fits[dists], legendtext = plot.legend, xlab = x_label)
  ppcomp(fits[dists], legendtext = plot.legend)
  dev.off()
  print(paste("Plot",filename,"was saved."))
}

# Function to generate GOF plot for fits on different datasets 
# Input:
#   path: path to save the plot
#   fits: List of fitdist objects
#   x_label: Label use for the density and cdf plots.
# Output:
#   plot: plot with a column for each fitdist object
plot_multivar_GOF <- 
  function(fits, titles, x_label = "data"){
    
    no_fits <- length(fits) # Number of fits depending on the number of datasets
    cols <- viridis(no_fits, begin = 0.2, end = 0.7) # Different colour schemes per dataset
    # Plot Settings
    par(mfcol = c(4,no_fits), mar=c(4, 2, 0, 0) + 0.1, oma = c(0, 2, 3, 0) + 0.1, mgp = c(2.1,0.7,0) )
    #par(mar=c(5.1,4.1,4.1,2.1))
    ## oma = A vector of the form c(bottom, left, top, right) giving the 
    # size of the outer margins in lines of text.
    ## mar = A numerical vector of the form c(bottom, left, top, right) 
    # which gives the number of lines of margin to be specified on the 
    # four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
    for (i in c(1:no_fits)){
      #print(i)
      denscomp(fits[[i]], addlegend = FALSE, main = "", xlab = x_label, fitcol = cols[i])
      qqcomp(fits[[i]], addlegend = FALSE, main = "", fitpch = 16, fitcol = cols[i], ylab = "", line01lty = 1)
      cdfcomp(fits[[i]], addlegend = FALSE, main = "", xlab = x_label, fitcol = cols[i], ylab = "", lines01 = TRUE)
      ppcomp(fits[[i]], addlegend = FALSE, main = "", fitpch = 16, fitcol = cols[i], ylab = "", line01lty = 1)
      #mtext(titles[i], side = 3, outer = TRUE, at = i/3 - 0.15, line = 1)
      mtext(titles[i], side = 3, outer = TRUE, at = i/(no_fits) - (1/(no_fits*2)), line = 1, col = cols[i])
    }
    mtext("Density", side = 2, outer = TRUE, at =  0.9, cex = 0.7)
    mtext("Empirical quantiles", side = 2, outer = TRUE, at  = 0.65, cex = 0.7)
    mtext("CDF", side = 2, outer = TRUE, at  = 0.40, cex = 0.8)
    mtext("Empirical probabilities", side = 2, outer = TRUE, at = 0.15, cex = 0.7)
    
  }
