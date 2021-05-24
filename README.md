# distrmultifit

The ```multi_fit.R``` builds upon the fitdistrplus R package to:

- Fit multiple distributions
- Compute the Corrected AIC, difference in AIC (information loss)  and Akaike weights for each fitted distribution.
- Rank the distributions based on the AIC metrics
- Generate Goodness-of-Fit (GOF) plots for the distributions specified.

A paper (accepted for publication) demonstrates the use of ```multi_fit.R``` for the application of identifying an appropriate distribution for a set of wall U-value measurements from English homes. However, the method can be applied to any case where multiple distributions are to be fitted and compared.  Use ```demo_beyond_normal.R``` to reproduce the published work or adjust it to fit your own analysis. The method is adapted from Delignette-Muller and Durang (2015).

If you are using this script for your own analysis, please acknowledge it by referencing the published paper using the following reference:

Petrou, G., Mavrogianni, A., Symonds, P. and Davies, M., (2021). 'Beyond Normal: Guidelines on How to Identify Suitable Model Input Distributions for Building Performance Analysis', *Building Simulation 2021*. Bruges, Belgium, 1-3rd Sept 2021.

## Dependencies

- fitdistrplus 
- ggplot2
- dplyr
- viridis

## References

Delignette-Muller, M.L., Dutang, C., 2015. fitdistrplus: An R Package for Fitting Distributions. Journal of Statistical Software 64, 1–34. https://doi.org/10.18637/jss.v064.i04

