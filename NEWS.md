# v0.1.1

## Breaking Changes
* We have reformulated the estimating equations used to derive standard errors. In estimation settings we accommodate, testing has not revealed any differences in standard error estimates between the previous and current estimating equations, but we do not assure this is the case for all possible situations.

# v0.0.1

* Compatible with R 4.2.3
* Introduces functionality for direct adjusted and design-informed standard errors accommodating covariance adjustment in the model-based setting
* Cluster-robust standard errors can only be estimated using the HC0 estimator
