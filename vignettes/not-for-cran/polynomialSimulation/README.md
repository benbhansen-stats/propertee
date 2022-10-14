### Polynomial Simulation Validation for `flexida`

The files in this folder comprise a workflow that validates the `flexida` standard
error estimation process. It reproduces simulations run in
[Sales and Hansen 2019](https://journals.sagepub.com/doi/full/10.3102/1076998619884904), calculates
standard errors using `flexida`, and checks whether these standard errors produce confidence intervals
with similar coverage rates to the confidence intervals in the manuscript.

#### Running Simulations
Simulation data is not stored in this folder, only the markdown document confirming the results.
To regenerate the simulation data, `cd` into this folder from a terminal window and run
`make polynomialSimulation.md NREPS=5000 clean`.
