### Polynomial Simulation Validation for `propertee`

The files in this folder comprise a workflow that validates the `propertee` standard
error estimation process. It reproduces simulations run in
[Sales and Hansen 2019](https://journals.sagepub.com/doi/full/10.3102/1076998619884904), calculates
standard errors using `propertee`, and checks whether these standard errors produce confidence intervals
with similar coverage rates to the confidence intervals in the manuscript. The results can be viewed in
[`polynomialSimulation.md`](polynomialSimulationl.md).

#### Running Simulations
Simulation data is not stored in this folder, only the markdown document with the results.
To regenerate the simulation data, `cd` into this folder from a terminal window and run
`make polynomialSimulation.md NREPS=5000 clean` (the `clean` addendum removes the artifacts generated
during the simulations, such as the simulation data. When debugging either the `runSimulations.R` or
`polynomialSimulation.Rmd` files, it may be helpful to remove this addendum and to run fewer simulations,
i.e. run the make command with `NREPS=10`.) Running the make command with `NREPS=0` does not re-run
simulations and only regenerates the markdown document if simulation data is available.
