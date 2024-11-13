polynomialSimulation.md: totalPolySim.Rdata polynomialSimulation.Rmd
	Rscript -e "rmarkdown::render('polynomialSimulation.Rmd')"

totalPolySim.Rdata totalPolySim-runtime.txt: runSimulations.R
	Rscript runSimulations.R $(NREPS)

clean:
	rm totalPolySim.Rdata totalPolySim-runtime.txt
