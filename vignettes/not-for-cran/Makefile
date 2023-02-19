%.html: %.Rmd
	@echo "$< -> $@"
	@Rscript -e "rmarkdown::render('$<')"

estfun_DA.pdf: estfun_DA.tex
	@echo "$< -> $@"
	@pdflatex $<
	@pdflatex $<
