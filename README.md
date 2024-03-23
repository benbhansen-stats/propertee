# propertee: **P**rognostic **R**egression **O**ffsets with **P**ropagation of **ER**rors, for **T**reatment **E**ffect **E**stimation


<!-- badges: start -->
[![R-build-check](https://github.com/benbhansen-stats/propertee/workflows/R-build-check/badge.svg)](https://github.com/benbhansen-stats/propertee/actions)
<!-- badges: end -->

## Overview

Propertee enables flexible direct adjustment with design-informed standard errors
and optional prior covariance adjustment.

Random trials often utilize units of assignment and blocking in assigning
treatment status as a way to simplify implementation. This design information
must be utilized in future analyses. Using Propertee, a user can generate a
Design object which will keep track of the design structure.

    des <- rct_design(txt ~ unit_of_assignment(teacher) + block(school), data = teacherdata)

(Also supported are observational studies (`obs_design`) and regression
discontinuity designs (`rdd_design` which requires a `forcing()` variable as
well.)

In order to pass the design information into the model using the `weights=`
argument, functions `ett()` and `ate()` will be used to convert the Design into
a numeric vector with the Design object as an attribute.

    lm(y ~ txt, data = studentdata, weights = ate(des))

Note that the Design is created with teacher level data (`teacherdata`), but the
analysis is carried out at the student level (`studentdata`); the `ate()` (and
its alternative `ett()`) will expand the weights appropriately.

Optionally, we can also include a covariance adjustment model through the
`cov_adj()` function.

    covadjmod <- lm(y ~ x1 + x2 + ..., data = studentdata, subset = !txt)
    lm(y ~ txt, studentdata, weights = ett(des),
       offset = cov_adj(covadjmod, data = studentdata)
    )

## Contributing

You may use RStudio to develop for propertee, by opening the `propertee.Rproj` file.
We suggest you ensure all required dependencies are installed by running

    devtools::install_deps(dependencies = TRUE)

We prefer changes that include unit tests demonstrating the problem or showing
how the new feature should be added. The test suite uses the
[testthat](https://github.com/hadley/test_that) package to write and run tests.
(Please ensure you have the latest version of testthat (or at least v0.11.0), as
older versions stored the tests in a different directory, and may not test
properly.) See the `tests/testthat` directory for examples. You can run the test
suite via Build -> Test Package.

New features should include inline [Roxygen](http://roxygen.org/) documentation.
You can generate all `.Rd` documents from the `Roxygen` code using Build ->
Document, or using Make as describe below.

Finally, you can use Build -> Build and Reload or Build -> Clean and Rebuild to
load an updated version of `propertee` in your current RStudio session.
Alternatively, to install the developed version permanently, use Build -> Build
Binary Version, followed by

    install.packages("../propertee_VERSION.tgz", repo=NULL)

You can revert back to the current CRAN version by

    remove.packages("propertee")
    install.packages("propertee")

If you prefer not to use RStudio, you can develop using Make.

- `make test`: Run the full test suite.
- `make document`: Update all documentation from Roxygen inline comments.
- `make interactive`: Start up an interactive session with `propertee` loaded.
  (`make interactive-emacs` starts the session inside emacs.)
- `make check`: Run `R CMD check` on the package
- `make build`: Build a binary package.
- `make vignette`: Builds any vignettes in `vignettes/` directory
- `make clean`: Removes files built by `make vignette`, `make document` or `make
   check`. Should not be generally necessary, but can be useful for debugging.

When your change is ready, make a pull request on github.

### White space changes

To ease searches of the commit history:

- Commit white space changes only when they occur on lines with substantive
  changes.
- Avoid committing trailing white spaces.

In **RStudio**, there are options to enable automatically removing white space
as the end of lines and trailing whitespaces in the Settings, Code -> Saving.

In **emacs**, you can remove white spaces at ends of lines with `M-x
delete-trailing-whitespace`. To do this automatically whenever you save, add the
following to your init file:

    (add-hook 'before-save-hook (lambda ()
                                 (delete-trailing-whitespace)))

To remove trailing lines when saving, you can also add this:

    (setq delete-trailing-lines t)

### Internal functions

Any internal functions (for our use only) should be prefaced with a "`.`" (e.g.
`.myfunc <- function()`). Internal functions should be documented using roxygen
as described above, and given the `@keywords internal` tag to ensure they do not
get indexed. (Generally internally functions should not be `@export`'d but
exceptions may arise.)

During this period of development, after documenting an internal function, add
it to the "_pkgdown.yml" file in the appropriate category. Once propertee goes
public, we will remove those.

### Referring to functions

When documentation refers to another function (internal to the package or
otherwise), please include the trailing `()`, as that will help **pkgdown**
provide an appropriate link (see [https://pkgdown.r-lib.org/articles/linking.html](https://pkgdown.r-lib.org/articles/linking.html)).

E.g. `\code{lm()}` or `\code{cov_adj()}` or `\code{lme4::lmer()}`.

### Vignettes or simulations

Vignettes and simulations should go in the
[/vignettes/](https://github.com/benbhansen-stats/propertee/tree/main/vignettes)
folder. Anything that should not be built by R check (especially anything that
builds very slowly, or introduces dependencies not specified in
[DESCRIPTION](https://github.com/benbhansen-stats/propertee/blob/main/DESCRIPTION))
should go in the
[/vignettes/not-for-cran](https://github.com/benbhansen-stats/propertee/tree/main/vignettes/not-for-cran).
(The not-for-cran folder is in the
[.Rbuildignore](https://github.com/benbhansen-stats/propertee/blob/main/.Rbuildignore)
file.)

Note that ALL .Rmd files in /vignettes/ get built during building of the
reference site. To exclude a .Rmd file, it needs to start with "_". E.g.
`myslowvignette.Rmd` -> `_myslowvignette.Rmd`.
