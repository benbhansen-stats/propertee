# propertee: **P**rognostic **R**egression **O**ffsets with **P**ropagation of **ER**rors, for **T**reatment **E**ffect **E**stimation

**propertee** enables flexible direct adjustment with
specification-informed standard errors and optional prior covariance
adjustment.

Random trials often utilize units of assignment and blocking in
assigning treatment status as a way to simplify implementation. Standard
errors and other inferential calculations should respect such aspects of
the study’s design. Using **propertee**, a user can generate a
“StudySpecification” object which will keep track of the design
specification.

``` R
spec <- rct_spec(txt ~ unit_of_assignment(teacher) + block(school), data = teacherdata)
```

(Also supported are observational studies
([`obs_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md))
and regression discontinuity specifications (`rdd_spec()` which requires
a
[`forcing()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
variable as well.)

In order to pass the design specification into the model using the
`weights=` argument, functions
[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
and
[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
will be used to convert the StudySpecification into a numeric vector
with the StudySpecification object as an attribute.

``` R
lm(y ~ txt, data = studentdata, weights = ate(spec))
```

Note that this StudySpecification is created with teacher level data
(`teacherdata`), but the analysis is carried out at the student level
(`studentdata`); the
[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
(and its alternatives
[`att()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md),
[`atc()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
and
[`ato()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md))
will expand the weights appropriately.

Optionally, we can also include a covariance adjustment model through
the
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
function.

``` R
covadjmod <- lm(y ~ x1 + x2 + ..., data = studentdata, subset = !txt)
lm(y ~ txt, studentdata, weights = ett(spec),
   offset = cov_adj(covadjmod, data = studentdata)
)
```

The effect estimation model is expected to be an lm object, but the
covariance model can be of type `glm` or `lmrob`, among others: with a
binary outcome, for example, it might be a logistic regression.

Having fit covariance and effect estimation models separately, you may
need standard errors that propagate variation from the former to the
latter. propertee’s “lmitt” objects are lm’s enhanced with information
from the StudySpecification and covariance adjustment model that
standard error calculations may need. Fitted lm’s of the accommodated
type may be coerced to lmitt using
[`as.lmitt()`](https://benbhansen-stats.github.io/propertee/reference/as_lmitt.md);
or you can skip this step by using
[`propertee::lmitt()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md)
to fit the lm in the first place:

``` R
lmitt(y ~ 1, studentdata, specificaton = spec, weights = "ett",
      offset = cov_adj(covadjmod, data = studentdata)
  )
```

For more, see the package’s vignettes.

## Contributing

You may use RStudio to develop for **propertee**, by opening the
`propertee.Rproj` file. We suggest you ensure all required dependencies
are installed by running

``` R
devtools::install_deps(dependencies = TRUE)
```

We prefer changes that include unit tests demonstrating the problem or
showing how the new feature should be added. The test suite uses the
[testthat](https://github.com/hadley/test_that) package to write and run
tests. (Please ensure you have the latest version of testthat (or at
least v0.11.0), as older versions stored the tests in a different
directory, and may not test properly.) See the `tests/testthat`
directory for examples. You can run the test suite via Build -\> Test
Package.

New features should include inline [Roxygen](http://roxygen.org/)
documentation. You can generate all `.Rd` documents from the `Roxygen`
code using Build -\> Document, or using Make as describe below.

Finally, you can use Build -\> Build and Reload or Build -\> Clean and
Rebuild to load an updated version of **propertee** in your current
RStudio session. Alternatively, to install the developed version
permanently, use Build -\> Build Binary Version, followed by

``` R
install.packages("../propertee_VERSION.tgz", repo=NULL)
```

You can revert back to the current CRAN version by

``` R
remove.packages("propertee")
install.packages("propertee")
```

If you prefer not to use RStudio, you can develop using Make.

- `make test`: Run the full test suite.
- `make document`: Update all documentation from Roxygen inline
  comments.
- `make interactive`: Start up an interactive session with **propertee**
  loaded. (`make interactive-emacs` starts the session inside emacs.)
- `make check`: Run `R CMD check` on the package
- `make build`: Build a binary package.
- `make vignette`: Builds any vignettes in `vignettes/` directory
- `make clean`: Removes files built by `make vignette`, `make document`
  or `make check`. Should not be generally necessary, but can be useful
  for debugging.

When your change is ready, make a pull request on github.

### White space changes

To ease searches of the commit history:

- Commit white space changes only when they occur on lines with
  substantive changes.
- Avoid committing trailing white spaces.

In **RStudio**, there are options to enable automatically removing white
space as the end of lines and trailing whitespaces in the Settings, Code
-\> Saving.

In **emacs**, you can remove white spaces at ends of lines with
`M-x delete-trailing-whitespace`. To do this automatically whenever you
save, add the following to your init file:

``` R
(add-hook 'before-save-hook (lambda ()
                             (delete-trailing-whitespace)))
```

To remove trailing lines when saving, you can also add this:

``` R
(setq delete-trailing-lines t)
```

### Referring to functions

When documentation refers to another function (internal to the package
or otherwise), please include the trailing `()`, as that will help
**pkgdown** provide an appropriate link (see
<https://pkgdown.r-lib.org/articles/linking.html>).

E.g. `\code{lm()}` or `\code{cov_adj()}` or `\code{lme4::lmer()}`.

### Vignettes or simulations

Vignettes and simulations should go in the
[/vignettes/](https://github.com/benbhansen-stats/propertee/tree/main/vignettes)
folder. Anything that should not be built by R check (especially
anything that builds very slowly, or introduces dependencies not
specified in
[DESCRIPTION](https://github.com/benbhansen-stats/propertee/blob/main/DESCRIPTION))
should go in the
[/vignettes/not-for-cran](https://github.com/benbhansen-stats/propertee/tree/main/vignettes/not-for-cran).
(The not-for-cran folder is in the
[.Rbuildignore](https://github.com/benbhansen-stats/propertee/blob/main/.Rbuildignore)
file.)

Note that ALL .Rmd files in /vignettes/ get built during building of the
reference site. To exclude a .Rmd file, it needs to start with “\_“.
E.g. `myslowvignette.Rmd` -\> `_myslowvignette.Rmd`.
