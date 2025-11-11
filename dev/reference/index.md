# Package index

## StudySpecification

### Creation of `StudySpecification` Objects

- [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`rct_specification()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`rd_specification()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`obs_specification()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`obsstudy_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  [`obsstudy_specification()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  :

  Generates a `StudySpecification` object with the given specifications.

- [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  [`unitid()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  [`cluster()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  [`uoa()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  [`block()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  :

  Special terms in `StudySpecification` creation formula

- [`as_rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/specificationconversion.md)
  [`as_obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/specificationconversion.md)
  [`as_rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/specificationconversion.md)
  :

  Convert `StudySpecification` between types

### Accessors/Replacers for `StudySpecification` Objects

- [`treatment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `treatment<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`units_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `units_of_assignment<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`clusters()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `clusters<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`unitids()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `unitids<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`blocks()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `blocks<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`has_blocks()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`forcings()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  [`` `forcings<-`() ``](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
  :

  Accessors and Replacers for `StudySpecification` objects

- [`has_binary_treatment()`](https://benbhansen-stats.github.io/propertee/dev/reference/has_binary_treatment.md)
  :

  Check whether treatment stored in a `StudySpecification` object is
  binary

### Summary and print methods

- [`show(`*`<StudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/show-StudySpecification-method.md)
  :

  Show a `StudySpecification`

- [`summary(`*`<StudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_summary.md)
  [`print(`*`<summary.StudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_summary.md)
  :

  Summarizing `StudySpecification` objects

### Structure of `StudySpecification` objects

- [`get_structure()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_structure.md)
  [`show(`*`<StudySpecificationStructure>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_structure.md)
  :

  `StudySpecification` Structure Information

- [`specification_table()`](https://benbhansen-stats.github.io/propertee/dev/reference/specification_table.md)
  [`stable()`](https://benbhansen-stats.github.io/propertee/dev/reference/specification_table.md)
  :

  Table of elements from a `StudySpecification`

- [`var_table()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_var_names.md)
  [`var_names()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_var_names.md)
  :

  Extract Variable Names from `StudySpecification`

### Utility Functions for `StudySpecification` Objects

- [`identify_small_blocks()`](https://benbhansen-stats.github.io/propertee/dev/reference/identify_small_blocks.md)
  : Identify fine strata

- [`specification_data_concordance()`](https://benbhansen-stats.github.io/propertee/dev/reference/specification_data_concordance.md)
  : Check for variable agreement within units of assignment

- [`identical_StudySpecifications()`](https://benbhansen-stats.github.io/propertee/dev/reference/identical_StudySpecifications.md)
  :

  Test equality of two `StudySpecification` objects

## Weights

Functions to create or interact with Weights

- [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`att()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`etc()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`atc()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`ato()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`olw()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`owt()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  [`pwt()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  : Generate Direct Adjusted Weights for Treatment Effect Estimation

### Working with `WeightedStudySpecification` objects

- [`weights(`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/weights-WeightedStudySpecification-method.md)
  :

  Extract Weights from `WeightedStudySpecification`

- [`` `+`( ``*`<WeightedStudySpecification>`*`,`*`<numeric>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `+`( ``*`<numeric>`*`,`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `-`( ``*`<WeightedStudySpecification>`*`,`*`<numeric>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `-`( ``*`<numeric>`*`,`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `*`( ``*`<WeightedStudySpecification>`*`,`*`<numeric>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `*`( ``*`<numeric>`*`,`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `/`( ``*`<WeightedStudySpecification>`*`,`*`<numeric>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  [`` `/`( ``*`<numeric>`*`,`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecificationOps.md)
  :

  `WeightedStudySpecification` Operations

- [`c(`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/c-WeightedStudySpecification-method.md)
  : Concatenate weights

- [`subset(`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecification.subset.md)
  [`` `[`( ``*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightedStudySpecification.subset.md)
  :

  `WeightedStudySpecification` subsetting

- [`show(`*`<WeightedStudySpecification>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/show-WeightedStudySpecification-method.md)
  :

  Show a `WeightedStudySpecification`

## Covariance Adjustment & SandwichLayer

- [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  :

  Covariance adjustment of `teeMod` model estimates

- [`as.SandwichLayer()`](https://benbhansen-stats.github.io/propertee/dev/reference/as.SandwichLayer.md)
  :

  Convert a `PreSandwichLayer` to a `SandwichLayer` with a
  `StudySpecification` object

- [`subset(`*`<PreSandwichLayer>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/PreSandwichLayer.subset.md)
  [`` `[`( ``*`<PreSandwichLayer>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/PreSandwichLayer.subset.md)
  :

  `PreSandwichLayer` and `SandwichLayer` subsetting

- [`show(`*`<PreSandwichLayer>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/PreSandwichLayer.show.md)
  :

  Show a `PreSandwichLayer` or `SandwichLayer`

- [`bread(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/bread.teeMod.md)
  :

  Extract bread matrix from a `teeMod` model fit

- [`estfun(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/estfun.teeMod.md)
  :

  Extract empirical estimating equations from a `teeMod` model fit

- [`estfun(`*`<lmrob>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/lmrob_methods.md)
  [`bread(`*`<lmrob>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/lmrob_methods.md)
  :

  Generate matrix of estimating equations for `lmrob()` fit

- [`estfun(`*`<glmrob>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/glmrob_methods.md)
  [`bread(`*`<glmrob>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/glmrob_methods.md)
  :

  Extract empirical estimating equations from a `glmbrob` model fit

## Estimating Model

Functions to carry out the treatment estimation accounting for the
StudySpecification

- [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  : Linear Model for Intention To Treat

- [`assigned()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  [`adopters()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  [`a.()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  [`z.()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  : Obtain Treatment from StudySpecification

- [`as.lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/as_lmitt.md)
  [`as.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/as_lmitt.md)
  :

  Convert `lm` object into `teeMod`

- [`confint(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/confint.teeMod.md)
  :

  Confidence intervals with standard errors provided by
  [`vcov.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/vcov.teeMod.md)

- [`show(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/show-teeMod-method.md)
  :

  Show a `teeMod`

- [`summary(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  [`print(`*`<summary.teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  :

  Summarizing `teeMod` objects

- [`vcov(`*`<teeMod>`*`)`](https://benbhansen-stats.github.io/propertee/dev/reference/vcov.teeMod.md)
  :

  Compute variance-covariance matrix for fitted `teeMod` model

- [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  [`.vcov_DB0()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  [`.vcov_DB()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  :

  Variance/Covariance for `teeMod` objects

## Data

- [`STARplus`](https://benbhansen-stats.github.io/propertee/dev/reference/STARplus.md)
  : STAR participants plus nonexperimental controls
- [`lsoSynth`](https://benbhansen-stats.github.io/propertee/dev/reference/lsoSynth.md)
  : Synthethic Regression Discontinuity Data
- [`schooldata`](https://benbhansen-stats.github.io/propertee/dev/reference/studentdata.md)
  [`studentdata`](https://benbhansen-stats.github.io/propertee/dev/reference/studentdata.md)
  : Student data
- [`simdata`](https://benbhansen-stats.github.io/propertee/dev/reference/simdata.md)
  : Simulated data
- [`michigan_school_pairs`](https://benbhansen-stats.github.io/propertee/dev/reference/michigan_school_pairs.md)
  : Intervention data from a pair-matched study of schools in Michigan
- [`GV_data`](https://benbhansen-stats.github.io/propertee/dev/reference/GV_data.md)
  : Cluster-randomized experiment data on voter turnout in cable system
  markets
