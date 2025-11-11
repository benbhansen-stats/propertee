# OLS effect and error estimation with an auxiliary sample and a separate, not-necessarily-linear covariance model

## Overview

The `propertee` package offers tools for enhancing evaluations of
treatments, policies, and interventions that respect the statistical
properties endowed by the study specification. One such offering is a
routine for covariance adjustment that allows researchers to model
exogenous variation in outcomes of interest using data from their study
as well as available auxiliary data. `propertee` accommodates linear,
generalized linear, and robust regression models, providing users
flexibility in functional form, fitting procedure, and fitting sample.
This vignette serves as a step-by-step walkthrough of how users can use
a prior covariance adjustment model fit to inform estimates of
intervention effects–as well as associated standard errors–with the
software in `propertee`.

Pane et al. (2014) mounted a large-scale cluster randomized trial in
seven states to study the effectiveness of Cognitive Tutor, an
online/in-person blended algebra learning program. Their report assessed
program effectiveness in terms of scores on a test administered as part
of the study. However, scores on prior and subsequent state achievement
tests are available for study schools as well as others in their states
and districts, and these could be used for a complementary, perhaps more
precise, assessment of the program’s effects. Here we demonstrate this
idea using state and district data from Michigan, where the Cognitive
Tutor study had a significant footprint and where rich school
achievement data are available for download from state websites.

Pane and coauthors kindly shared with us the names, pre-randomization
pairings and treatment assignments of their study’s 14 Michigan schools,
but were not at liberty to make this information public. To create a
functionally similar case study while maintaining the participating
schools’s anonymity, we optimally pair-matched them to schools from a
large nearby county in Michigan, Oakland. The pseudo-RCT considered in
this vignette replaces each Michigan Cognitive Tutor study school with
the Oakland County school it was paired to, otherwise inheriting from
the actual RCT salient specification characteristics, such as the
composition of school pairs and triples within which randomization was
conducted.

Valid analysis of an RCT calls for careful attention to such
characteristics, both for selecting a compatible estimator and for
correctly implementing it. By introducing a dedicated S4 structure for
such characteristics, the `StudySpecification`, along with
StudySpecification-aware functions for such tasks as inverse probability
weighting and effect estimation with optional stratum fixed effects,
`propertee` helps the analyst stay on top of implementation details. Its
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md),
a specialized [`predict()`](https://rdrr.io/r/stats/predict.html),
decouples effect estimation from covariance adjustment while continuing
to track what’s necessary for valid standard error estimation,
significantly broadening the range of estimators that are compatible
with a given `StudySpecification`.

## Data

To run this vignette, first download the necessary data. We will use
school-level averages of student performance on the Michigan Merit
Examination (MME) in 2014 to measure intervention effects, and we will
use school-level averages of scores on the 2012 and 2013 tests as
covariates in our covariance adjustment model. These scores can be
downloaded as a zipped file from the [Michigan Department of Education
website](https://mischooldata.org/historical-assessment-data-files/).
After unzipping that file, convert the resulting .xls file to a .csv to
facilitate the use of base R commands for loading it into an R session.

We also use school-level characteristics from the [Common Core of
Data](https://nces.ed.gov/ccd/files.asp#Fiscal:2,LevelId:7,SchoolYearId:28,Page:1)
in the covariance adjustment model. Click on the link provided here, and
download the “Flat File” for the 2013-2014 Public Elementary/Secondary
School Universe Survey under “Data File”. We can use base R commands to
load in the unzipped .txt file.

The last necessary file can be be loaded from the `propertee` package by
calling `data(michigan_school_pairs)`. This dataframe tracks which
schools were paired together in the study and which schools were
assigned to intervention and control.

## Package Installation

This vignette requires the installation of three package in addition to
`propertee`: `httr` and `readxl`, which we’ll use to import the Michigan
schools data; and `robustbase`, providing functionality for
outlier-robust regression. The vignette uses `robustbase` to demonstrate
how `propertee` can handle alternatives to ordinary least squares for
covariance adjustment modeling.

## Walkthrough

After loading the installed packages and reading in the downloaded data
files, we clean the MME scores and school characteristics datasets. The
scores data has rows corresponding to state-, intermediate school
district (ISD)-, district-, and campus-wide averages. In addition to
averages taken over all students in these subpopulations, some rows
correspond to averages taken within substrata formed by gender,
race/ethnicity, learning ability, or economic background. In the
provided cleaning script (`get_and_clean_external_data.R`), we create
two cleaned datasets, one for an analysis of the marginal effect of the
intervention and one for an analysis of the heterogeneity of the
intervention effect. Both datasets keep only rows corresponding to
schools where MME scores were reported campus-wide or for the particular
substratum in each of 2012, 2013, and 2014.

The school characteristics data spans the universe of public schools in
the United States, so to clean it for this vignette, we first limit it
to schools relevant to the study. The MME is taken almost exclusively by
11th graders and, as the name suggests, only taken by students in
Michigan, so we first subset the data to schools in Michigan serving
11th graders. Then, we perform feature generation, creating derived
covariates such as demographic breakdowns by gender, race/ethnicity, and
free- or reduced-price lunch eligibility at the school level and in the
11th grade specifically. (The provided cleaning script performs these
steps also.)

``` r
if (!require("robustbase")) library(robustbase)
if (!require("readxl")) library(readxl)
if (!require("httr")) library(httr)
if (!require("propertee")) library(propertee)

extdataURLs <- list(
CCD="https://nces.ed.gov/ccd/data/zip/sc132a_txt.zip",
MME="https://www.michigan.gov/cepi/-/media/Project/Websites/cepi/MiSchoolData/historical/Historical_Assessments/2011-2014MME.zip"
)
data(michigan_school_pairs)
michigan_school_pairs <-
     subset(michigan_school_pairs, !is.na(blk), c(schoolid,blk,z))
source("get_and_clean_external_data.R")
```

#### Creating the `StudySpecification` Object

The first step in estimating intervention effects using `propertee` is
to create a `StudySpecification` object. This will store the information
from `michigan_school_pairs` in a way that will allow for quick
calculation of inverse probability of assignment weights that attend to
the pair-matched structure of the study.

In studies where units of observation (eg. individual students) within
units of assignment (eg. schools or districts) are used to estimate
intervention effects, the data structure should reflect this nesting.
For example, if we had student-level scores data as the observed data
within schools or districts, it is important to specify how these units
of assignment are allocated to different treatment or control
conditions. In such studies, the `StudySpecification` object would also
facilitate the definition of a vector of assignment indicators at the
unit of observation level we could use for estimating the intervention
effect.

Should more than one variable be needed to identify the unit of
assignment, block, or forcing, they can be included. For example,
perhaps `schoolidk` may be unique within district, but potentially not
unique across districts. Then we’d use something like
`block(districtid, schoolidk)` in the `_spec` function.

``` r
spec <- rct_spec(z ~ unitid(schoolid) + block(blk), michigan_school_pairs)
```

In this
[`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
call, the lefthand side indicates the assignment variable, and the
righthand side indicates the unit of assignment and, if applicable,
variable that identify matched sets or strata. There are also
[`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
and
[`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
constructors, for regression discontinuity and for observational
studies/quasiexperimental specifications, respectively.

The `StudySpecification` object’s `structure` slot lists units of
assignment, their allocations to conditions and, if applicable, blocks
within which they were allocated.

``` r
spec@structure
```

    ##    z   schoolid blk
    ## 1  1 6307005976   A
    ## 2  0 6315004226   A
    ## 3  1 6316006171   B
    ## 4  0 6329004340   B
    ## 5  1 6314002317   C
    ## 6  0 6324009415   C
    ## 7  1 6318000385   D
    ## 8  0 6320001204   D
    ## 9  0 6301004608   E
    ## 10 0 6305000291   E
    ## 11 1 6328002123   E
    ## 12 1 6326003242   F
    ## 13 0 6326005819   F
    ## 14 1 6327000710   F

#### Fitting the Covariance Adjustment Model

We now fit the covariance adjustment model. The `propertee` package will
generate predictions from this regression to explain residual variation
of the outcomes in the study. Often, this produces more accurate and
precise effect estimates. As mentioned earlier, `propertee` accommodates
a host of fitting procedures for estimating this model, and the
regression may leverage data from available auxiliary sources. We
demonstrate this flexibility by fitting two covariance adjustment models
for each analysis, one with least squares and one with robust
regression. We fit these models to a sample including all schools in the
study and all schools in Oakland County whose outcomes and covariates
are measured in the data we’ve downloaded. The exact specification for
this school-level model is provided in Equation 1.

$$\begin{aligned}
\text{Average\_Score\_14} & {= \beta_{0} + \beta_{1}\text{Total\_Enrollment\_14} + \beta_{2}\text{Title\_I\_Status\_14}} \\
 & {\quad + \beta_{3}\text{Magnet\_Status\_14} + \beta_{4}\text{Charter\_Status\_14} + \beta_{5}\text{Education\_Type\_14}} \\
 & {\quad + \beta_{6}\text{Perc\_Female\_14} + \beta_{7}\text{Perc\_Native\_14} + \beta_{8}\text{Perc\_Asian\_14}} \\
 & {\quad + \beta_{9}\text{Perc\_Hispanic\_14} + \beta_{10}\text{Perc\_Black\_14} + \beta_{11}\text{Perc\_White\_14}} \\
 & {\quad + \beta_{12}\text{Perc\_Pacific\_Islander\_14} + \beta_{13}\text{Perc\_Econ\_Disadvantaged\_14}} \\
 & {\quad + \beta_{14}\text{Perc\_Female\_Grade11\_14} + \beta_{15}\text{Perc\_Native\_Grade11\_14}} \\
 & {\quad + \beta_{16}\text{Perc\_Asian\_Grade11\_14} + \beta_{17}\text{Perc\_Hispanic\_Grade11\_14}} \\
 & {\quad + \beta_{18}\text{Perc\_Black\_Grade11\_14} + \beta_{19}\text{Perc\_White\_Grade11\_14}} \\
 & {\quad + \beta_{20}\text{Perc\_Pacific\_Islander\_Grade11\_14} + \beta_{21}\text{Average\_Score\_13}} \\
 & {\quad + \beta_{22}\text{Average\_Score\_12}}
\end{aligned}$$

``` r
coname <- "OAKLAND COUNTY"
RESPONSE_COL <- "Average.Scale.Score.2014"
MODELING_COLS <- c(
  "TOTAL_ENROLLMENT", setdiff(CCD_CAT_COLS, "TYPE"),
  setdiff(colnames(analysis1data)[grepl("_PERC$", colnames(analysis1data))],
          c("MALE_PERC", "TR_PERC", "MALE_G11_PERC", "TR_G11_PERC")),
  paste0("Average.Scale.Score.", c(2013, 2012))
)
```

``` r
not_missing_resp <- !is.na(analysis1data[[RESPONSE_COL]])
not_missing_covs <- rowSums(is.na(analysis1data[, MODELING_COLS])) == 0
county_ix <- analysis1data$CONAME == coname
county_camod_dat <- analysis1data[not_missing_resp & not_missing_covs,]

camod_form <- as.formula(
  paste0(RESPONSE_COL, "~", paste(MODELING_COLS, collapse = "+")))
lm_county_camod <- lm(camod_form, county_camod_dat,
                      weights = county_camod_dat$Total.Tested.2014)
```

``` r
set.seed(650)
rob_county_camod <- robustbase::lmrob(
  camod_form, county_camod_dat, weights = county_camod_dat$Total.Tested.2014,
  control = robustbase::lmrob.control(max.it = 500L))
```

#### Estimating Marginal Intervention Effects

With the `StudySpecification` object created and the covariance
adjustment model fit, we’re prepared to evaluate the intervention.
`propertee` supports the calculations of inverse probability of
assignment weights, which can be used to estimate either the average
intervention effect (ATE) or the average effect for the treated (ETT).
These weights can be combined with additional unit weights to reflect
varying sizes of units in the sample. In this analysis and the one that
follows, we estimate the student-level ATE by calculating inverse
probability of assignment weights for each school using the
[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
function, then multiplying those weights by the number of students at
the corresponding school who took the test. We incorporate the
prognostic model using the
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
function, which generates model predictions for each school and
identifies overlap between the prognostic sample and the study sample.

``` r
study1data <-
     merge(michigan_school_pairs, analysis1data, by = "schoolid", all.x = TRUE)
ip_wts <- propertee::ate(spec, data = study1data) * study1data$Total.Tested.2014
lm_ca <- propertee::cov_adj(lm_county_camod, newdata = study1data, specification = spec)
```

Next, we estimate the intervention effect by using the
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
function. The only major difference between
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
and the base [`lm()`](https://rdrr.io/r/stats/lm.html) function in the
requirement to specify a `StudySpecification` object. In this
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
function, we pass the inverse probability of assignment weights to the
`weights` argument and the adjusted predictions to the `offset`
argument.

``` r
main_effect_fmla <- as.formula(paste0(RESPONSE_COL, "~1"))
lm_ca_effect <- propertee::lmitt(
  main_effect_fmla, specification = spec, data = study1data, weights = ip_wts,
  offset = lm_ca
)
```

The summary of a fitted
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
model, which is called a `teeMod`, shows the estimated intervention
effect and the estimated standard error that has propagated uncertainty
from the covariance adjustment regression.

``` r
summary(lm_ca_effect, vcov.type = "HC0")
```

    ## 
    ## Call:
    ## lmitt.formula(main_effect_fmla, specification = spec, data = study1data, weights = ip_wts, offset = lm_ca)
    ## 
    ##  Treatment Effects :
    ##                                       Estimate Std. Error t value Pr(>|t|)    
    ## z.                                      0.4732     1.0903   0.434    0.664    
    ## Average.Scale.Score.2014:(Intercept) 1116.2153     3.9336 283.764   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Std. Error calculated via type "HC0"

Since no actual intervention was implemented in the pseudo-RCT, the
result of no effect is as expected.

To explore different variance estimation techniques, we can specify a
different `vcov.type` in the
[`summary()`](https://rdrr.io/r/base/summary.html) function.

``` r
summary(lm_ca_effect, vcov.type = "HC1")
```

    ## 
    ## Call:
    ## lmitt.formula(main_effect_fmla, specification = spec, data = study1data, weights = ip_wts, offset = lm_ca)
    ## 
    ##  Treatment Effects :
    ##                                       Estimate Std. Error t value Pr(>|t|)    
    ## z.                                      0.4732     1.1143   0.425    0.671    
    ## Average.Scale.Score.2014:(Intercept) 1116.2153     3.9336 283.764   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Std. Error calculated via type "HC1"

If parts of the auxiliary sample (here, Oakland County other than the 14
study schools) follow a different pattern of association between
covariates and response, the covariance adjustment might wind up doing
more harm than good. To increase robustness to such contamination, we
can use robust linear regression for generating predictions.

``` r
rob_ca <- propertee::cov_adj(rob_county_camod, newdata = study1data, specification = spec)
rob_ca_effect <- propertee::lmitt(
  main_effect_fmla, specification = spec, data = study1data, weights = ip_wts,
  offset = rob_ca
)
summary(rob_ca_effect, vcov.type = "HC1")
```

    ## 
    ## Call:
    ## lmitt.formula(main_effect_fmla, specification = spec, data = study1data, weights = ip_wts, offset = rob_ca)
    ## 
    ##  Treatment Effects :
    ##                                       Estimate Std. Error t value Pr(>|t|)    
    ## z.                                      0.4446     1.1304   0.393    0.694    
    ## Average.Scale.Score.2014:(Intercept) 1116.2153     3.9336 283.764   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Std. Error calculated via type "HC1"

#### Estimating Heterogeneous Intervention Effects

We can estimate average intervention effects conditional on different
race/ethinicity groups by using the second cleaned dataset. `propertee`
allows us to interpret the heterogenous effect estimates as the average
effect of the intervention on the average MME score for students within
each race/ethnicity group. The process for estimating these
heterogeneous effects follows a similar user experience to the
estimation of the marginal effect, but with two notable exceptions.

##### Exception 1: Formula Specification

The first exception is general to heterogeneous effect estimation with
the `propertee` package. When estimating these effects, users need to
modify the formula passed to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md).
Instead of specifying a constant (i.e. 1) on the right-hand side, users
should provide a formula that specifies the subgroup variable on the
right-hand side.

The first part of the code prepares the data by filtering out any rows
with missing values in the response variable or the covariates. We also
isolate data for the Oakland County.

``` r
not_missing_resp <- !is.na(analysis2data[[RESPONSE_COL]])
not_missing_covs <- rowSums(is.na(analysis2data[, MODELING_COLS])) == 0
county_ix <- analysis2data$CONAME == coname
county_mod_camod_dat <- analysis2data[not_missing_resp & not_missing_covs,]
```

Next, we update the model formula to include the subgroup variable
`DemographicGroup`, which will be used to estimate heterogeneous
effects. We then use this updated formula to fit a weighted linear
regression model.

``` r
mod_camod_form <- update(camod_form, . ~ . + factor(DemographicGroup))
lm_county_mod_camod <- lm(mod_camod_form, county_mod_camod_dat,
                          weights = county_mod_camod_dat$Total.Tested.2014)
```

Then, we prepare the study dataset by merging the cleaned dataset with
`michigan_school_pairs` to ensure we have school-level data that allows
for the estimation of treatment effects based on school and demographic
group.

``` r
study2data <- merge(michigan_school_pairs, analysis2data, by = "schoolid", all.x = TRUE)
study2data <- study2data[study2data$DemographicGroup %in%
                           c("White", "Black or African American"),]
```

Now, we compute the inverse probability weights using the
[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
function, which estimates the average effect treatment. Then, we adjust
for covariates using
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md).

``` r
ip_wts <- propertee::ate(spec, data = study2data) * study2data$Total.Tested.2014
lm_mod_ca <- propertee::cov_adj(lm_county_mod_camod, newdata = study2data,
                                specification = spec, by = "uniqueid")
```

Finally, we specify a formula for estimating heterogeneous effects and
fit the model using the
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
function.

``` r
mod_effect_fmla <- as.formula(paste0(RESPONSE_COL, "~ DemographicGroup"))
lm_ca_mod_effect <- propertee::lmitt(mod_effect_fmla, specification = spec,
                                     data = study2data, weights = ip_wts,
                                     offset = lm_mod_ca)
summary(lm_ca_mod_effect, vcov.type = "CR1", cluster = "schoolid")
```

    ## 
    ## Call:
    ## lmitt.formula(mod_effect_fmla, specification = spec, data = study2data, weights = ip_wts, offset = lm_mod_ca)
    ## 
    ##  Treatment Effects :
    ##                                                                    Estimate
    ## `z._DemographicGroupBlack or African American`                        0.326
    ## z._DemographicGroupWhite                                              1.337
    ## Average.Scale.Score.2014:DemographicGroupBlack or African American 1090.477
    ## Average.Scale.Score.2014:DemographicGroupWhite                     1115.743
    ##                                                                    Std. Error
    ## `z._DemographicGroupBlack or African American`                          1.990
    ## z._DemographicGroupWhite                                                1.290
    ## Average.Scale.Score.2014:DemographicGroupBlack or African American      2.347
    ## Average.Scale.Score.2014:DemographicGroupWhite                          2.716
    ##                                                                    t value
    ## `z._DemographicGroupBlack or African American`                       0.164
    ## z._DemographicGroupWhite                                             1.036
    ## Average.Scale.Score.2014:DemographicGroupBlack or African American 464.646
    ## Average.Scale.Score.2014:DemographicGroupWhite                     410.767
    ##                                                                    Pr(>|t|)    
    ## `z._DemographicGroupBlack or African American`                         0.87    
    ## z._DemographicGroupWhite                                               0.30    
    ## Average.Scale.Score.2014:DemographicGroupBlack or African American   <2e-16 ***
    ## Average.Scale.Score.2014:DemographicGroupWhite                       <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Std. Error calculated via type "CR1"

##### Second Exception: Overlapping Rows

The second exception arises when units of assignment (in this case,
schools) contribute multiple observations to the heterogeneous effect
estimation. This can happen because schools may have students in
multiple race/ethnicity groups, leading to overlap in the data.The
`StudySpecification` object does not provide enough information to
uniquely identify these rows, which causes an issue for standard error
calculations that must determine the exact overlap between the
covariance adjustment and effect estimation samples. To alleviate this
issue, both dataframes must have a column that uniquely identifies each
row. If the two dataframes have overlapping rows, these unique
identifiers should match up.

#### Conclusion

This vignette has demonstrated how to use the `propertee` package to
enhance the evaluation of treatment effects by incorporating covariance
adjustment models. The following key concepts and commands were covered:

- Create a `StudySpecification` object which encodes the design,
  including the unit of assignment, treatment status of each unit of
  assignment, and optionally block information. This is done using the
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  or optionally
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  and
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  functions.
- Fit both least squares and robust regression models to adjust for
  covariates and enhance the precision of treatment effect estimates.
- Use the
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  function to process the covariate adjustment model.
- Fit a model using the
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  function to estimate treatment effect that accounts for the
  specification information and the covariate adjustment by including
  the `weights` and `offset` arguments in the function.
