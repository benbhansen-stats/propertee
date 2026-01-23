# Regression Discontinuity Design

## Data and StudySpecification

``` r
data(lsoSynth)
```

The data for this example were randomly simulated using the `synthpop`
package in `R` based on data originally collected by [Lindo, Sanders,
and Oreopoulos (2010; hereafter
LSO)](https://www.aeaweb.org/articles?id=10.1257/app.2.2.95). (The
original real data can be found
[here](https://www.openicpsr.org/openicpsr/project/113751/version/V1/view);
code to simulate the fake data can be found
[here](https://github.com/benbhansen-stats/propertee/blob/main/data-raw/makeRDDdata.r).)

At a major public university, students with a cumulative grade point
average (GPA) below a pre-specified cutoff are placed on academic
probation (AP). Students on AP are given extra support, and threatened
with suspension if their GPAs do not improve. See LSO for details. Along
with LSO, we will consider only students at the end of their first year
at the university. The university consisted of three campuses, and they
AP cutoff varied by campus. To simplify matters, we centered each
student’s first year GPA on their campus’s cutoff, creating a new
variable $R_{i} \equiv GPA_{i} - c_{campus{\lbrack i\rbrack}}$, where
$GPA_{i}$ is student $i$’s first-year GPA, and
$c_{campus{\lbrack i\rbrack}}$ is the AP cutoff for student $i$’s
campus. A student $i$ is placed on AP if $R_{i} < 0$ and avoids AP if
$R_{i} \geq 0$.

We will attempt to estimate the average effect of AP placement on
`nextGPA`, students’ subsequent GPA (either over the summer or in the
following academic semester). This figure plots each the average
`nextGPA` for students with distinct `R` values (i.e. centered
first-year GPA). The cutoff for AP, 0, is denoted with a dotted line.
The the sizes of the points are proportional to the natural log of the
numbers of students with each unique value of `R`.

``` r

figDat <- aggregate(lsoSynth[,c('nextGPA','lhsgrade_pct')],by=list(R=lsoSynth$R),
                    FUN=mean,na.rm=TRUE)
figDat$n <- as.vector(table(lsoSynth$R))

with(figDat,
     plot(R,nextGPA,cex=log(n+1)/3,main='Average Subsequent GPA\n as a function of 1st-Year GPA,\n Centered at AP Cutoff'))
abline(v=0,lty=2)
```

![](RDD_files/figure-html/rddPlots-1.png)

## Background

### Regression Discontinuity Design

What characterizes discontinuity (RD) design, including the LSO study,
is that treatment is assigned as a function of a numeric “running
variable” $R$, along with a prespecified cutoff value $c$, so that
treatment is assigned to subjects $i$ for whom $R_{i} > c$, or for whom
$R_{i} < c$. If $Z_{i}$ denote’s $i$’s treatment assignment, we write
$Z_{i} = \mathbf{1}\{ R_{i} < c\}$ or $\mathbf{1}\{ R_{i} < c\}$, where
$\mathbf{1}\{ x\}$ is the indicator function–equal to 1 when $x$ is true
and 0 otherwise. In the LSO study, $i$ indexes students, centered
first-year GPA \$R_i \$ is the running variable, and the treatment under
study is AP placement. Then \$ Z_i={R_i\<0} \$. (In “fuzzy” RD design,
$R$’s value relative to $c$ doesn’t completely determine treatment, but
merely affects the probability of treatment, so subjects with, say,
$R_{i} > c$ are more likely to be treated than those with $R_{i} < c$;
in fuzzy RD design, $Z_{i} = \mathbf{1}\{ R_{i} > c\}$ is typically
modeled as an instrument for treatment receipt.)

RD design hold a privileged place in causal inference, because unlike
most other observational designs, the mechanism for treatment assignment
is known. That is, $R$ is the only confounder. On the other hand, since
$Z$ is completely determined by $R$, there are no subjects with the same
values for $R$ but different values for $Z$, so common observational
study techniques, such as matching on $R$, are impossible. Instead, it
is necessary to adjust for $R$ with modeling–typically regression.

### Analyzing RD Design with ANCOVA

Traditionally, the typical way to analyze data from RD design is with
ANCOVA, by fitting a regression model such as
$$Y_{i} = \beta_{0} + \beta_{1}R_{i} + \beta_{2}Z_{i} + \epsilon_{i}$$
where $Y_{i}$ is the outcome of interest measured for subject $i$ (in
our example `nextGPA`), and $\epsilon_{i}$ is a random regression error.
Then, the regression coefficient for $Z$, $\beta_{2}$, is taken as an
estimate of the treatment effect. Common methodological updates to the
ANCOVA approach include an interaction term between $Z_{i}$ and $R_{i}$,
and the substitution semi-parametric regression, such as local linear or
polynomial models, for linear ordinary least squares.

Under suitable conditions, the ANCOVA model is said to estimate a “Local
Average Treatment Effect” (LATE). If $Y_{1}$ and $Y_{0}$ are the
potential outcomes for $Y$, then the LATE is defined as
$$\lim\limits_{r\rightarrow c^{+}}E\left( Y_{1}|R = r \right) - \lim\limits_{r\rightarrow c^{-}}E\left( Y_{0}|R = r \right)$$
or equivalently
$$\lim\limits_{\Delta\rightarrow 0^{+}}E\left( Y_{1} - Y_{0}|R \in (c - \Delta,c + \Delta) \right)$$
where $E$ denotes expectation–that is, the LATE is the limit of average
treatment effects for subjects with $R$ in ever-shrinking regions around
$c$.

Among other considerations, the LATE target suggests that data analysts
restrict their attention to subjects with $R$ falling within a bandwidth
$b > 0$ of $c$, e.g. only fitting the model to subjects $i$ with within
a “window of analysis” $\mathcal{W} = \{ i:R_{i} \in (c - b,c + b)\}$. A
number of methods have been proposed to select $b$, including cross
validation, non-parametric modeling of the second derivative of
$f(r) = E\left( Y|R = r \right)$, and specification tests.

### The **propertee** Approach to RD Design

The **propertee** approach to RD design breaks the data analysis into
three steps:

1.  Conduct design tests and choose a bandwidth $b > 0$, along with,
    possibly, other data exclusions, resulting in an analysis sample
    $W$.
2.  Fit a covariance model $Y_{i0} = g\left( R_{i},x_{i};\beta \right)$,
    modeling $Y_{i0}$ as a function of running variable $R_{i}$ and
    (optionally) other covariates $\mathbf{x}_{i}$. Fitting a covariance
    model to only the control subjects, as in other design, would entail
    extrapolation of a model fit to subjects with $R_{i} \in (c - b,c)$
    to subjects with $R_{i} \in (c,c + b)$, or vice-versa, which is
    undesirable. Instead, we fit the covariance model to the full
    analysis sample $\mathcal{W}$, including both treated and untreated
    subjects. However, since we are interested in modeling $Y_{0}$, and
    not $Y_{1}$ or $Y$, so rather than fitting $g( \cdot )$, we fit an
    extended model
    $$\widetilde{g}\left( R_{i},x_{i},Z_{i};\beta,\gamma \right) = g\left( R_{i},x_{i};\beta \right) + \gamma Z_{i}$$
    including a term for treatment assignment. The estimate for $\gamma$
    is, in essence, a provisional estimate of the treatment effect.
3.  Let
    ${\widehat{Y}}_{i0} = g\left( R_{i},x_{i};\widehat{\beta} \right)$,
    using only the model for $Y_{0}$—not including the term for
    $Z$—along with $\widehat{\beta}$ estimated in step 2. Then estimate
    the average treatment effect for subjects in $W$ using the
    difference in means estimator:
    $$d_{\widehat{\beta}} = \frac{\sum\limits_{i \in W}Z_{i}\left( Y_{i} - {\widehat{Y}}_{i0} \right)}{\sum\limits_{i \in W}Z_{i}} - \frac{\sum\limits_{i \in W}\left( 1 - Z_{i} \right)\left( Y_{i} - {\widehat{Y}}_{i0} \right)}{\sum\limits_{i \in W}\left( 1 - Z_{i} \right)}$$

The estimator $d_{\widehat{\beta}}$ will be consistent if the model
$g\left( R_{i},x_{i};\beta_{0} \right)$, with $\beta_{0}$ as the
probability limit for $\widehat{\beta}$, successfully removes all
confounding due to $R$,
i.e. $Y_{0} - g\left( R_{i},x_{i};\beta_{0} \right)\bot\!\!\!\bot Z$ for
$i \in \mathcal{W}$.

## Analyzing an RD design in **propertee**

### Determining the window of analysis

There is an extensive literature on determining the appropriate window
of analysis for an RD design See, for instance, Imbens and Kalyanaraman
(2012) and Sales and Hansen (2020). This stage of the analysis is beyond
the scope of this vignette, and does not require the propertee package.
For the purpose of this example, we will focus our analysis on
$\mathcal{W} = \{ i:R_{i} \in \lbrack - 0.5,0.5\rbrack\}$.

### Initializing the RD Design Object

In general, **propertee** specification objects require users to specify
a unit of assignment variable, that corresponds to the level of
treatment assignment–in other words, which units were assigned between
conditions. This is necessary both for combining results across
different models or datasets, and for estimating correct
specification-based standard errors. In `lsoSynth`, there are no
clusters, and individual students are assigned to academic probation on
an individual basis. The dataset does not include an ID variable, but
any variable that takes a unique value for each row will do. We will use
row-names.

``` r
lsoSynth$id <- rownames(lsoSynth)
```

Defining an RD design requires, at minimum, identifying the running
variable(s) $R$, as well as how $R$ determines treatment assignment. In
our example,

``` r
lsoSynth$Z <- lsoSynth$R<0

lsoSynthW <- subset(lsoSynth,abs(R)<=0.5)

#lsoSynthW$id <- 1:nrow(lsoSynthW)
spec <- rd_spec(Z ~ forcing(R) + unitid(id), data=lsoSynth, subset=abs(lsoSynth$R) <= 0.5)
```

### Modeling $Y_{C}$ as a function of $R$

We will consider two potential models $\widetilde{g}( \cdot )$ of
$Y_{C}$ as a function of $R$. First, the standard OLS model:

``` r
### this doesn't work:
g1 <- lm(nextGPA ~ R + Z, data = lsoSynth, weights = ett(spec)
```

``` r
##this works, but it's annoying to enter in the subset expression a 2nd time:
g1 <- lm(nextGPA ~ R + Z, data = lsoSynth, subset = abs(R) <= 0.5)
```

The second is a bounded-influence polynomial model of the type
recommended in Sales & Hansen (2020):

``` r
g2 <-
if(requireNamespace("robustbase", quietly = TRUE)){
  robustbase::lmrob(nextGPA~poly(R,5)+Z,data=lsoSynthW)
} else g1
```

\[Put something here evaluating them?\]

### Estimating Effects

``` r

yhat1 <- predict(g1,data.frame(R=forcings(spec)[[1]],Z=FALSE))
yhat2 <- predict(g2,data.frame(R=forcings(spec)[[1]],Z=FALSE))

plot(yhat1,yhat2)
```

![](RDD_files/figure-html/yhat-1.png)

``` r
### method 1:

mean(lsoSynthW$nextGPA[lsoSynthW$Z]-yhat1[lsoSynthW$Z])-
  mean(lsoSynthW$nextGPA[!lsoSynthW$Z]-yhat1[!lsoSynthW$Z])
#> [1] 0.2553202

#### method 2:
coef(lm(nextGPA~Z, offset=yhat1,data=lsoSynthW))['ZTRUE']
#>     ZTRUE 
#> 0.2553202
```

``` r
### method 1:

mean(lsoSynthW$nextGPA[lsoSynthW$Z]-yhat2[lsoSynthW$Z])-
  mean(lsoSynthW$nextGPA[!lsoSynthW$Z]-yhat2[!lsoSynthW$Z])
#> [1] 0.2295554

#### method 2:
coef(lm(nextGPA~Z, offset=yhat2,data=lsoSynthW))['ZTRUE']
#>     ZTRUE 
#> 0.2295554
```
