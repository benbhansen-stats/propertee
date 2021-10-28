---
title: 'Data structures to propagate error from earlier regressions used as inputs to subsequent regressions'
author: 'Ben Hansen'
date: 'October 2021'
---

## Context

Suppose an estimator chain beginning with fitting a covariance model,
`covmod` say, to a covariance sample, and ending with a model along the
lines of `lm(promotion ~ treat, data=qe_spl, weights=ate(des),
offset=cov_adj(covmod))` or `lm(promotion ~ treat*demographics,
data=qe_spl, weights=ate(des), offset=cov_adj(covmod))` --- a *direct
adjustment* model.  (Much or all of what's said here may apply also
with end models of more general forms, for example `lm(promotion ~
treat + cov_adj(covmod), data=qe_spl, weights=ate(des))` or
`lm(promotion ~ treat * cov_adj(covmod), <...>)`. But for now our focus
is covariance modeling followed by direct adjustment.)  The covariance
sample and quasiexperimental sample may be disjoint, identical or
overlapping.  This note describes computational infrastructure for
handling the second two cases, where errors from the fitting of the
covariance model need to be propagated into standard errors reported
for `treat` coefficients.

The direct adjustment model may be fit in a different context than the covariance model, using only the artifacts of that modeling that are stored in `cov_adj(covmod)`.  We determine this list, but we want it to be a minimal list, for storage, privacy and regulatory reasons.  It's imperative that it not contain student- or patient-level data, directly or implicitly, and it's preferable that it not contain data at all.


## Proposal

(Class & function names are tentative.)

- S4 classes `SandwichLayer` and `VeganLayer`, carrying certain extracts of the fit of a linear or asymptotically linear model. When the predictions of this model, model A, are used as independent variables in the fitting of a subsequent linear or asymptotically linear model, model B, objects and data carried in the `SandwichLayer` derived from model A will suffice to estimate the variance of errors of estimation in model B's parameters that is propagated down from estimation of model A. For variance propagation under the design-based perspective, it suffices to use a `VeganLayer`, i.e. a `SandwichLayer` minus the meat.^[A `SandwichLayer` object deriving from model A carries estimates of model A's "bread" and "meat" matrices, along with additional material relating to predictions and estimating functions of model A that will serve to connect it to the next "deck" of a multi-decker sandwich.  `VeganLayer` has all of this minus the meat.] The `SandwichLayer` and `VeganLayer` classes do not use the `.Data` slot.
- An S3 generic maker of `VeganLayers`, with methods for `glm`, ... A separate `add_meat()` function accepting a `glm` or other fitted model and a `VeganLayer`, transforming the `VeganLayer` into a `SandwichLayer`.  
- S4 classes `SandwichLayerWithPreds`/`VeganLayerWithPreds`, class extensions of `SandwichLayer` and `VeganLayer`.  A {`Sandwich`/`Vegan`}`LayerWithPreds` is a {`Sandwich`/`Vegan`}`Layer` with a `.Data` slot added on, to be occupied by a numeric vector (so that it extends the numeric base type). 
- A maker of `SandwichLayerWithPreds` and `VeganLayerWithPreds` objects, a function of the {`Sandwich`/`Vegan`}`Layer` and a data frame carrying variables needed for prediction from model A.  This presumes that the {`Sandwich`/`Vegan`}`Layer` carry what's needed to generate model A's predictions for a new data set.^[Predictions from model A aren't themselves needed to connect model A's sandwich layer to the next deck of the sandwich.  But their derivatives (with respect to each of model A's parameters) are required, and once we're taking responsibility for those the predictions themselves aren't much more to add.  Plus, by bundling these in the {`Sandwich`/`Vegan`}`Layer` without also encasing microdata in there, storing a {`Sandwich`/`Vegan`}`Layer` becomes a means of storing what's needed for predictions from model A while respecting privacy and data sharing agreement constraints.]
- `cov_adj()` invokes the maker function described at the last step, returning a {`Sandwich`/`Vegan`}`LayerWithPreds`, the predictions of which are on the response scale.  (`optmatch::scores()` might be changed to do the same, but returning the predictions on the scale of the linear predictor.) 


As a result, in the model frame stored with model B, the column of predictions deriving from model A has type {`Sandwich`/`Vegan`}`LayerWithPreds`, and carries the needed info from fitting of model A to calculated sandwich SEs for the multi-decker sandwich.  Two more thing will facilitate folding our calculations into existing R workflows:

- `as.multi_decker_lm()`, a conversion function accepting an object `x` inheriting from `lm`,  identifying {`Sandwich`/`Vegan`}`LayerWithPreds` objects among its predictor variables, assembling the associated {`Sandwich`/`Vegan`}`Layer` objects into a list that's tacked on to `x`.  The class attribute ot `x` is then prepended with `multi_decker_lm`, and the result is returned.
- Method `vcov.multi_decker_lm()`


## Basis in known extensions of Huber-White setup to chained estimators

With reference to the formulas for stacked estimating equations of
Carroll, Ruppert, Stefanski and Crainiceanu (2006, p.373), the
covariance model has estimating equations $\phi(\tilde{\mathbf{Y}},
\alpha)$ with parameters $\alpha$, and Fisher information and estimating-equation covariance matrices $A_{11}$ and $B_{11}$ respectively; while the direct adjustment model's
are $\psi(\tilde{\mathbf{Y}}, \tau, \alpha)$, the `treat` coefficients
being $\tau$, with sandwich components $A_{22}$ and $B_{22}$. (The symbols "$\phi()$", "$\psi()$" and "$\alpha$" are
used as Carroll et al use them, while our "$\tau$" corresponds to
their "$\mathcal{B}$".)

The leading $n^{-1}$ factor should be struck from the expression for $\mathrm{var}(\hat{\mathcal{B}})$.  Carroll et al's formulas for $A_{11}$, $A_{21}$ and $A_{22}$ then apply, although design-based standard errors call for a different calculations of $B_{11}$, $B_{12}$ and $B_{22}$. With the understanding that we don't apply Carroll et al's formulas for B matrices in terms of estimating functions^[Carroll et al's B-formulas are for model-based SEs without clustering; even when we're calculating model-based and not design-based SEs, we'll need to aggregate (by summing) the elementwise estimating equation contributions to the cluster level, to form cluster-wise estimating equation contributions, before applying their formulas.], we can then take $i$ to range over elements not clusters.    Denote elements represented in covariance and quasiexperimental samples by $\mathcal{C}$ and $\mathcal{Q}$ respectively.

## Required materials for SE calculations

To estimate variances and covariances of $\tau$, we'll need to assemble the following materials. 

1. Sufficient information about $\mathcal{C}$ and $\mathcal{Q}$ to identify their intersection $\mathcal{C}\cap\mathcal{Q}$, as is needed to estimate $B_{21}$; 
2. Vectors of estimating
functions $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C}\cap\mathcal{Q}\}$
and $\{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, \hat{\alpha}): i \in \mathcal{C}\cap\mathcal{Q}\}$, as are needed for $B_{21}$, along with element-cluster containment relationships; 
3. For the quasiexperimental sample $\mathcal{Q}$, matrices
$$\nabla_{\alpha} \{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, {\alpha}): i \in \mathcal{Q}\} \vert_{\alpha=\hat\alpha},$$ corresponding to $A_{21}$; 
4. Estimates of the direct adjustment model's A matrix $n A_{22} = \nabla_\tau \sum_{i \in \mathcal{Q}} \mathbb{E}[\psi(\tilde{\mathbf{Y}}_i, \tau, \hat{\alpha})] \vert_{\tau=\hat\tau}$, i.e. its unscaled Fisher information w.r.t. $\tau$ only, and unscaled B matrix $n B_{22} = \mathrm{Cov}[\sum_{i \in \mathcal{Q}} \psi(\tilde{\mathbf{Y}}_i, {\tau}, {\alpha})]_{({\tau}, {\alpha})=(\hat{\tau}, \hat{\alpha})}$; 
5. The covariance model's observed information matrix/empirical A matrix estimate $n \hat{A}_{11} = \sum_{i \in \mathcal{C}} \nabla_\alpha  [\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$; and
6. for covariance estimation in the conventional "model-based" setup only, estimates of the covariance model's B matrix $nB_{11} = \mathrm{Cov}[\sum_{i \in \mathcal{C}}\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$ (a "clustered" covariance estimate). 

In (5), observed information is preferred to "observed expected" information,
$\sum_{i \in \mathcal{C}} \nabla_\alpha \mathbb{E} [\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$, because observed information is agnostic as to whether expectation is calculated with conditioning on potential outcomes, ie the finite population perspective, or with conditioning on treatment assignment, the model based perspective.   In the special case of quantile regression^[Strictly speaking, the estimating equations of a quantile regression aren't differentiable in $\alpha$. All directional derivatives will exist, and I'm expecting this to be enough for our purposes, but I haven't thought it through.], observed information isn't ordinarily used in standard error calculations, and it may take some doing to get.    

Regarding (6), $B_{11}$ is not needed for design-based standard errors, as in this setting observations outside of the quasiexperimental sample do not contribute to the covariance model's B matrix.  Only quasiexperimental sample observations do, and we'll have access to these when the direct adjustment model is fit. As we also have $\{\phi(\tilde{\mathbf{Y}}_i; 
\hat{\alpha}): i \in \mathcal{Q}\}$ and $\{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, \hat{\alpha}): i \in \mathcal{Q}\}$, we can use these materials to estimate $B_{12}$ and $B_{22}$.


## Software implementation comments on 1--6 above, including contents of {`Sandwich`/`Vegan`}`Layer` objects

1\. A {`Sandwich`/`Vegan`}`Layer` object must identify $\mathcal{Q}$, by a key variable table defaulting to the row names of whatever data contributed to the fit of the model being summarized in this Layer (as an *n* by 1 data frame), for consistency with non-default situations where there may be multiple key variables.  

When there's a cluster/elements distinction, `covmod` will likely be aware of elements but not clusters. That's one reason that $\mathcal{Q}$ is identified at the element level; if it were a collection of clusters, it wouldn't be possible to identify it from `covmod` alone. When we get to calculating covariances, e.g. $B_{21}$, we'll use the design object to identify element-cluster containment associations.  

2\. The {`Sandwich`/`Vegan`}`Layer` classes also carry an `@estfun_maker` slot for a specialized function, a function of data frames returning numeric matrices of estimating functions that align with those data frames.

Rather than passing $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C}\}$ forward along with `covmod`, it is preferable^[This will help to address RQ2 of the project proposal.] to pass forward what's needed to reconstruct $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C} \cap \mathcal{Q}\}$ from the combination of a person-record free summary of `covmod` with microdata on elements within clusters in $\mathcal{Q}$.  I propose that we serve these needs with a function that applies to fitted covariance models (or reductions of them that we might hide away and pass around as attributes of other things) and returns a function.  This returned function itself accepts a data frame as its argument, e.g. the quasiexperimental sample $\mathcal{Q}$, and returns $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{Q}\}$, the empirical estimating function of the covariance model, as returned by `sandwich::estfun()`, except evaluated over $\mathcal{Q}$, not whatever sample the covariance model was fit to.

Element-cluster associations are needed, to organize these element-wise estimating function contributions prior to estimation of $B_{21}$; it's assumed that there'll be a Design object present that we can gather this material from. 

3\. The {`Sandwich`/`Vegan`}`LayerWithPreds` class extends {`Sandwich`/`Vegan`}`Layer` not ony with the inclusion of an `@.Data` slot carrying an numeric vector of predictions, but also an `@prediction_gradient` slot for a numeric matrix.  This matrix has as many rows as there are entries in the `.Data` vector, and as many columns as there are estimating equations.

The `@prediction_gradient` slot carries $\{\nabla_\alpha f_{\alpha}(\tilde{\mathbf{Y}}_i)|_\alpha=\hat\alpha, i \in \mathcal{Q}\}$, where $f_\alpha(\mathbf{Y})$ represents the prediction for data $\mathbf{Y}$ from a fitted model of `class(cov_mod)` with parameters $\alpha$. For $\psi()$'s that use only "predictions" of the covariance model, as ours does, the first derivative of the covariance model predictions as applied to data in $\mathcal{Q}$ will provide sufficient information from the covariance model to complete the calculation of  $\{\nabla_{\alpha}  \psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, {\alpha}) \vert_{\alpha=\hat\alpha}: i \in \mathcal{Q}\}$. 

(To calculate these gradients in `cov_adj()`, there are two options.

a. Expand the role of `scores()` to include furnishing first derivatives w.r.t. variables figuring in the prediction.  (Perhaps when given additional argument `deriv=TRUE`; `predict.smooth.spline()` does something like this.)
b.  Expect to find a fast `predict` method for the model supplied as an argument to `cov_adj()`; use it to figure numerical derivatives of the prediction.

I understand that `geex` does (b).  I'm attracted to (a).  The two approaches aren't necessarily incompatible, in that one could have dedicated `scores()` methods for specific model classes plus a generic method that figures the derivatives numerically.)

4\. The {`Sandwich`/`Vegan`}`Layer` classes aren't implicated in (4).  We can take extant calculations of a direct adjustment model's information matrix, with the proviso that we keep track of whatever scaling those calculations may have applied. For design-based SEs we'll need our own B matrix calculation.  For model-based SEs we can plug in to extant software, again with attention to how they're scaling things. (I have yet to consider whether use of HC0--3 etc for $B_{22}$ calls for corresponding adjustment to estimation of $B_{21}$ and/or $B_{11}$, nor whether heuristics animating these adjustments make sense in this context.)

5\. {`Sandwich`/`Vegan`}`Layer` classes carry a slot for the covariance models's observed information matrix.  The rows and columns of this matrix align with the columns of estimating functions that are produced by the function residing in the `@estfun_maker` slot.

We can take extant calculations of a covariance model's information matrix, with the proviso that we keep track of whatever scaling those calculations may have applied.

6\. The `VeganLayer` class carries a slot for an estimated B matrix.  `SandwichLayer` does not carry this slot, and that's the whole of the distinction between the `VeganLayer` and `SandwichLayer` classes.

We can take extant calculations of a covariance model's B matrix, with the provisos that we keep track of whatever scaling those calculations may have applied and we look to the study design for clustering.  (I have yet to consider whether use of HC0--3 etc for $B_{11}$ calls for corresponding adjustment to estimation of $B_{21}$ and/or $B_{22}$, nor whether heuristics animating these adjustments make sense in this context.)

