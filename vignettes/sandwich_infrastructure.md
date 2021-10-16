---
title: 'Enumeration of data objects required for standard errors of prior covariance model-assisted effect estimation'
author: 'Ben Hansen'
date: 'October 2021'
---

## Context

Suppose an estimator chain beginning with fitting a covariance model,
`cmod` say, to a covariance sample, and ending with a model along the
lines of `lm(promotion ~ treat, data=qe_spl, weights=ate(des),
offset=cov_adj(cmod))` or `lm(promotion ~ treat*demographics,
data=qe_spl, weights=ate(des), offset=cov_adj(cmod))` --- a *direct
adjustment* model.  (Much or all of what's said here may apply also
with end models of more general forms, for example `lm(promotion ~
treat + cov_adj(cmod), data=qe_spl, weights=ate(des))` or
`lm(promotion ~ treat * cov_adj(cmod), <...>)`. But for now our focus
is covariance modeling followed by direct adjustment.)  The covariance
sample and quasiexperimental sample may be disjoint, identical or
overlapping.  This note describes computational infrastructure for
handling the second two cases, where errors from the fitting of the
covariance model need to be propagated into standard errors reported
for `treat` coefficients.

The direct adjustment model may be fit in a different context than the covariance model, with only the artifacts of that modeling that are stored in `cov_adj(cmod)`.  We determine this list, but we want it to be a minimal list, for storage, privacy and regulatory reasons.  It's imperative that it not contain student- or patient-level data, directly or implicitly, and it's preferable that it not contain data at all.

With reference to the formulas for stacked estimating equations of
Carroll, Ruppert, Stefanski and Crainiceanu (2006, p.373), the
covariance model has estimating equations $\phi(\tilde{\mathbf{Y}},
\alpha)$ with parameters $\alpha$, and Fisher information and estimating-equation covariance matrices $A_{11}$ and $B_{11}$ respectively; while the direct adjustment model's
are $\psi(\tilde{\mathbf{Y}}, \tau, \alpha)$, the `treat` coefficients
being $\tau$, with sandwich components $A_{22}$ and $B_{22}$. (The symbols "$\phi()$", "$\psi()$" and "$\alpha$" are
used as Carroll et al use them, while our "$\tau$" corresponds to
their "$\mathcal{B}$".)

The leading $n^{-1}$ factor should be struck from the expression for $\mathrm{var}(\hat{\mathcal{B}})$.  Carroll et al's formulas for $A_{11}$, $A_{21}$ and $A_{22}$ then apply, although design-based standard errors call for a different calculations of $B_{11}$, $B_{12}$ and $B_{22}$. In these expressions, $i$ is taken to range over clusters, element-wise estimating equation contributions being summed within clusters to form cluster-wise estimating equation contributions.  Denote clusters represented in covariance and quasiexperimental samples by $\mathcal{C}$ and $\mathcal{Q}$ respectively.

## Materials to have assembled when SE calculations are to be performed

To estimate variances and covariances of $\tau$, we'll need:

1. Sufficient information about $\mathcal{C}$ and $\mathcal{Q}$ to identify their intersection $\mathcal{C}\cap\mathcal{Q}$, as is needed to estimate $B_{21}$; 
2. Vectors of estimating
functions $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C}\cap\mathcal{Q}\}$
and $\{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, \hat{\alpha}): i \in \mathcal{C}\cap\mathcal{Q}\}$, as are needed for $B_{21}$; 
3. For the quasiexperimental sample $\mathcal{Q}$, matrices
$$\nabla_{\alpha} \{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, {\alpha}): i \in \mathcal{Q}\} \vert_{\alpha=\hat\alpha},$$ corresponding to $A_{21}$; 
4. Estimates of the direct adjustment model's A matrix $n A_{22} = \nabla_\tau \sum_{i \in \mathcal{Q}} \mathbb{E}[\psi(\tilde{\mathbf{Y}}_i, \tau, \hat{\alpha})] \vert_{\tau=\hat\tau}$, i.e. its unscaled Fisher information w.r.t. $\tau$ only, and unscaled B matrix $n B_{22} = \mathrm{Cov}[\sum_{i \in \mathcal{Q}} \psi(\tilde{\mathbf{Y}}_i, {\tau}, {\alpha})]_{({\tau}, {\alpha})=(\hat{\tau}, \hat{\alpha})}$; 
5. The covariance model's observed information matrix/empirical A matrix estimate $n \hat{A}_{11} = \sum_{i \in \mathcal{C}} \nabla_\alpha  [\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$; and
6. for covariance estimation in the conventional "model-based" setup only, estimates of the covariance model's B matrix $nB_{11} = \mathrm{Cov}[\sum_{i \in \mathcal{C}}\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$ (a "clustered" covariance estimate). 

In (5), observed information is preferred to "observed expected" information,
$\sum_{i \in \mathcal{C}} \nabla_\alpha \mathbb{E} [\phi(\tilde{\mathbf{Y}}_i; {\alpha})]_{\alpha=\hat\alpha}$, because observed information is agnostic as to whether expectation is calculated with conditioning on potential outcomes, ie the finite population perspective, or with conditioning on treatment assignment, the model based perspective. 

Regarding (6), $B_{11}$ is not needed for design-based standard errors, as in this setting observations outside of the quasiexperimental sample do not contribute to the covariance model's B matrix.  Only quasiexperimental sample observations do, and we'll have access to these when the direct adjustment model is fit. As we also have $\{\phi(\tilde{\mathbf{Y}}_i; 
\hat{\alpha}): i \in \mathcal{Q}\}$ and $\{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, \hat{\alpha}): i \in \mathcal{Q}\}$, we can use these materials to estimate $B_{12}$ and $B_{22}$.


## Software implementation comments on 1--6 above

1\. When there's a cluster/elements distinction, `cmod` will likely be aware of elements not clusters. That's natural for the direct adjustment model also.  We can make this an expectation: when we get to calculating covariances, we'll use the design object to identify element-cluster containment associations.

2\. Rather than passing $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C}\}$ forward along with `cmod`, it would be preferable to pass forward what's needed to reconstruct $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{C} \cap \mathcal{Q}\}$ from the combination of a person-record free summary of `cmod` with microdata on elements within clusters in $\mathcal{Q}$.  I propose that we serve these needs with a function that applies to fitted covariance models (or reductions of them that we might hide away and pass around as attributes of other things) and returns a function.  This returned function itself accepts a data frame as its argument, e.g. the quasiexperimental sample, and returns $\{\phi(\tilde{\mathbf{Y}}_i; \hat{\alpha}): i \in \mathcal{Q}\}$, the empirical estimating function of the covariance model, as returned by `sandwich::estfun()`, except evaluated over $\mathcal{Q}$, not whatever sample the covariance model was fit to.

3\. This leaves, $\nabla_{\alpha} \{\psi(\tilde{\mathbf{Y}}_i, \hat{\tau}, {\alpha}): i \in \mathcal{Q}\} \vert_{\alpha=\hat\alpha}$.  For $\psi()$'s that use only "predictions" of the covariance model, as ours does, the first derivative of the covariance model predictions as applied to data in $\mathcal{Q}$ will provide sufficient information from the covariance model to complete this calculation. Specifically, the derivative of whatever's returned by `optmatch::scores(cmod)`, or `predict(cmod, newdata=Q)`.  Our options include:

a. Expand the role of `scores()` to include furnishing first derivatives w.r.t. variables figuring in the prediction.  (Perhaps when given additional argument `deriv=TRUE`; `predict.smooth.spline()` does something like this.)
b.  Expect a fast predict method, and use it to figure numerical derivatives ofthe prediction.

I understand that `geex` does (b).  I'm attracted to (a).  The two approaches aren't necessarily incompatible, in that one could have dedicated `scores()` methods for specific model classes plus a generic method that figures the derivatives numerically. 

4\. We can take extant calculations of a direct adjustment model's information matrix, with the proviso that we keep track of whatever scaling those calcuations may have applied. For design-based SEs we'll need our own B matrix calculation.  For model-based SEs we can plug in to extant software, again with attention to how they're scaling things. (I have yet to consider whether use of HC0--3 etc for $B_{22}$ calls for corresponding adjustment to estimation of $B_{21}$ and/or $B_{11}$, nor whether heuristics animating these adjustments make sense in this context.)

5\. We can take extant calculations of a covariance model's information matrix, with the proviso that we keep track of whatever scaling those calculations may have applied.

6\. We can take extant calculations of a covariance model's B matrix, with the provisos that we keep track of whatever scaling those calculations may have applied and we look to the study design for clustering.  (I have yet to consider whether use of HC0--3 etc for $B_{11}$ calls for corresponding adjustment to estimation of $B_{21}$ and/or $B_{22}$, nor whether heuristics animating these adjustments make sense in this context.)

