# Cluster-randomized experiment data on voter turnout in cable system markets

This dataset is a toy example derived from a cluster-randomized field
experiment that evaluates the effect of “Rock the Vote” TV
advertisements on voter turnout rate. The original study included 23,869
first-time voters across 85 cable television markets in 12 states. These
markets were grouped into matched sets based on their past voter turnout
rates and then randomly assigned to either a treatment or control
condition. This toy dataset is constructed by randomly sampling 10% of
individuals from selected cable television markets in the original
dataset.

## Usage

``` r
GV_data
```

## Format

A `data.frame` with 248 rows and 7 columns.

- age Age of participant

- vote_04 Outcome variable indicating whether participant voted

- tv_company Cable system serving participant's residential area

- treatment Binary variable denoting treatment assignment

- pairs A numeric indicator for the strata or matched pair group to
  which a cable system belongs (1-3)

- population_size Total population size of residential area served by
  cable system

- sample_size Number of individuals sampled from the cable system
  cluster

## Source

<https://isps.yale.edu/research/data/d005>

## Details

The original dataset was drawn from a randomized controlled trial in
which 85 cable system areas were first grouped into 40 matched sets
based on historical voter turnout. Within each matched set, one cable
system area was randomly assigned to the treatment condition, while the
others served as controls.

This toy dataset includes a subset of the original replication data,
specifically individuals from matched sets 1–3, which encompass 7 of the
85 cable system areas. Within these selected clusters, a 10% random
sample of individuals was taken.

The fuller Green-Vavreck dataset that this derives from bears a Creative
Commons BY-NC-ND license (v3.0) and is housed in Yale University's
Institution for Social and Policy Studies (ID: D005).

## References

Green, Donald P. & Lynn Vavreck (2008) "Analysis of Cluster-Randomized
Experiments: A Comparison of Alternative Estimation Approaches."
Political Analysis 16(2):138-152.
