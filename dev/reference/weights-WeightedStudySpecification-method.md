# Extract Weights from `WeightedStudySpecification`

A `WeightedStudySpecification` object contains a numeric vector with a
few additional slots, this extracts only the numeric vector.

## Usage

``` r
# S4 method for class 'WeightedStudySpecification'
weights(object, ...)
```

## Arguments

- object:

  a `WeightedStudySpecification` object

- ...:

  Ignored

## Value

A numeric `vector` of the weights
