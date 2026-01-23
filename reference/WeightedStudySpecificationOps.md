# `WeightedStudySpecification` Operations

Algebraic operators on `WeightedStudySpecification` objects and numeric
vectors. `WeightedStudySpecification`s do not support addition or
subtraction.

## Usage

``` r
# S4 method for class 'WeightedStudySpecification,numeric'
e1 + e2

# S4 method for class 'numeric,WeightedStudySpecification'
e1 + e2

# S4 method for class 'WeightedStudySpecification,numeric'
e1 - e2

# S4 method for class 'numeric,WeightedStudySpecification'
e1 - e2

# S4 method for class 'WeightedStudySpecification,numeric'
e1 * e2

# S4 method for class 'numeric,WeightedStudySpecification'
e1 * e2

# S4 method for class 'WeightedStudySpecification,numeric'
e1/e2

# S4 method for class 'numeric,WeightedStudySpecification'
e1/e2
```

## Arguments

- e1, e2:

  `WeightedStudySpecification` or `numeric` objects

## Value

a `WeightedStudySpecification` object

## Details

These are primarily used to either combine weights via multiplication,
or to invert weights. Addition and subtraction are not supported and
will produce errors.
