# (Internal) Modeling weights with an accompanying StudySpecification

(Internal) Modeling weights with an accompanying StudySpecification

## Details

`@target` is used for calculation purpose; defining what weight to
calculate `@weightAlias` is only to store the alias used in creation of
the weights in case we want to report it later.

## Slots

- `.Data`:

  numeric vector of modeling weights

- `StudySpecification`:

  a StudySpecification

- `target`:

  character string, e.g. "ate"

- `weightAlias`:

  alias for target appearing in an originating call

- `dichotomy`:

  formula describing a treatment/comparison dichotomy
