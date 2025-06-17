data(simdata)
inputs <- character(0)
storage <- character(0)
has_binary <- logical(0)

## numeric 0/1 treatment
inputs <- c(inputs, "Numeric 0/1 only")
spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## numeric non-0/1
inputs <- c(inputs, "Numeric NOT 0/1 only")
spec <- rct_spec(o ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## factor 0/1 levels
simdata$foo <- as.factor(simdata$z)
inputs <- c(inputs, "factor with only 0/1 levels")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## ordinal 0/1 levels
simdata$foo <- as.ordered(simdata$z)
inputs <- c(inputs, "ordinal with only 0/1 levels")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## factor non-0/1
simdata$foo <- as.factor(simdata$o)
inputs <- c(inputs, "factor NOT 0/1 levels only")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## ordinal non-0/1
simdata$foo <- as.ordered(simdata$o)
inputs <- c(inputs, "ordinal NOT 0/1 levels only")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## factor/ordinal with string levels
simdata$foo <- factor(ifelse(simdata$z == 1, "treatment", "control"))
inputs <- c(inputs, "factor with string levels (2 levels)")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

simdata$foo <- ordered(ifelse(simdata$o <= 1, "low", ifelse(simdata$o <= 3, "medium", "high")))
inputs <- c(inputs, "ordinal with string levels (3 levels)")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## character (2 levels)
simdata$foo <- ifelse(simdata$z == 1, "yes", "no")
inputs <- c(inputs, "character (2 levels)")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## character (more than 2 levels)
simdata$foo <- ifelse(simdata$o <= 1, "A", ifelse(simdata$o <= 2, "B", ifelse(simdata$o <= 3, "C", "D")))
inputs <- c(inputs, "character (more than 2 levels)")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## logical
simdata$foo <- as.logical(simdata$z)
inputs <- c(inputs, "logical")
spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

## conditionals (e.g. rct_spec(o > 2 ~ ...))
inputs <- c(inputs, "conditional (o > 2)")
spec <- rct_spec(o > 2 ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

inputs <- c(inputs, "conditional (z == 1)")
spec <- rct_spec(z == 1 ~ cluster(uoa1, uoa2), data = simdata)
storage <- c(storage, class(treatment(spec)[, 1])[1])
has_binary <- c(has_binary, has_binary_treatment(spec))

results <- cbind(inputs, storage, has_binary)
knitr::kable(results)
