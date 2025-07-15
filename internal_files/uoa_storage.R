data(simdata)
inputs <- character(0)
uoa1_storage <- character(0)
uoa2_storage <- character(0)

## numeric + numeric
inputs <- c(inputs, "numeric + numeric")
spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## numeric + factor
simdata$foo2 <- as.factor(simdata$uoa2)
inputs <- c(inputs, "numeric + factor")
spec <- rct_spec(z ~ uoa(uoa1, foo2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## factor + numeric
simdata$foo1 <- as.factor(simdata$uoa1)
inputs <- c(inputs, "factor + numeric")
spec <- rct_spec(z ~ uoa(foo1, uoa2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## factor + factor
inputs <- c(inputs, "factor + factor")
spec <- rct_spec(z ~ uoa(foo1, foo2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## numeric + character
simdata$char2 <- as.character(simdata$uoa2)
inputs <- c(inputs, "numeric + character")
spec <- rct_spec(z ~ uoa(uoa1, char2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## character + numeric
simdata$char1 <- as.character(simdata$uoa1)
inputs <- c(inputs, "character + numeric")
spec <- rct_spec(z ~ uoa(char1, uoa2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## character + character
inputs <- c(inputs, "character + character")
spec <- rct_spec(z ~ uoa(char1, char2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## ordered + factor
simdata$ord1 <- as.ordered(simdata$uoa1)
inputs <- c(inputs, "ordered + factor")
spec <- rct_spec(z ~ uoa(ord1, foo2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## factor + ordered
simdata$ord2 <- as.ordered(simdata$uoa2)
inputs <- c(inputs, "factor + ordered")
spec <- rct_spec(z ~ uoa(foo1, ord2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## logical + numeric
simdata$log1 <- as.logical(simdata$uoa1 %% 2)
inputs <- c(inputs, "logical + numeric")
spec <- rct_spec(z ~ uoa(log1, uoa2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## numeric + logical
simdata$log2 <- as.logical(simdata$uoa2 %% 2)
inputs <- c(inputs, "numeric + logical")
spec <- rct_spec(z ~ uoa(uoa1, log2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

## character + factor
inputs <- c(inputs, "character + factor")
spec <- rct_spec(z ~ uoa(char1, foo2), data = simdata)
uoa_data <- units_of_assignment(spec)
uoa1_storage <- c(uoa1_storage, class(uoa_data[, 1])[1])
uoa2_storage <- c(uoa2_storage, class(uoa_data[, 2])[1])

results <- cbind(inputs, uoa1_storage, uoa2_storage)
knitr::kable(results)
