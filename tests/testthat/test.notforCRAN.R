test_that("#137 issues with grouped_df", {
  if (suppressMessages(require(tidyverse, quietly = TRUE))) {
    data(simdata)

    # The failure we were seeing was due to the continuous moderator being renamed
    # through tidyverse. This should no longer error.
    mod <- lmitt(y ~ o, data = group_by(simdata, z), specification = z ~ cluster(uoa1, uoa2))
    expect_s4_class(mod, "teeMod")
  } else {
    # Avoid warnings about empty tests
    expect_true(TRUE)
  }

})


test_that("#174 non-data-frame input handled", {
  if (suppressMessages(require(tidyverse, quietly = TRUE))) {
    data(simdata)

    expect_error(rct_spec(z ~ uoa(uoa1, uoa2), data = y ~ x),
                 "Failed to convert")

    sd_tibble <- as_tibble(simdata)

    spec1 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
    spec2 <- rct_spec(z ~ uoa(uoa1, uoa2), data = sd_tibble)
    spec1@call <- spec2@call
    expect_identical(spec1, spec2)
  } else {
    # Avoid warnings about empty tests
    expect_true(TRUE)
  }
})
