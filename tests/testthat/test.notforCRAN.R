test_that("#137 issues with grouped_df", {
  suppressMessages(library(tidyverse))
  data(simdata)

  # The failure we were seeing was due to the continuous moderator being renamed
  # through tidyverse. This should no longer error.
  mod <- lmitt(y ~ o, data = group_by(simdata, z), design = z ~ cluster(uoa1, uoa2))
  expect_s4_class(mod, "teeMod")
})


test_that("#174 non-data-frame input handled", {
  suppressMessages(library(tidyverse))
  data(simdata)

  expect_error(rct_design(z ~ uoa(uoa1, uoa2), data = y ~ x),
               "Failed to convert")

  sd_tibble <- as_tibble(simdata)

  des1 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  des2 <- rct_design(z ~ uoa(uoa1, uoa2), data = sd_tibble)
  des1@call <- des2@call
  expect_identical(des1, des2)
})
