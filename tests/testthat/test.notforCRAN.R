test_that("#137 issues with grouped_df", {
  suppressMessages(library(tidyverse))
  data(simdata)

  # The failure we were seeing was due to the continuous moderator being renamed
  # through tidyverse. This should no longer error.
  mod <- lmitt(y ~ o, data = group_by(simdata, z), design = z ~ cluster(uoa1, uoa2))
  expect_s4_class(mod, "teeMod")
})
