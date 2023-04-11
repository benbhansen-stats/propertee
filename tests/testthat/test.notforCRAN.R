test_that("expand.model.frame changes", {
  # Since the .expand.model.frame.DA is a clone of stats::expand.model.frame,
  # this "test" captures if stats::expand.model.frame changes. See issue #102.

  fcode <- capture.output(flexida:::.expand.model.frame.DA)
  scode <- capture.output(stats::expand.model.frame)

  # Remove environment flags
  scode <- scode[!grepl("bytecode:", scode)]
  fcode <- fcode[!grepl("bytecode:", fcode)]
  scode <- scode[-length(scode)]
  fcode <- fcode[-length(fcode)]

  # Strip changes
  fcode <- fcode[!grepl("JE\\ addition", fcode)]
  fcode <- sub("\\#.*", "", fcode)
  fcode <- gsub("stats::", "", fcode)

  # Collapse to remove all whitespace differences
  fcode <- paste(fcode, collapse = "")
  scode <- paste(scode, collapse = "")

  # Strip preceeding or trailing spaces
  fcode <- gsub("[[:space:]]+", "", fcode)
  scode <- gsub("[[:space:]]+", "", scode)

  # For debugging purposes, this will do a character-wise comparison:
  ## print(Reduce(`==`, strsplit(c(fcode, scode), "")))

  if (fcode != scode) {
    stop("stats::expand.model.frame may have been changed! See Issue #102")
  }

  expect_true(TRUE) # To avoid empty test warnings
})