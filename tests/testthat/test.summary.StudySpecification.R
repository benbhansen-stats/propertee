test_that("summary.StudySpecification", {
  data(simdata)

  specrct <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  specrd  <-  rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                       data = simdata)
  specobs <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  expect_s3_class(summary(specrct), "summary.StudySpecification")
  expect_s3_class(summary(specrd) , "summary.StudySpecification")
  expect_s3_class(summary(specobs), "summary.StudySpecification")

  expect_s4_class(summary(specrct)$StudySpecification, "StudySpecification")
  expect_s4_class(summary(specrd)$StudySpecification , "StudySpecification")
  expect_s4_class(summary(specobs)$StudySpecification, "StudySpecification")

  expect_output(print(summary(specrct)), "Randomized")
  expect_output(print(summary(specrd)) , "Discontinuity")
  expect_output(print(summary(specobs)), "Observational")

  specrct <- rct_spec(o ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  expect_output(print(summary(specrct)), "...")
  expect_output(print(summary(specrct)), "excluded")
  expect_identical(specrct, summary(specrct)$StudySpecification)
  expect_s3_class(summary(specrct)$treatment_table, "table")

  expect_output(print(summary(specrct)),
                "group excluded")
  expect_output(print(summary(specrct), max_unit_print = 2),
                "groups excluded")

})
