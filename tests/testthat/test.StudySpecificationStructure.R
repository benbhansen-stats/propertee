test_that("get_structure() fn", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  ss <- get_structure(spec)

  expect_true(is(ss, "StudySpecificationStructure"))
  expect_identical(data.frame(ss), spec@structure)
  expect_identical(ss@StudySpecification, spec)

  ss2 <- get_structure(spec)
  # No dichotomy so should be identical
  expect_identical(data.frame(ss), data.frame(ss2))
  expect_identical(ss2@StudySpecification, spec)


  spec <- rct_spec(o ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ss <- get_structure(spec)

  expect_true(is(ss, "StudySpecificationStructure"))
  expect_identical(data.frame(ss), spec@structure)
  expect_identical(ss@StudySpecification, spec)
})

test_that("show.StudySpecification.Structure", {

  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ss <- get_structure(spec)
  expect_output(show(ss), "Treatment")
  expect_output(show(ss), "Cluster")
  expect_output(show(ss), "Block")

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)

  cap <- capture_output(show(get_structure(spec)))
  expect_match(cap, "Treatment")
  expect_match(cap, "Unit of Assignment", fixed = TRUE)
  expect_no_match(cap, "Block")

  spec <- rd_spec(z ~ unitid(uoa1, uoa2) + forcing(force), data = simdata)

  cap <- capture_output(show(get_structure(spec)))

  expect_match(cap, "Treatment")
  expect_match(cap, "Unit ID", fixed = TRUE)
  expect_no_match(cap, "Block")
  expect_match(cap, "Forcing")

})
