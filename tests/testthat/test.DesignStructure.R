test_that("get_structure() fn", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  ss <- get_structure(des)

  expect_true(is(ss, "DesignStructure"))
  expect_identical(data.frame(ss), des@structure)
  expect_false(ss@binary)
  expect_identical(ss@Design, des)

  ss2 <- get_structure(des, binary = TRUE)
  # No dichotomy so should be identical
  expect_identical(data.frame(ss), data.frame(ss2))
  expect_true(ss2@binary)
  expect_identical(ss2@Design, des)


  des <- rct_design(o ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ss <- get_structure(des)

  expect_true(is(ss, "DesignStructure"))
  expect_identical(data.frame(ss), des@structure)
  expect_false(ss@binary)
  expect_identical(ss@Design, des)

  expect_error(get_structure(des, binary = TRUE),
               "No binary treatment")

  dichotomy(des) <- o <= 3 ~ .

  ss2 <- get_structure(des, binary = TRUE)
  expect_identical(ss[, -1], ss2[, -1])
  expect_identical(ss2[, 1], treatment(des, binary = TRUE)[, 1])
  expect_true(ss2@binary)
  expect_identical(ss2@Design, des)
})

test_that("show.Design.Structure", {

  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ss <- get_structure(des)
  expect_output(show(ss), "Treatment")
  expect_output(show(ss), "Cluster")
  expect_output(show(ss), "Block")

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)

  cap <- capture_output(show(get_structure(des)))
  expect_match(cap, "Treatment")
  expect_match(cap, "Unit of Assignment", fixed = TRUE)
  expect_no_match(cap, "Block")

  des <- rd_design(z ~ unitid(uoa1, uoa2) + forcing(force), data = simdata)

  cap <- capture_output(show(get_structure(des)))

  expect_match(cap, "Treatment")
  expect_match(cap, "Unit ID", fixed = TRUE)
  expect_no_match(cap, "Block")
  expect_match(cap, "Forcing")

})
