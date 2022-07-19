test_that("structure() fn", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  ss <- structure(des)

  expect_true(is(ss, "DesignStructure"))
  expect_identical(data.frame(ss), des@structure)
  expect_false(ss@binary)
  expect_identical(ss@Design, des)

  ss2 <- structure(des, binary = TRUE)
  # No dichotomy so should be identical
  expect_identical(data.frame(ss), data.frame(ss2))
  expect_true(ss2@binary)
  expect_identical(ss2@Design, des)


  des <- rct_design(o ~ cluster(cid1, cid2) + block(bid), data = simdata)
  ss <- structure(des)

  expect_true(is(ss, "DesignStructure"))
  expect_identical(data.frame(ss), des@structure)
  expect_false(ss@binary)
  expect_identical(ss@Design, des)

  expect_error(structure(des, binary = TRUE),
               "No binary treatment")

  dichotomy(des) <- o <= 3 ~ .

  ss2 <- structure(des, binary = TRUE)
  expect_identical(ss[, -1], ss2[, -1])
  expect_identical(ss2[, 1], treatment(des, binary = TRUE)[, 1])
  expect_true(ss2@binary)
  expect_identical(ss2@Design, des)
})

test_that("show.Design.Structure", {

  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  ss <- structure(des)
  expect_output(show(ss), "Treatment")
  expect_output(show(ss), "Cluster")
  expect_output(show(ss), "Block")

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)

  cap <- capture_output(show(structure(des)))
  expect_match(cap, "Treatment")
  expect_match(cap, "Unit of Assignment", fixed = TRUE)
  expect_no_match(cap, "Block")

  des <- rd_design(z ~ unitid(cid1, cid2) + forcing(force), data = simdata)

  cap <- capture_output(show(structure(des)))

  expect_match(cap, "Treatment")
  expect_match(cap, "Unit ID", fixed = TRUE)
  expect_no_match(cap, "Block")
  expect_match(cap, "Forcing")

})