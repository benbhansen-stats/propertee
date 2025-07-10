test_that("StudySpecification formula checking", {
  expect_true(.check_spec_formula(y ~ cluster(x)))
  expect_true(.check_spec_formula(y ~ cluster(x, z, q, r)))
  expect_true(.check_spec_formula(y ~ unitid(x)))
  expect_true(.check_spec_formula(y ~ unitid(x, z, q, r)))
  expect_true(.check_spec_formula(y ~ unit_of_assignment(x)))
  expect_true(.check_spec_formula(y ~ unit_of_assignment(x, z, q, r)))
  expect_true(.check_spec_formula(y ~ uoa(x)))
  expect_true(.check_spec_formula(y ~ uoa(x, z, q, r)))

  expect_error(.check_spec_formula(~ cluster(x)),
               "treatment")

  expect_error(.check_spec_formula(y ~ cluster(x) + unitid(z)),
               "Only one of")

  expect_error(.check_spec_formula(y ~ cluster(x) + cluster(z)),
               "Only one instance of `cluster")

  expect_error(.check_spec_formula(y ~ unitid(x) + unitid(z)),
               "Only one instance of `unitid")

  expect_error(.check_spec_formula(y ~ unit_of_assignment(x) +
                                       unit_of_assignment(z)),
               "Only one instance of `unit_of")

  expect_error(.check_spec_formula(y ~ uoa(x) + uoa(z)),
               "Only one instance of `unit_of")

  expect_true(.check_spec_formula(y ~ cluster(x) + block(z)))
  expect_true(.check_spec_formula(y ~ cluster(x) + block(z, a, b, c)))

  expect_error(.check_spec_formula(y ~ cluster(x) + block(z) + block(q)),
               "only one block")

  expect_error(.check_spec_formula(y ~ cluster(x) + forcing(z)),
               "only allowed")

  expect_true(.check_spec_formula(y ~ cluster(x) + forcing(z),
                                 allow_forcing = TRUE))

  expect_error(.check_spec_formula(y ~ cluster(x) + forcing(z) + forcing(q),
                                  allow_forcing = TRUE),
               "only one forcing")

})

test_that("binary treatment and dichotomy", {
  expect_error(has_binary_treatment(1),
               "must be")

  spec1 <- obs_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  spec2 <- obs_spec(o ~ uoa(uoa1, uoa2), data = simdata)

  expect_true(has_binary_treatment(spec1))

  expect_false(has_binary_treatment(spec2))
})

test_that("identical_StudySpecifications function", {
  data(simdata)

  spec1 <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  spec2 <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  spec3 <- rct_spec(dose ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  expect_true(identical_StudySpecifications(spec1, spec2))
  expect_false(identical_StudySpecifications(spec1, spec3))
  expect_false(identical_StudySpecifications(spec2, spec3))

})

test_that("identify_small_blocks", {
  spec_data <- data.frame(uid = letters[1:12],
                         bid = rep(LETTERS[1:3], each = 4),
                         a = c(rep(rep(c(0, 1), each = 2), 2), 0, rep(1, 3)))
  spec <- rct_spec(a ~ unitid(uid) + block(bid), spec_data)

  small_blocks <- identify_small_blocks(spec)
  expect_equal(length(small_blocks), 3)
  expect_equal(names(small_blocks), LETTERS[1:3])
  expect_true(all.equal(small_blocks, c(FALSE, FALSE, TRUE), check.attributes = FALSE))
})

test_that(".make_uoa_cluster_df errors", {
  expect_error(.make_uoa_cluster_df("spec"), "valid `StudySpecification` object")

  data(simdata)
  simdata_copy <- simdata
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata_copy)

  simdata_copy$uoa2[simdata_copy$uoa1 == 1 & simdata_copy$uoa2 == 2] <- 1
  expect_warning(.make_uoa_cluster_df(spec), "Some units of assignment")

  simdata_copy <- NULL
  expect_error(.make_uoa_cluster_df(spec), "specification data")

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  expect_error(.make_uoa_cluster_df(spec, "id"), "id column in the specification data")
})

test_that(".make_uoa_cluster_df", {
  spec_data <- data.frame(id = letters[1:20],
                         uid1 = rep(letters[1:4], each = 5),
                         uid2 = rep(letters[1:4], each = 5),
                         bid = rep(LETTERS[1:2], each = 10),
                         a = rep(rep(c(0, 1), each = 5), 2))

  # test specification creation using `cluster()`, `unitid()`, and `unit_of_assignment()`
  cluster_spec <- rct_spec(a ~ cluster(uid1), spec_data)
  unitid_spec <- rct_spec(a ~ unitid(uid1), spec_data)
  uoa_spec <- rct_spec(a ~ unit_of_assignment(uid1), spec_data)

  expect_true(all.equal(uc_df <- .make_uoa_cluster_df(cluster_spec),
                        .make_uoa_cluster_df(unitid_spec),
                        check.attributes = FALSE))
  expect_true(all.equal(uc_df,
                        .make_uoa_cluster_df(uoa_spec),
                        check.attributes = FALSE))
  expect_equal(uc_df[,1,drop=FALSE], clusters(cluster_spec))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, letters[1:4],
                        check.attributes = FALSE))
  
  # test when there's a subset
  subset_spec <- rct_spec(a ~ unitid(uid1, uid2), spec_data, subset = bid == "A")
  expect_equal(.make_uoa_cluster_df(subset_spec, c("uid1", "uid2")),
               data.frame(uid1 = letters[1:2], uid2 = letters[1:2],
                          cluster = paste(letters[1:2], letters[1:2], sep = "_")))

  # test cluster arg specification when level is coarser than assignment level
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_spec, "bid")),
               c(4, 2))
  expect_equal(uc_df[,1,drop=FALSE], clusters(cluster_spec))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, rep(LETTERS[1:2], each = 2),
                        check.attributes = FALSE))

  # test cluster arg specification when level is finer than assignment level
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_spec, "id")),
               c(20, 2))
  expect_true(all.equal(uc_df[,1], rep(clusters(cluster_spec)[,1], each = 5),
                        check.attributes = FALSE))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, letters[1:20],
                        check.attributes = FALSE))

  # test cluster arg specification with multiple columns
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_spec, c("uid2", "bid"))),
               c(4, 2))
  expect_true(all.equal(uc_df$cluster,
                        paste(letters[1:4], rep(LETTERS[1:2], each = 2),
                              sep = "_"),
                        check.attributes = FALSE))

  # test specification with multiple unit of assignment columns
  cluster_spec <- rct_spec(a ~ cluster(uid1, uid2), spec_data)
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_spec)),
               c(4, 3))
  expect_equal(uc_df[,1:2,drop=FALSE], clusters(cluster_spec))
  expect_equal(colnames(uc_df), c("uid1", "uid2", "cluster"))
  expect_true(all.equal(uc_df$cluster,
                        paste(letters[1:4], letters[1:4], sep = "_"),
                        check.attributes = FALSE))
})
