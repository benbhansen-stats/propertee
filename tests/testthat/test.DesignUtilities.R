test_that("Design formula checking", {
  expect_true(.check_design_formula(y ~ cluster(x)))
  expect_true(.check_design_formula(y ~ cluster(x, z, q, r)))
  expect_true(.check_design_formula(y ~ unitid(x)))
  expect_true(.check_design_formula(y ~ unitid(x, z, q, r)))
  expect_true(.check_design_formula(y ~ unit_of_assignment(x)))
  expect_true(.check_design_formula(y ~ unit_of_assignment(x, z, q, r)))
  expect_true(.check_design_formula(y ~ uoa(x)))
  expect_true(.check_design_formula(y ~ uoa(x, z, q, r)))

  expect_error(.check_design_formula(~ cluster(x)),
               "treatment")

  expect_error(.check_design_formula(y ~ x),
               "cluster")

  expect_error(.check_design_formula(y ~ cluster(x) + unitid(z)),
               "Only one of")

  expect_error(.check_design_formula(y ~ cluster(x) + cluster(z)),
               "Only one instance of `cluster")

  expect_error(.check_design_formula(y ~ unitid(x) + unitid(z)),
               "Only one instance of `unitid")

  expect_error(.check_design_formula(y ~ unit_of_assignment(x) +
                                       unit_of_assignment(z)),
               "Only one instance of `unit_of")

  expect_error(.check_design_formula(y ~ uoa(x) + uoa(z)),
               "Only one instance of `unit_of")

  expect_true(.check_design_formula(y ~ cluster(x) + block(z)))
  expect_true(.check_design_formula(y ~ cluster(x) + block(z, a, b, c)))

  expect_error(.check_design_formula(y ~ cluster(x) + block(z) + block(q)),
               "only one block")

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z)),
               "only allowed")

  expect_true(.check_design_formula(y ~ cluster(x) + forcing(z),
                                 allow_forcing = TRUE))

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z) + forcing(q),
                                  allow_forcing = TRUE),
               "only one forcing")

})

test_that("binary treatment and dichotomy", {
  expect_error(is_dichotomized(1),
               "must be")
  expect_error(has_binary_treatment(1),
               "must be")
  expect_error(is_binary_or_dichotomized(1),
               "must be")

  des1 <- obs_design(z ~ uoa(cid1, cid2), data = simdata)
  des2 <- obs_design(o ~ uoa(cid1, cid2), data = simdata)
  des3 <- obs_design(o ~ uoa(cid1, cid2), data = simdata,
                     dichotomy = o > 3 ~ o == 1)
  des4 <- obs_design(z ~ uoa(cid1, cid2), data = simdata,
                     dichotomy = z == 1 ~ z == 0)

  expect_false(is_dichotomized(des1))
  expect_true(has_binary_treatment(des1))
  expect_true(is_binary_or_dichotomized(des1))

  expect_false(is_dichotomized(des2))
  expect_false(has_binary_treatment(des2))
  expect_false(is_binary_or_dichotomized(des2))

  expect_true(is_dichotomized(des3))
  expect_false(has_binary_treatment(des3))
  expect_true(is_binary_or_dichotomized(des3))

  expect_true(is_dichotomized(des4))
  expect_true(has_binary_treatment(des4))
  expect_true(is_binary_or_dichotomized(des4))

  # Adding afterwards
  dichotomy(des2) <- o >= 2 ~ .
  expect_true(is_dichotomized(des2))
  expect_false(has_binary_treatment(des2))
  expect_true(is_binary_or_dichotomized(des2))

  # Removing
  dichotomy(des3) <- NULL
  expect_false(is_dichotomized(des3))
  expect_false(has_binary_treatment(des3))
  expect_false(is_binary_or_dichotomized(des3))
})

test_that("identical_Designs function", {
  data(simdata)

  des1 <- rct_design(dose ~ cluster(cid1, cid2), data = simdata)
  des2 <- rct_design(dose ~ cluster(cid1, cid2), data = simdata,
                     dichotomy = dose > 200 ~ dose <= 200)
  des3 <- rct_design(dose ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_true(identical_Designs(des1, des2))
  expect_false(identical_Designs(des1, des3))
  expect_false(identical_Designs(des2, des3))

  expect_false(identical_Designs(des1, des2, TRUE))
  expect_false(identical_Designs(des1, des3, TRUE))
  expect_false(identical_Designs(des2, des3, TRUE))

  # Was hitting an error if the Design passed in via `design=` contained a
  # dichotomy. `des2` above has a dichotomy.
  expect_no_error(lmitt(y ~ 1, design = des2, data = simdata, weights = "ate"))

})

test_that("identify_small_blocks", {
  des_data <- data.frame(uid = letters[1:12],
                         bid = rep(LETTERS[1:3], each = 4),
                         a = c(rep(rep(c(0, 1), each = 2), 2), 0, rep(1, 3)))
  des <- rct_design(a ~ unitid(uid) + block(bid), des_data)

  small_blocks <- identify_small_blocks(des)
  expect_equal(length(small_blocks), 3)
  expect_equal(names(small_blocks), LETTERS[1:3])
  expect_true(all.equal(small_blocks, c(FALSE, FALSE, TRUE), check.attributes = FALSE))
})

test_that(".make_uoa_cluster_df errors", {
  expect_error(.make_uoa_cluster_df("des"), "valid `Design` object")

  data(simdata)
  simdata_copy <- simdata
  des <- rct_design(z ~ cluster(cid1, cid2), simdata_copy)

  simdata_copy$cid2[simdata_copy$cid1 == 1 & simdata_copy$cid2 == 2] <- 1
  expect_warning(.make_uoa_cluster_df(des), "Some units of assignment")

  simdata_copy <- NULL
  expect_error(.make_uoa_cluster_df(des), "design data")

  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  expect_error(.make_uoa_cluster_df(des, "id"), "id column in the design data")
})

test_that(".make_uoa_cluster_df", {
  des_data <- data.frame(id = letters[1:20],
                         uid1 = rep(letters[1:4], each = 5),
                         uid2 = rep(letters[1:4], each = 5),
                         bid = rep(LETTERS[1:2], each = 10),
                         a = rep(rep(c(0, 1), each = 5), 2))

  # test design creation using `cluster()`, `unitid()`, and `unit_of_assignment()`
  cluster_des <- rct_design(a ~ cluster(uid1), des_data)
  unitid_des <- rct_design(a ~ unitid(uid1), des_data)
  uoa_des <- rct_design(a ~ unit_of_assignment(uid1), des_data)

  expect_true(all.equal(uc_df <- .make_uoa_cluster_df(cluster_des),
                        .make_uoa_cluster_df(unitid_des),
                        check.attributes = FALSE))
  expect_true(all.equal(uc_df,
                        .make_uoa_cluster_df(uoa_des),
                        check.attributes = FALSE))
  expect_equal(uc_df[,1,drop=FALSE], clusters(cluster_des))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, letters[1:4],
                        check.attributes = FALSE))

  # test cluster arg specification when level is coarser than assignment level
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_des, "bid")),
               c(4, 2))
  expect_equal(uc_df[,1,drop=FALSE], clusters(cluster_des))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, rep(LETTERS[1:2], each = 2),
                        check.attributes = FALSE))

  # test cluster arg specification when level is finer than assignment level
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_des, "id")),
               c(20, 2))
  expect_true(all.equal(uc_df[,1], rep(clusters(cluster_des)[,1], each = 5),
                        check.attributes = FALSE))
  expect_equal(colnames(uc_df), c("uid1", "cluster"))
  expect_true(all.equal(uc_df$cluster, letters[1:20],
                        check.attributes = FALSE))

  # test cluster arg specification with multiple columns
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_des, c("uid2", "bid"))),
               c(4, 2))
  expect_true(all.equal(uc_df$cluster,
                        paste(letters[1:4], rep(LETTERS[1:2], each = 2),
                              sep = "_"),
                        check.attributes = FALSE))

  # test design with multiple unit of assignment columns
  cluster_des <- rct_design(a ~ cluster(uid1, uid2), des_data)
  expect_equal(dim(uc_df <- .make_uoa_cluster_df(cluster_des)),
               c(4, 3))
  expect_equal(uc_df[,1:2,drop=FALSE], clusters(cluster_des))
  expect_equal(colnames(uc_df), c("uid1", "uid2", "cluster"))
  expect_true(all.equal(uc_df$cluster,
                        paste(letters[1:4], letters[1:4], sep = "_"),
                        check.attributes = FALSE))
})
