test_that("equivalent of lm and lmitt calls", {
  data(simdata)
  simdata$o <- as.factor(simdata$o)
  des <- obs_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)

  ### Dimensions:
  # Weights Yes/No
  # Subgroup Yes/No
  # Absorb Yes/No
  # CovAdj Yes/No
  # Subset No/In `Design`/in `lm`

  test_coeffs <- function(lmcoef, lmittcoef, subgroup) {
    length <- ifelse(subgroup, 4, 2)
    expect_equal(length(lmittcoef), length)
    lmittcoef <- lmittcoef[!grepl("(Intercept)", names(lmittcoef))]
    expect_true(all(round(lmittcoef, 4) %in% round(lmcoef, 4)))
  }

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des), data = simdata, weights = ate(des))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o, data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o, data = simdata, weights = ett(des))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata,
               weights = ett(des))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE,
                     weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o  + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o + as.factor(bid), data = simdata,
               weights = ate(des))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE,
                     weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des), data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(), data = simdata, offset = cov_adj(cmod),
               weights = ate(des))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des,
                     offset = cov_adj(cmod), weights = ate())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o, data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o, data = simdata,
               offset = cov_adj(cmod), weights = ett(des))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, offset = cov_adj(cmod),
                     weights = ett())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE,
               offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ett(des))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE,
               offset = cov_adj(cmod), weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(des):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ate(des))
  mod_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE,
                     offset = cov_adj(cmod), weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, TRUE)


})
