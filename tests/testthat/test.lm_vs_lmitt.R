# Ensures for each combination of usages below, lm and lmitt produce the same
# point esimate.

### Dimensions:
# Weights Yes/No
# Subgroup Yes/No
# Absorb Yes/No
# CovAdj Yes/No
# Subset No/In `StudySpecification`/in `lm`

# To avoid a bunch of messages being printed:
save_options <- options()
options("propertee_message_on_unused_blocks" = FALSE)

test_that("equivalent of lm and lmitt calls - no subset", {
  data(simdata)
  simdata$o <- as.factor(simdata$o)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)


  test_coeffs <- function(lmcoef, lmittcoef, n_sg_levels, cov_adj) {
    ncoef <- if (cov_adj) 4 else 3
    if (n_sg_levels > 0 & !cov_adj) {
      ncoef <- ncoef + 3 * (n_sg_levels - 1)
    } else if (n_sg_levels > 0 & cov_adj) {
      ncoef <- ncoef + 4 * (n_sg_levels - 1)
    }
    expect_equal(length(lmittcoef), ncoef)
    lmittcoef <- lmittcoef[!grepl("(Intercept)", names(lmittcoef))]
    expect_true(all(round(lmittcoef[!grepl("^(y|cov_adj):", names(lmittcoef))], 4) %in% round(lmcoef, 4)))
  }

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec), data = simdata, weights = ate(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata, weights = ett(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               weights = ett(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
                     weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o  + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |   No   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               weights = ate(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec), data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(), data = simdata, offset = cov_adj(cmod),
               weights = ate(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec,
                     offset = cov_adj(cmod), weights = ate())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata,
               offset = cov_adj(cmod), weights = ett(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, offset = cov_adj(cmod),
                     weights = ett())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
               offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ett(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
               offset = cov_adj(cmod), weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |  Yes   |   No   |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ate(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     offset = cov_adj(cmod), weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)


})

test_that("subset in specification", {
  data(simdata)
  simdata$o <- as.factor(simdata$o)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata,
                    subset = !(simdata$uoa1 == 5 & simdata$uoa2 == 2))
  cmod <- lm(y ~ x, data = simdata)

  test_coeffs <- function(lmcoef, lmittcoef, n_sg_levels, cov_adj) {
    ncoef <- if (cov_adj) 4 else 3
    if (n_sg_levels > 0 & !cov_adj) {
      ncoef <- ncoef + 3 * (n_sg_levels - 1)
    } else if (n_sg_levels > 0 & cov_adj) {
      ncoef <- ncoef + 4 * (n_sg_levels - 1)
    }
    expect_equal(length(lmittcoef), ncoef)
    lmittcoef <- lmittcoef[!grepl("(Intercept)", names(lmittcoef))]
    expect_true(all(round(lmittcoef[!grepl("^(y|cov_adj):", names(lmittcoef))], 4) %in% round(lmcoef, 4)))
  }



  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec), data = simdata, weights = ate(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata, weights = ett(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               weights = ett(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
                     weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o  + as.factor(bid), data = simdata)
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE)
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |   No   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               weights = ate(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), FALSE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |   No   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec), data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |   No   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(), data = simdata, offset = cov_adj(cmod),
               weights = ate(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec,
                     offset = cov_adj(cmod), weights = ate())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |   No   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata, offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |   No   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o, data = simdata,
               offset = cov_adj(cmod), weights = ett(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, offset = cov_adj(cmod),
                     weights = ett())
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   No     |  Yes   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
               offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   No     |  Yes   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec) + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ett(spec))
  mod_lmitt <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE,
               offset = cov_adj(cmod), weights = "ett")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, 0, TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   No    |   Yes    |  Yes   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     offset = cov_adj(cmod))
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)

  ## | Weights | Subgroup | Absorb | CovAdj | Subset |
  ## |---------|----------|--------|--------|--------|
  ## |   Yes   |   Yes    |  Yes   |  Yes   | StudySpecification |

  mod_lm <- lm(y ~ assigned(spec):o + o + as.factor(bid), data = simdata,
               offset = cov_adj(cmod), weights = ate(spec))
  mod_lmitt <- lmitt(y ~ o, data = simdata, specification = spec, absorb = TRUE,
                     offset = cov_adj(cmod), weights = "ate")
  test_coeffs(mod_lm$coeff, mod_lmitt$coeff, length(levels(simdata$o)), TRUE)
})

options(save_options)
