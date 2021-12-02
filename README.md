# flexida

[![R-build-check](https://github.com/benbhansen-stats/flexida/actions/workflows/r.yml/badge.svg)](https://github.com/benbhansen-stats/flexida/actions/workflows/r.yml)

## Overview

Flexida enables flexible direct adjustment with design-informed standard errors
and optional prior covariance adjustment.

Random trials often utilize clustering and blocking in assigning treatment
status as a way to simplify implementation. This design information must be
utilized in future analyses. Using Flexida, a user can generate a Design object
which will keep track of the design structure.

    des <- RCT_Design(txt ~ cluster(teacher) + block(school), data = teacherdata)

(Also supported are observational studies (`Obs_Design`) and regression
discontinuity designs (`RDD_Design` which requires a `forcing()` variable as
well.)

In order to pass the design information into the model using the weights=
argument, functions ett() and ate() will be used to convert the Design into a
numeric vector with the Design object as an attribute.

    lm(y ~ x, data = studentdata, weights = ate(des))

Note that the Design is created with teacher level data (`teacherdata`), but the
analysis is carried out at the student level (`studentdata`); the `ate()` (and
its alternative `ett()`) will expand the weights appropriately.

Alternatively, the `ittestimate` function offers a more customizeable treatment
effect estimation, including the use of a separate covariate adjustment model.

    covadjmod <- lm(y ~ x1 + x2 + ..., data = studentdata)
    ittest <- ittestimate(des, studentdata,
                          outcome = "y",
                          target = "ett",
                          covAdjModel = covadjmod)
