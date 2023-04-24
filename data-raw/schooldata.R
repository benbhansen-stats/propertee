set.seed(3)
nschool <- 139
nstudent <- 8713

schooldata <- data.frame(schoolid = seq_len(nschool),
                         treatment = rbinom(nschool, 1, .2),
                         state = sample(seq_len(5), nschool, replace = TRUE,
                                        prob = c(.2, .25, .1, .35, .2)))

studentdata <- data.frame(id = seq_len(nstudent),
                          schoolid = sample(seq_len(nschool),
                                           size = nstudent,
                                           replace = TRUE,
                                           prob = runif(nschool, .2, .8)),
                          grade = sample(c(3:5), nstudent,
                                          replace = TRUE),
                          gpa = round(4*(1 - rbeta(nstudent, 2, 5)), 1))

merged <- .merge_preserve_order(studentdata, schooldata, by = "schoolid")

merged$math <- 68*(1-rbeta(nstudent, 2, 3)) + 32#runif(nstudent, 23, 32)
merged$math <- merged$math + 2*(merged$grade == 4) +
                .5*merged$gpa +
                2.5*(merged$grade == 5) +
                2.3*(merged$grade == 3 & merged$treatment == 1) +
                3.9*(merged$grade == 4 & merged$treatment == 1) +
                1.4*(merged$grade == 5 & merged$treatment == 1)

merged$math <- round(merged$math)
merged$math <- pmin(merged$math, 100)
merged$math <- pmax(merged$math, 0)

studentdata$math <- merged$math

des <- obs_design(treatment ~ uoa(schoolid), data = schooldata)

mod1 <- lmitt(math ~ as.factor(grade), design = des, data = studentdata)
mod2 <- lmitt(math ~ as.factor(grade), design = des, data = merged)
mod3 <- lm(math ~ as.factor(grade):treatment + as.factor(grade), data = merged)

stopifnot(all(mod1$coef == mod2$coef))
stopifnot(all.equal(mod3$coef[4:6], mod1$coef,
                    check.attributes = FALSE,
                    scale = 1, tol = 1e-10)) # slight rounding discrepancies

usethis::use_data(schooldata, overwrite = TRUE)
usethis::use_data(studentdata, overwrite = TRUE)
