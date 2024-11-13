nschool <- 139
nstudent <- 8713

schoolids <- sapply(seq_len(nschool), function(x) {
  set.seed(x)
  paste(LETTERS[sample(1:26, 10, replace = TRUE)], collapse = "")
})
stopifnot(length(unique(schoolids)) == nschool)

set.seed(5)

schooldata <- data.frame(schoolid = schoolids,
                         treatment = rbinom(nschool, 1, .2),
                         state = sample(seq_len(5), nschool, replace = TRUE,
                                        prob = c(.2, .25, .1, .35, .2)),
                         pct_disadvantage = round(runif(nschool, .01, .99), 2))

studentdata <- data.frame(studentid = seq_len(nstudent),
                          schoolid = sample(schoolids,
                                           size = nstudent,
                                           replace = TRUE,
                                           prob = runif(nschool, .2, .8)),
                          grade = sample(c(3:5), nstudent,
                                          replace = TRUE),
                          gpa = round(4*(1 - rbeta(nstudent, 2, 5)), 1))

merged <- merge(schooldata, studentdata, by = "schoolid")

merged$math <- 68*(1-rbeta(nstudent, 2, 3)) + 32#runif(nstudent, 23, 32)
merged$math <- merged$math + 2*(merged$grade == 4) +
                .5*merged$gpa +
                2.5*(merged$grade == 5) +
                2.3*(merged$grade == 3 & merged$treatment == 1) +
                3.9*(merged$grade == 4 & merged$treatment == 1)

merged$math <- round(merged$math)
merged$math <- pmin(merged$math, 100)
merged$math <- pmax(merged$math, 0)

studentdata <- merge(studentdata, merged[, c("studentid", "math")], by = "studentid")

spec <- obs_spec(treatment ~ uoa(schoolid), data = schooldata)

mod1 <- lmitt(math ~ as.factor(grade), specification = spec, data = studentdata)
mod2 <- lmitt(math ~ as.factor(grade), specification = spec, data = merged)
mod3 <- lm(math ~ as.factor(grade):treatment + as.factor(grade), data = merged)

stopifnot(all.equal(mod1$coef, mod2$coef,
                    check.attributes = FALSE,
                    scale = 1, tol = 1e-10))
stopifnot(all.equal(mod1$coef, mod3$coef[4:6],
                    check.attributes = FALSE,
                    scale = 1, tol = 1e-10))

usethis::use_data(schooldata, studentdata, overwrite = TRUE)
