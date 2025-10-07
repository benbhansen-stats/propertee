data(schooldata)
data(studentdata)
spec <- rct_spec(treatment ~ unitid(schoolid), schooldata)
tm <- lmitt(math ~ 1, spec, studentdata, weights = "ate")
show(tm)
summary(tm)

# with covariance adjustment
covadj_mod <- lm(math ~ gpa, studentdata)
tm <- lmitt(math ~ 1, spec, studentdata, weights = "ate",
            offset = cov_adj(covadj_mod))
show(tm)
summary(tm)