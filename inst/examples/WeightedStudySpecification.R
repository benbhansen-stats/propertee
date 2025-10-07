data(schooldata)
data(studentdata)
spec <- rct_spec(treatment ~ unitid(schoolid) + block(state), schooldata)
wt <- ate(spec, data = schooldata)
wt.subset <- wt[seq_len(4)]
identical(subset(wt, c(rep(TRUE, 4), rep(FALSE, nrow(schooldata)-4))),
          wt.subset)
show(wt.subset)
