data(schooldata)
data(studentdata)
spec <- rct_spec(treatment ~ unitid(schoolid) + block(state), schooldata)
wt <- ate(spec, data = schooldata)
wt.subset <- wt[seq_len(4)]
identical(subset(wt, c(rep(TRUE, 4), rep(FALSE, nrow(schooldata)-4))),
          wt.subset)
show(wt.subset)
weights(wt.subset)
wt2 <- ett(spec, data = schooldata)
wt2.subset <- wt2[seq_len(4)]
can_add <- tryCatch(wt.subset + wt2.subset, error = function(e) FALSE)
can_add
can_multiply <- tryCatch(wt.subset * wt2.subset, error = function(e) FALSE)
can_multiply
can_divide <- tryCatch(wt.subset / wt2.subset, error = function(e) FALSE)
can_divide
other.wt <- c(1, 1, 0, 1)
can_multiply2 <- tryCatch(wt.subset * other.wt, error = function(e) FALSE)
can_multiply2
can_divide2 <- tryCatch(other.wt / wt.subset, error = function(e) FALSE)
can_divide2
