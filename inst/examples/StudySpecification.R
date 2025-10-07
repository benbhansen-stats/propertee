data(simdata)
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
                 data = simdata)
show(spec)
summary(spec)

data(schooldata)
spec <- obs_spec(treatment ~ unit_of_assignment(schoolid) + block(state),
                 data = schooldata)
show(spec)
summary(spec)