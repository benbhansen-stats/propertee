data("STARplus")

##' A prognostic model fitted to experimental + non-experimental controls
y0hat_read <- lm(read_post_yr1 ~ gender*dob + race,
                 data = STARplus,
                 subset = cond_at_entry!="small")

STARspec <- rct_spec(cond_at_entry ~ unit_of_assignment(stdntid) +
                         block(grade_at_entry, school_at_entry),
                     subset=!is.na(grade_at_entry),# excludes non-experimentals
                     data = STARplus)
ett_wts    <- ett(STARspec, data = STARplus,
                  dichotomy= cond_at_entry =="small" ~.)

ett_read <- lm(read_post_yr1 ~ assigned(dichotomy= cond_at_entry =="small" ~.),
               ### expect warning about NAs generated here:
               offset = cov_adj(y0hat_read),
               data = STARplus,
               weights = ett_wts)
coef(ett_read)
ett_read |> as.lmitt() |> vcov()

ate_read <- lmitt(read_post_yr1 ~ 1, STARspec, STARplus,
                  dichotomy= cond_at_entry =="small" ~.,
                  offset = cov_adj(y0hat_read),
                  weights = "ate")
show(ate_read)
coef(ate_read)
vcov(ate_read)

ate_read_loc <-
    lmitt(read_post_yr1 ~ race, STARspec, STARplus,
          dichotomy= cond_at_entry =="small" ~.,
          offset = cov_adj(y0hat_read, newdata = STARplus),
          weights = "ate")
show(ate_read_loc)
