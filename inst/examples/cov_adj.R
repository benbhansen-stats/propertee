data(STARdata)
STARdata$studentID <- rownames(STARdata)

y0hat_read <- lm(readk ~ gender + ethnicity + lunchk +
                     ladderk + experiencek + tethnicityk,
                 data = STARdata, subset = stark !="small")

STARspec <- rct_spec(stark ~ unit_of_assignment(studentID),
                     data = STARdata)
ett_wts    <- ett(STARspec, data = STARdata,
                  dichotomy= stark =="small" ~.)

ett_read <- lm(readk ~ assigned(dichotomy= stark =="small" ~.),
               ### expect warning about NAs generated here:
               offset = cov_adj(y0hat_read, newdata = STARdata),
               data = STARdata,
               weights = ett_wts)
coef(ett_read)
ett_read |> as.lmitt() |> vcov()

ate_read <- lmitt(readk ~ 1, STARspec, STARdata,
                  dichotomy= stark =="small" ~.,
                  offset = cov_adj(y0hat_read, newdata = STARdata),
                  weights = "ate")
show(ate_read)
coef(ate_read)
vcov(ate_read)

ate_read_loc <-
    lmitt(readk ~ schoolk, STARspec, STARdata,
          dichotomy= stark =="small" ~.,
          offset = cov_adj(y0hat_read, newdata = STARdata),
          weights = "ate")
show(ate_read_loc)
