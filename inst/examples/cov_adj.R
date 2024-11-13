data(STARdata)
STARdata$treatment <- STARdata$stark == "small"
STARdata$treatment[is.na(STARdata$treatment)] <- FALSE
STARdata$studentid <- as.character(seq_len(nrow(STARdata)))

covariance_y0_read <- lm(readk ~ gender + ethnicity + lunchk +
                                 ladderk + experiencek + tethnicityk,
                         data = STARdata, subset = !treatment)

STARdata_spec <- rct_spec(treatment ~ cluster(studentid), data = STARdata)
STARdata_ate    <- ate(STARdata_spec, data = STARdata)
STARdata_ett    <- ett(STARdata_spec, data = STARdata)

ett_read <- lm(readk ~ treatment,
               offset = cov_adj(covariance_y0_read, newdata = STARdata),
               data = STARdata,
               weights = STARdata_ett)
coef(ett_read)

ate_read <- lm(readk ~ treatment,
               offset = cov_adj(covariance_y0_read, newdata = STARdata),
               data = STARdata,
               weights = STARdata_ate)
coef(ate_read)

ate_read_eth <- lm(readk ~ treatment * ethnicity,
                   offset = cov_adj(covariance_y0_read, newdata = STARdata),
                   data = STARdata,
                   weights = STARdata_ate)
coef(ate_read_eth)
