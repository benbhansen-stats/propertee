library(AER, quietly = TRUE)
data(STAR)
STAR$treatment <- STAR$stark == "small"
STAR$treatment[is.na(STAR$treatment)] <- FALSE
STAR$studentid <- as.character(1:nrow(STAR))

covariance_y0_read <- lm(readk ~ gender + ethnicity + birth + lunchk +
                                 ladderk + experiencek + tethnicityk,
                         data = STAR, subset = !treatment)

STAR_design <- RCT_Design(treatment ~ cluster(studentid), data = STAR)
STAR_ate    <- ate(STAR_design)
STAR_ett    <- ett(STAR_design)

ett_read <- lm(readk ~ treatment,
               offset = cov_adj(covariance_y0_read),
               data = STAR,
               weights = STAR_ett)
coef(ett_read)

ate_read <- lm(readk ~ treatment,
               offset = cov_adj(covariance_y0_read),
               data = STAR,
               weights = STAR_ate)
coef(ate_read)

ate_read_eth <- lm(readk ~ treatment*ethnicity,
                   offset = cov_adj(covariance_y0_read),
                   data = STAR,
                   weights = STAR_ate)
coef(ate_read_eth)
