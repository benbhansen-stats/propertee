library(AER, quietly = TRUE)
data(STAR)
STAR$treatment <- STAR$stark == "small"
STAR$treatment[is.na(STAR$treatment)] <- FALSE
STAR$studentid <- as.character(1:nrow(STAR))

covariance_y0_read <- lm(readk ~ gender + ethnicity + birth + lunchk +
                                 ladderk + experiencek + tethnicityk,
                         data = STAR, subset = !treatment)

STAR_design <- RCT_Design(I(1 * treatment) ~ cluster(studentid), data = STAR) 
STAR_ate    <- ate(STAR_design)
STAR_ett    <- ett(STAR_design)

ett_read <- lm(readk ~ treatment, 
               data = STAR,
               weights = cov_adj(STAR_ett, covariance_y0_read))
coef(ett_read)

ate_read <- lm(readk ~ treatment, 
               data = STAR,
               weights = cov_adj(STAR_ate, covariance_y0_read))
coef(ate_read)

ate_read_eth <- lm(readk ~ treatment*ethnicity,
                   data = STAR,
                   weights = cov_adj(STAR_ate, covariance_y0_read))
coef(ate_read_eth)
