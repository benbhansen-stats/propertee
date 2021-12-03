# Compare lm_robust to flexida results

library(tidyverse)
library(flexida)
library(estimatr)
library(blkvar)
library(here)

test_dat <- read_csv(here("tests/testthat","test_dat.csv"),col_types="iicci")
## This is a cluster-randomized design (cid) within blocks (bid) (actually the
## result of full matching at neighborhood level with outcomes measured using a survey
## within neighborhoods).
head(test_dat)
test_dat$uid <- 1:nrow(test_dat)
with(test_dat,table(z,bid))
with(test_dat,table(z,cid))

des <- RCT_Design(z ~ cluster(cid) + block(bid), data = test_dat)
summary(des)
it_ate <- ittestimate(des, test_dat, "y", "ate")
it_ett <- ittestimate(des, test_dat, "y", "ett")
lm_ate <- lm(y~z,data=test_dat,weights=ate(des,data=test_dat,clusterIds=list("cid"="cid")))
lm_ett <- lm(y~z,data=test_dat,weights=ett(des,data=test_dat,clusterIds=list("cid"="cid")))

lr_cr2 <- lm_robust(y~z,fixed_effects=~bid,clusters=cid,data=test_dat,se_type="CR2")
lr_cr0 <- lm_robust(y~z,fixed_effects=~bid,clusters=cid,data=test_dat,se_type="CR0")
lr_stata <- lm_robust(y~z,fixed_effects=~bid,clusters=cid,data=test_dat,se_type="stata")

results <-  list(
        it_ate=filter(tidy(it_ate),term=="Design_Treatment1"),
        it_ett=filter(tidy(it_ett),term=="Design_Treatment1"),
        lm_ate=filter(tidy(lm_ate),term=="z"),
        lm_ett=filter(tidy(lm_ett),term=="z"),
        lr_cr2=filter(tidy(lr_cr2),term=="z"),
        lr_cr0=filter(tidy(lr_cr0),term=="z"),
        lr_stata=filter(tidy(lr_stata),term=="z")
)

compare_results <- bind_rows(results)
compare_results$version <- names(results)

## Large differences across lm_robust and the two different flexida approaches
compare_results


# Step back. Just pretend that this is block-randomized.

test_dat <- test_dat %>% group_by(cid) %>% mutate(n_clus=n(),n_trt_clus=sum(z),
n_ctrl_clus=sum(1-z))

test_dat <- test_dat %>% group_by(bid) %>% mutate(n_indiv_block=n(),
    p_z_indiv_block=mean(z),
    mean_clus_size_block=mean(n_clus),
    nbwt = (z / p_z_indiv_block) + ((1 - z) / (1 - p_z_indiv_block)),
    hbwt = nbwt * (p_z_indiv_block * (1 - p_z_indiv_block))
    ) %>% ungroup()

## We know that difference_in_means does the block-size weighting
tidy(difference_in_means(y~z,blocks=bid,data=test_dat))
## Verify that this does block-size weighting
tidy(lm(y~z,data=test_dat,weights=nbwt))
tidy(lm_robust(y~z,data=test_dat,weights=nbwt,se_type="HC2"))

## This next should be precision weighting
tidy(lm_robust(y~z,fixed_effects=~bid,data=test_dat))
tidy(lm_robust(y~z+bid,data=test_dat)) %>% filter(term=="z")
tidy(lm_robust(y~z,data=test_dat,weights=hbwt))

## Now trying flexida
des_blocked <- RCT_Design(z~block(bid)+unitid(uid),data=test_dat)
tidy(lm(y~z,data=test_dat,weights=ate(des_blocked,data=test_dat)))
tidy(lm(y~z,data=test_dat,weights=ett(des_blocked,data=test_dat)))

## Using the blkvar package to see whether ate() and ett() are doing something
## in there.
compare_methods(Yobs=test_dat$y,Z=test_dat$z,B=test_dat$bid)


## Stopping here because, as should have
## been clear from the code, ate() and ett() just return indices and not
## weights!


