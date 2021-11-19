set.seed(2)
simdata <- data.frame(cid1 = rep(1:5, each = 10),
                      cid2 = rep(rep(1:2, times = c(4, 6)),
                                 times = 5),
                      bid = rep(1:3, times = c(20, 14, 16)),
                      force = rep(rnorm(10, mean = 5),
                                  times = rep(c(4, 6),
                                              times = 5)),
                      z = rep(sample(0:1, 10, TRUE),
                              times = rep(c(4, 6),
                                          times = 5)),
                      o = as.ordered(rep(sample(1:4, 10, TRUE),
                                         times = rep(c(4, 6),
                                                     times = 5))),
                      dose = rep(rep(c(50, 250, 100, 200, 300), 2),
                                 times = rep(c(4, 6),
                                             times = 5)),
                      x = rnorm(1:50),
                      y = rnorm(1:50))

usethis::use_data(simdata, overwrite = TRUE)
