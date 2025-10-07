data(schooldata)
data(studentdata)
covadj_mod <- lm(math ~ gpa + grade, studentdata[seq_len(15),])
suppressWarnings(psl <- cov_adj(covadj_mod))
show(psl)
identical(subset(psl, c(rep(TRUE, 4), rep(FALSE, 11))),
          psl[seq_len(4)])