library(robustbase)
data(studentdata)
robreg <- lmrob(math ~ gpa + factor(grade), studentdata)
bread(robreg)
head(estfun(robreg))
studentdata$prof <- studentdata$math > 85
robreg <- glmrob(prof ~ grade + gpa, studentdata, family = stats::binomial())
bread(robreg)
head(estfun(robreg))