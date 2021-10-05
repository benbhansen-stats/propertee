### this isn't intended as part of the package, but will create the fake data from the package
### that means I get to use tidyverse :)
library(tidyverse)
library(synthpop)

set.seed(613)

## first, load in the data
dat <- foreign::read.dta('data_for_analysis.dta')%>%
    mutate(
        R=round(dist_from_cut,2),
        hsgrade_pct=ifelse(hsgrade_pct==100,99.5,hsgrade_pct),
        lhsgrade_pct=qlogis(hsgrade_pct/100))%>%
    select(
        R,lhsgrade_pct,nextGPA,probation_year1,
        lhsgrade_pct,totcredits_year1,male,loc_campus1,loc_campus2,
        bpl_north_america,english,age_at_entry
    )%>%
    na.omit()

### plot the original RD
dat%>%
    group_by(R)%>%
    summarize(nextGPA=mean(nextGPA),n=n())%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,nextGPA,size=n,color=A))+
    geom_point()+geom_smooth()

### remove the treatment effect estimated in S&H (2020)
### we'll add it back in later
dat <- mutate(dat,nextGPA=ifelse(probation_year1==1,nextGPA-0.24,nextGPA))

## now let's look again:
dat%>%
    group_by(R)%>%
    summarize(nextGPA=mean(nextGPA),n=n())%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,nextGPA,size=n,color=A))+
    geom_point()+geom_smooth()


### look at # at each value of R
dat%>%
    group_by(R)%>%
    summarize(logn=log(n()))%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,logn,color=A))+
    geom_point()+geom_smooth()




## also plot lhsgrade_pct for good measure:
dat%>%
    group_by(R)%>%
    summarize(lhsgrade_pct=mean(lhsgrade_pct),n=n())%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,lhsgrade_pct,size=n,color=A))+
    geom_point()+geom_smooth()


### create synthetic data
faux <- syn(dat,minnumlevels=6)

### built-in comparison tool
compare(faux,dat)

### the actual fake data:
fdat <- faux$syn%>%
    mutate_if(is.factor,~as.numeric(as.character(.)))

### now plot the new data, side by side with old:
## nextGPA
bind_rows(
    mutate(dat,type='real'),mutate(fdat,type='fake'))%>%
    group_by(R,type)%>%
    summarize(nextGPA=mean(nextGPA),n=n())%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,nextGPA,size=n,color=A))+
    geom_point()+geom_smooth()+facet_wrap(~type,ncol=2)

## mccrary data
bind_rows(
    mutate(dat,type='real'),mutate(fdat,type='fake'))%>%
    group_by(R,type)%>%
    summarize(n=n(),logn=log(n()))%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,n,color=A))+
    geom_point()+geom_smooth()+facet_wrap(~type,ncol=2)#+xlim(-0.1,0.1)

## lhsgrade_pct

bind_rows(
    mutate(dat,type='real'),mutate(fdat,type='fake'))%>%
    group_by(R,type)%>%
    summarize(lhsgrade_pct=mean(lhsgrade_pct),n=n())%>%
    mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
    ggplot(aes(R,lhsgrade_pct,size=n,color=A))+
    geom_point()+geom_smooth()+facet_wrap(~type,ncol=2)


#### now check/fix mccrary test violation
par(mfrow=c(1,2))
rdd::DCdensity(dat$R,-0.005, bin=0.01,verbose=TRUE)
rdd::DCdensity(fdat$R,-0.005, bin=0.01,verbose=TRUE)

### now for donut design
dat <- subset(dat,R!=0)
fdat <- subset(fdat,R!=0)


par(mfrow=c(1,2))
rdd::DCdensity(dat$R,-0.005, bin=0.01,verbose=TRUE)
rdd::DCdensity(fdat$R,-0.005, bin=0.01,verbose=TRUE)
par(mfrow=c(1,1))
## ### hmmmm donut doesn't work quite as well with fake data
## sum(fdat$R==-0.01)
## sum(fdat$R==-0.02)

## fdat <- rbind(
##     subset(fdat,R!=-0.01 & R!=-0.02),
##     subset(fdat,R==-0.01)[sample(1:sum(fdat$R==-0.01),sum(dat$R==-0.01)),],
##     subset(fdat,R==-0.02)[sample(1:sum(fdat$R==-0.02),sum(dat$R==-0.02)),]
## )

## rdd::DCdensity(fdat$R,-0.005, bin=0.01,verbose=TRUE)
## ## that's more like it

## ### Now add the treatment effect back in:
## fdat <- mutate(fdat,nextGPA=ifelse(probation_year1==1,nextGPA+0.24,nextGPA))

## ## and plot
## fdat%>%
##     group_by(R)%>%
##     summarize(nextGPA=mean(nextGPA),n=n())%>%
##     mutate(A=ifelse(R<=0,'AP','non-AP'))%>%
##     ggplot(aes(R,nextGPA,size=n,color=A))+
##     geom_point()+geom_smooth()


save(fdat,file='../data/fakeRDD.RData')
