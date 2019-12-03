rm(list=ls())

#load simulated data
load("~/Documents/nBox/git_projects/NItrialsimulation/manuscript_graphs/fullset.Rdata")
fullset=

library(geepack)
pp<- as.data.frame(rbind(d20[d20$r==1 & d20$int==1,],d20[d20$r==0 & d20$int==0,]))# perprotocol population
E.out<-glm(pp$int~pp$ds,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
ps<-predict(E.out, type="response")
sptw<-pp$int*mean(pp$int)/ps+(1-pp$int)*(1-mean(pp$int))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
((geeglm(recurrence~int, family=binomial(link="identity"), weight=sptw, id=pt.id, data=pp)))$coefficients['int']# fit linear binomial model to the weighted data

library(MatchIt); library(Matching)
m<- glm(int~ds, family = binomial(), data=d) 
pscore<- m$fitted.values
mout<-matchit(int~ds, data=d, method='nearest', distance = "mahalanobis", replace = TRUE)
match<-Match(Tr=d$int, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
matched<-d[unlist(match[c('index.treated','index.control')]),]
y_trt<-matched$recurrence[matched$int==1]
y_con<-matched$recurrence[matched$int==0]
diffy<-y_trt-y_con
(t.test(diffy))$estimate

y<- d20$recurrence
x<- d20$int
iv<-d20$r
data<-data.frame(y,x,iv)
asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"iv"], vcov="iid")
(summary(asmm))$ coefficients [2,1]
