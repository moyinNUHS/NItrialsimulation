rm(list=ls())

r<- c(rep(1,10),rep(0,10))
pt.id<-1:20
ds<-rep(1:10,2)
phy.id<-c(2,2,1,1,1,2,2,1,1,2,2,1,2,2,1,1,1,2,2,2)
int<-c(1,1,1,1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0)
recurrence<-c(0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1)
compliance<-c(1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1)
d<-data.frame(cbind(pt.id, r, ds, phy.id, int, recurrence))

t.test((d[d$int==1,])$ds,(d[d$int==0,])$ds)
t.test((d[d$recurrence==1,])$ds,(d[d$recurrence==0,])$ds)

(a<-nrow(d[phy.id==1 & int==1,]))
(b<-nrow(d[phy.id==1 & int==0,]))
(c<-nrow(d[phy.id==2 & int==1,]))
(j<-nrow(d[phy.id==2 & int==0,]))
chisq.test(x=matrix(c(a,c,b,j),2,2), simulate.p.value = TRUE, B=1000)

(e<-nrow(d[phy.id==1 & recurrence==1,]))
(f<-nrow(d[phy.id==1 & recurrence==0,]))
(g<-nrow(d[phy.id==2 & recurrence==1,]))
(h<-nrow(d[phy.id==2 & recurrence==0,]))
chisq.test(x=matrix(c(e,g,f,h),2,2), simulate.p.value = TRUE, B=1000)

#repeat each row 20 times 
d20<-d[rep(seq_len(nrow(d)), each=20),]
nrow(d20)

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
