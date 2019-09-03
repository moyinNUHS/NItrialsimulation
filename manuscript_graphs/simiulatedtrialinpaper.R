simdata<-data.frame(
  id=1:20,
  rand=c(rep(1, 10), rep(0, 10)), 
  severe=      c(1,2,3,4,5,6,7,8,9,10, #exp
                 1,2,3,4,5,6,7,8,9,10), #control
  physician=   c(0,0,0,0,1,0,0,1,1,1, #exp
                 0,0,0,0,1,1,1,1,1,1),#control
  intervention=c(1,1,1,1,1,1,1,0,0,0,#exp
                 1,0,0,1,0,0,0,0,0,0),#control
  death=       c(0,0,0,0,1,0,1,1,1,1,#exp
                 0,0,0,0,0,1,1,1,1,1)#control
)
simdata$id=as.factor(simdata$id)
simdata$rand=as.factor(simdata$rand)
simdata$severe=as.numeric(simdata$severe)
simdata$physician=ifelse(simdata$physician==1, 'Good','Bad')
simdata$intervention=as.factor(simdata$intervention)
simdata$death=as.numeric(simdata$death)

## intention to treat 
pz1.value = mean(simdata[which(simdata$rand==1),]$death)
pz0.value = mean(simdata[which(simdata$rand==0),]$death)
eff.itt = pz1.value-pz0.value

## per protocol 
pp = simdata[which(simdata$intervention==simdata$rand),] # perprotocol population
p.experiment.vector= pp[which(pp$rand==1),]$death
p.experiment.value= mean(p.experiment.vector)
p.stdcare.vector= pp[which(pp$rand==0),]$death
p.stdcare.value= mean(p.stdcare.vector) 
eff.pp = p.experiment.value-p.stdcare.value

## inverse probability weights on per protocol patients 
pp=as.data.frame(pp)
colnames(pp)=c('id','randomisation','confounder','phy','intervention','outcome')
ipwmodel=glm(intervention~ confounder+phy,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
score=predict(ipwmodel, type="response")
weight= ifelse(pp$intervention==1,1/score, 1/(1-score)) #create weights
outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
eff.mpp=coef(outcomemodel)[2]


# iv with 2 stage regression
asmm <- gmm(simdata$death ~ simdata$intervention, x=simdata$rand, vcov="iid")
eff.iv=(summary(asmm))$coefficients [2,1]

eff.itt
eff.pp
eff.mpp
eff.iv

