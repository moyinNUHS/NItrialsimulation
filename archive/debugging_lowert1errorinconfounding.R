#######################################
########works for non confounding######
#######################################

n=505
p.experiment=0.4
p.stdcare=0.3
p.alt = 0.5
nIterations = 1000
NImargin = 0.1
interval  = 1
comply.experiment=1
comply.stdcare=1

#significance 
 z <- qnorm(0.975)
 
 getoutcome <- function(outcome0, outcome1, outcome2, intervention){
   outcome = c()
   outcome.matrix = matrix(c(outcome0, outcome1, outcome2), ncol = 3)
   for (i in 1:length(intervention)){
     outcome[i] = outcome.matrix[i, (intervention[i]+1)]
   }
   return(unlist(outcome))
 }

 estimates<-sd<-c()
 for (i in 1:1000) {
   id=seq(1,(2*n), by=1) #create participant id  
   randomisation=c(rep(1,n), rep(0,n)) #randomisation
   confounder=rbeta(n=(2*n),shape1=2,shape2=2) #confounder beta distribution ranging 0-1 
   
   #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
   outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))  #probability of outcome if intervention = 0
   outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment,1-p.experiment))  #probability of outcome if intervention = 1
   outcome2 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.alt, 1-p.alt)) #probability of outcome if intervention = alternate to both experimental and stdcare
   
   #INTERVENTION dependent on adherence

  experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(0, 1))
  stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1, 0))
  
   intervention = c(experiment.intervention,stdcare.intervention)
   
   #ACTUAL OUTCOMES depend on intervention
   outcome = getoutcome(outcome0=outcome0, outcome1=outcome1, outcome2 = outcome2, intervention=intervention)
   
   simdata=matrix(data=c(id, randomisation, confounder, intervention, outcome), nrow=(2*n))
   
   ## intention to treat 
   pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
   pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
   eff.itt = pz1.value-pz0.value
   sd.itt = sqrt(pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n)
   
   ## per protocol 
   pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
   p.experiment.vector= pp[which(pp[,2]==1), ][,5]
   p.experiment.value= mean(p.experiment.vector)
   p.stdcare.vector= pp[which(pp[,2]==0),][,5]
   p.stdcare.value= mean(p.stdcare.vector) 
   eff.pp = p.experiment.value-p.stdcare.value
   sd.pp = sqrt(p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector))
   
   
   ## inverse probability weights on per protocol patients 
   pp=as.data.frame(pp)
   colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
   ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
   score=predict(ipwmodel, type="response")
   weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
   outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
   eff.mpp=coef(outcomemodel)[2]
   sd.mpp=summary(outcomemodel)$coefficients[2,2]
   
   # iv with 2 stage regression
   asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
   eff.iv = summary(asmm)$coefficients[2,1] 
   sd.iv=summary(asmm)$coefficients [2,2]
   
   estimates[[i]] = (c(eff.iv, eff.itt, eff.mpp, eff.pp))
   sd[[i]] = c(sd.iv, sd.itt, sd.mpp, sd.pp)
 }
 est.dfnc=matrix(unlist(estimates), byrow=T,ncol = 4)
 sd.dfnc=matrix(unlist(sd), byrow = T, ncol = 4)
 ub.dfnc = est.df + z*sd.df
 mean(ub.dfnc[,1])
 z.iv = quantile(ub.dfnc[,1], probs = 0.025)
 z.itt = quantile(ub.dfnc[,2], probs = 0.025)
 z.mpp = quantile(ub.dfnc[,3], probs = 0.025)
 z.pp = quantile(ub.dfnc[,4], probs = 0.025)

 getoutcome<-function(vector.outcome1, vector.outcome0, intervention){
   outcome<-c()
   for (i in 1:length(vector.outcome0)){
     if (intervention[i]==1) {outcome[[i]]=vector.outcome1[i]} else {outcome[[i]]=vector.outcome0[i]}
   }
   return(unlist(outcome))
 }
 estimates<-sd<-c()
 p.stdcare=0.4
#simulate and derive treatment effect 
for(l in 1:nIterations) { 
    
    id=seq(1,(2*n), by=1) #create participant id  
    randomisation=c(rep(1,n), rep(0,n)) #randomisation
    confounder = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment,1-p.experiment))  #probability of outcome if intervention = 1
    outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))       #probability of outcome if intervention = 0
    
    #INTERVENTION dependent on compliance 
    experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1-comply.experiment,comply.experiment))
    stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(comply.stdcare,1-comply.stdcare))
    intervention = c(experiment.intervention,stdcare.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome1, outcome0, intervention)
    
    simdata = matrix(data=c(id,randomisation,confounder,intervention,outcome), nrow=(2*n))
    
    ## intention to treat 
    pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
    pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
    eff.itt = pz1.value-pz0.value
    sd.itt = sqrt(pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n)

    ## per protocol 
    pp = simdata[which(simdata[,2]==simdata[,4]),] # perprotocol population
    p.experiment.vector = pp[which(pp[,2]==1),][,5]
    p.experiment.value = mean(p.experiment.vector)
    p.stdcare.vector = pp[which(pp[,2]==0),][,5]
    p.stdcare.value = mean(p.stdcare.vector) 
    eff.pp = p.experiment.value-p.stdcare.value
    sd.pp = sqrt(p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector))
    
    ## inverse probability weights on per protocol patients 
    pp = as.data.frame(pp)
    colnames(pp) = c('id','randomisation','confounder','intervention','outcome')
    ipwmodel = glm(intervention~confounder,family=binomial(link="logit"), data = pp) #calculate denominators used in inverse probability weights
    score = predict(ipwmodel, type="response")
    weight = pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
    outcomemodel = glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
    eff.mpp = coef(outcomemodel)[2]
    sd.mpp = sqrt(diag(vcovHC(outcomemodel)))[2]

    # iv with 2 stage regression
    asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
    eff.iv = summary(asmm)$coefficients[2,1]
    sd.iv = summary(asmm)$coefficients [2,2]
    
    estimates[[l]] = (c(eff.iv, eff.itt, eff.mpp, eff.pp))
    sd[[l]] = c(sd.iv, sd.itt, sd.mpp, sd.pp)
  }

  est.df.altnc=matrix(unlist(estimates), byrow=T,ncol = 4)
  sd.df.altnc=matrix(unlist(sd), byrow = T, ncol = 4)
  ub.df.altnc = est.df.altnc + z*sd.df.altnc
  
  mean(ub.df.altnc[,1]<z.iv)
  
#######################################
#######################################
#######################################
#######################################
#######################################
##########works for confounding########
#######################################
n=505
p.experiment=0.4
p.stdcare=0.3
p.alt = 0.5
nIterations = 1000
NImargin = 0.1
interval  = 1
 
#significance 
z <- qnorm(0.975)

getoutcome<-function(vector.outcome1, vector.outcome0, intervention){
  outcome<-c()
  for (i in 1:length(vector.outcome0)){
    if (intervention[i]==1) {outcome[[i]]=vector.outcome1[i]} else {outcome[[i]]=vector.outcome0[i]}
  }
  return(unlist(outcome))
}

estimates<-sd<-c()
for (i in 1:1000) {
  id = seq(1,(2*n), by=1) #create participant id  
  randomisation = c(rep(1,n), rep(0,n)) #randomisation
  confounder = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
  
  #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
  
    #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
    shape2 = runif(1)
    shape1 = shape2*p.experiment/(1-p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    outcome1 = rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome 
    
    shape1 = shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    outcome0 = rbinom(2*n, 1, prob = p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 

  d = matrix(data=c(id, randomisation, confounder), nrow=(2*n))
  d.ordered = matrix(cbind(d[order(d[,3]),], outcome0, outcome1),ncol=5)#order confounder in ascending order 
  d.grouped = rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),])
  
  #INTERVENTION dependent on randomisation and confounders
    shape2 = runif(1, min=2, max=10)
    shape1 = shape2*comply.experiment/(1-comply.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention = rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1 = shape2*(1-comply.stdcare)/(1-(1-comply.stdcare)) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention = rbinom(n, 1, prob=comply.stdcare.ind)
    
  intervention<-c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome = getoutcome(vector.outcome1=d.grouped[,5], vector.outcome0=d.grouped[,4], intervention)
  
  simdata = matrix(cbind(d.grouped[,c(-4,-5)],intervention,outcome), ncol=5)
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
  eff.itt = pz1.value-pz0.value
  sd.itt = sqrt(pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n)
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
  p.experiment.vector= pp[which(pp[,2]==1), ][,5]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,5]
  p.stdcare.value= mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  sd.pp = sqrt(p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector))
  
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
  ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  sd.mpp=summary(outcomemodel)$coefficients[2,2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
  eff.iv = summary(asmm)$coefficients[2,1] 
  sd.iv=summary(asmm)$coefficients [2,2]
  
  estimates[[i]] = (c(eff.iv, eff.itt, eff.mpp, eff.pp))
  sd[[i]] = c(sd.iv, sd.itt, sd.mpp, sd.pp)
}
est.dfc=matrix(unlist(estimates), byrow=T,ncol = 4)
sd.dfc=matrix(unlist(sd), byrow = T, ncol = 4)
ub.dfc = est.df + z*sd.df
mean(ub.df[,1])
z.iv = quantile(ub.df[,1], probs = 0.025)
z.itt = quantile(ub.df[,2], probs = 0.025)
z.mpp = quantile(ub.df[,3], probs = 0.025)
z.pp = quantile(ub.df[,4], probs = 0.025)

getoutcome<-function(vector.outcome1, vector.outcome0, intervention){
  outcome<-c()
  for (i in 1:length(vector.outcome0)){
    if (intervention[i]==1) {outcome[[i]]=vector.outcome1[i]} else {outcome[[i]]=vector.outcome0[i]}
  }
  return(unlist(outcome))
}
p.stdcare=0.4
#simulate and derive treatment effect 
for(l in 1:nIterations) { 
  tryCatch({
    
    id=seq(1,(2*n), by=1) #create participant id  
    randomisation=c(rep(1,n), rep(0,n)) #randomisation
    confounder = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment,1-p.experiment))  #probability of outcome if intervention = 1
    outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))       #probability of outcome if intervention = 0
    
    #INTERVENTION dependent on compliance 
    experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1-comply.experiment,comply.experiment))
    stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(comply.stdcare,1-comply.stdcare))
    intervention = c(experiment.intervention,stdcare.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome1, outcome0, intervention)
    
    simdata = matrix(data=c(id,randomisation,confounder,intervention,outcome), nrow=(2*n))
    
    ## intention to treat 
    pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
    pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
    eff.itt = pz1.value-pz0.value
    sd.itt = sqrt(pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n)

    ## per protocol 
    pp = simdata[which(simdata[,2]==simdata[,4]),] # perprotocol population
    p.experiment.vector = pp[which(pp[,2]==1),][,5]
    p.experiment.value = mean(p.experiment.vector)
    p.stdcare.vector = pp[which(pp[,2]==0),][,5]
    p.stdcare.value = mean(p.stdcare.vector) 
    eff.pp = p.experiment.value-p.stdcare.value
    sd.pp = sqrt(p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector))
    
    ## inverse probability weights on per protocol patients 
    pp = as.data.frame(pp)
    colnames(pp) = c('id','randomisation','confounder','intervention','outcome')
    ipwmodel = glm(intervention~confounder,family=binomial(link="logit"), data = pp) #calculate denominators used in inverse probability weights
    score = predict(ipwmodel, type="response")
    weight = pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
    outcomemodel = glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
    eff.mpp = coef(outcomemodel)[2]
    sd.mpp = sqrt(diag(vcovHC(outcomemodel)))[2]

    # iv with 2 stage regression
    asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
    eff.iv = summary(asmm)$coefficients[2,1]
    sd.iv = summary(asmm)$coefficients [2,2]
    
    estimates[[l]] = (c(eff.iv, eff.itt, eff.mpp, eff.pp))
    sd[[l]] = c(sd.iv, sd.itt, sd.mpp, sd.pp)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
}

# mean of power from iterated data 
est.df.altc=matrix(unlist(estimates), byrow=T,ncol = 4)
sd.df.altc=matrix(unlist(sd), byrow = T, ncol = 4)
ub.df.altc = est.df.altc + z*sd.df.altc
mean(ub.df.altc[,1]<z.iv)
mean(ub.df.altc[,2]<z.itt)
mean(ub.df.altc[,3]<z.mpp)
mean(ub.df.altc[,4]<z.pp)

plot(density(ub.df.altc[,1]))
lines(density(ub.dfc[,1]))
################################################
################################################

##debugging 
###compare SE of the 2 types of non adherence 
analysis <- function(simdata, type, n, z, NImargin, confounding.type){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
  p.experiment.vector= pp[which(pp[,2]==1), ][,5]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,5]
  p.stdcare.value= mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
  ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
  eff.iv = summary(asmm)$coefficients[2,1] 
  
  if (type=="bias"){
    
    return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
    
  } else { #type 1 error and power 
    
    ## intention to treat 
    var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
    CI.itt =  eff.itt + z * sqrt(var.eff.itt)
    itt = CI.itt < NImargin
    itt.se = sqrt(var.eff.itt)
    
    ## per protocol 
    var.eff.pp =   p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
    CI.pp =  eff.pp + z * sqrt(var.eff.pp)
    ppp =  CI.pp < NImargin
    ppp.se = sqrt(var.eff.pp)
    
    ## inverse probability weights
    se = summary(outcomemodel)$coefficients[2,2]
    CI.mpp = eff.mpp + z * se
    mpp =  CI.mpp < NImargin
    mpp.se = se
    
    # iv with 2 stage regression
    se = summary(asmm)$coefficients [2,2]
    CI.iv = eff.iv + z * se
    iv =  CI.iv < NImargin
    iv.se = se
    
    return(c(iv.se, itt.se, mpp.se, ppp.se))
  }
}
type1.nonconfounding <- function(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop, cross.over){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, nonadhere.pop=nonadhere.pop, cross.over = cross.over, i=i)
        .estimate[[l]] = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    # mean of se from iterated data 
    estimate[[i]] = colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE))
  }
  
  se.df = matrix(unlist(estimate), ncol = 4, byrow = T)
  return(se.df)
  
}
type1.confounding <- function(n, p.experiment, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop, i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]] = analysis(simdata=simdata, type='type1',n=n, z=z, NImargin=NImargin, confounding.type = 'confounding')
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of type 1 error from iterated data 
    estimate[[i]] = colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE),na.rm = TRUE)
  }
  
  #put into dataframe- each column contains n iterations of upper bound of 96% CI from each method
  se.df = matrix(unlist(estimate), ncol = 4, byrow = T)
  
  return(se.df)
  
}

se.nonconfounding = type1.nonconfounding(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop = 'both', cross.over=T)
se.confounding = type1.confounding(n, p.experiment, p.alt, NImargin, confounder.intervention='Increase likelihood', confounder.outcome='Decrease likelihood', interval, nIterations, cross.over=T, nonadhere.pop='both')
##SEs are the same for both! 

###compare estimates for the 2 types of non adherence 
analysis <- function(simdata, type, n, z, NImargin, confounding.type){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
  p.experiment.vector= pp[which(pp[,2]==1), ][,5]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,5]
  p.stdcare.value= mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
  ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
  eff.iv = summary(asmm)$coefficients[2,1] 
  
  if (type=="bias"){
    
    return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
    
  } else { #type 1 error and power 
    
    ## intention to treat 
    var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
    CI.itt =  eff.itt + z * sqrt(var.eff.itt)
    itt = CI.itt < NImargin
    itt.se = sqrt(var.eff.itt)
    
    ## per protocol 
    var.eff.pp =   p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
    CI.pp =  eff.pp + z * sqrt(var.eff.pp)
    ppp =  CI.pp < NImargin
    ppp.se = sqrt(var.eff.pp)
    
    ## inverse probability weights
    se = summary(outcomemodel)$coefficients[2,2]
    CI.mpp = eff.mpp + z * se
    mpp =  CI.mpp < NImargin
    mpp.se = se
    
    # iv with 2 stage regression
    se = summary(asmm)$coefficients [2,2]
    CI.iv = eff.iv + z * se
    iv =  CI.iv < NImargin
    iv.se = se
    
    return(c(iv.se, itt.se, mpp.se, ppp.se))
  }
}
type1.nonconfounding <- function(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop, cross.over){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate = c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, nonadhere.pop=nonadhere.pop, cross.over = cross.over, i=i)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    # mean of estimates from iterated data 
    estimate[[i]] = colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE))
    # mean of upper bounds from iterated data 
    ub[[i]] = colMeans(matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE))
  }
  
  estimates.df = matrix(unlist(estimate), ncol = 4, byrow = T)
  ub.df = matrix(unlist(ub), ncol = 4, byrow = T)
  return(list(estimates.df, ub.df))
  
}
type1.confounding <- function(n, p.experiment, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate= c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop, i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of estimates from iterated data 
    estimate[[i]] = colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE))
    # mean of upper bounds from iterated data 
    ub[[i]] = colMeans(matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE))
  }
  
  estimates.df = matrix(unlist(estimate), ncol = 4, byrow = T)
  ub.df = matrix(unlist(ub), ncol = 4, byrow = T)
  return(list(estimates.df, ub.df))

}

estimates.nonconfounding = type1.nonconfounding(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop='both', cross.over=T)
estimates.confounding = type1.confounding(n, p.experiment, p.alt, NImargin, confounder.intervention='Increase likelihood', confounder.outcome='Decrease likelihood', interval, nIterations, cross.over=T, nonadhere.pop='both')

View(estimates.nonconfounding[1][[1]])
estimates.confounding[1][[1]]

View(estimates.nonconfounding[2][[1]])
estimates.confounding[2][[1]]
###both estimates and the upper bounds are the same (to 2 decimal points)!!! 

###compare distribution of the upper bounds for the 2 types of non adherence 
type1.nonconfounding <- function(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop, cross.over){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate = c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, nonadhere.pop=nonadhere.pop, cross.over = cross.over, i=i)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    # mean of estimates from iterated data 
    estimate[[i]] = matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE)
    # mean of upper bounds from iterated data 
    ub[[i]] = matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE)
  }
  return(list(estimate, ub))
  
}
type1.confounding <- function(n, p.experiment, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate= c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop, i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of estimates from iterated data 
    estimate[[i]] = matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE)
    # mean of upper bounds from iterated data 
    ub[[i]] = matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE)
  }

  return(list(estimate, ub))
  
}

distributions.nonconfounding = type1.nonconfounding(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop='both', cross.over=T)
distributions.confounding = type1.confounding(n, p.experiment, p.alt, NImargin, confounder.intervention='Increase likelihood', confounder.outcome='Decrease likelihood', interval, nIterations, cross.over=T, nonadhere.pop='both')

##plot distributions of the upper bounds
ub.iv.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,1] #distribution of the upper bounds when compliance is 100%, for IV 
ub.iv.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,1]
dat = data.frame(dens = c(ub.iv.fullcomply.nc, ub.iv.fullcomply.c), lines = rep(c("iv.nc", "iv.c"), each = nIterations) )
a=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.itt.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,2] #distribution of the upper bounds when compliance is 100%, for ITT
ub.itt.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,2]
dat = data.frame(dens = c(ub.itt.fullcomply.nc, ub.itt.fullcomply.c), lines = rep(c("itt.nc", "itt.c"), each = nIterations) )
b=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.ipw.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,3] #distribution of the upper bounds when compliance is 100%, for IPW
ub.ipw.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,3]
dat = data.frame(dens = c(ub.ipw.fullcomply.nc, ub.ipw.fullcomply.c), lines = rep(c("ipw.nc", "ipw.c"), each = nIterations) )
c=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)
  

ub.pp.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,4] #distribution of the upper bounds when compliance is 100%, for PP
ub.pp.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,4]
dat = data.frame(dens = c(ub.pp.fullcomply.nc, ub.pp.fullcomply.c), lines = rep(c("pp.nc", "pp.c"), each = nIterations) )
d=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ggarrange(a, b, c, d, ncol=2, nrow = 2)

##plot distributions of the estimates 
ub.iv.fullcomply.nc = distributions.nonconfounding[[1]][[length(interval)]][,1] #distribution of the upper bounds when compliance is 100%, for IV 
ub.iv.fullcomply.c = distributions.confounding[[1]][[length(interval)]][,1]
dat = data.frame(dens = c(ub.iv.fullcomply.nc, ub.iv.fullcomply.c), lines = rep(c("iv.nc", "iv.c"), each = nIterations) )
a=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('risk difference estimates')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.itt.fullcomply.nc = distributions.nonconfounding[[1]][[length(interval)]][,2] #distribution of the upper bounds when compliance is 100%, for ITT
ub.itt.fullcomply.c = distributions.confounding[[1]][[length(interval)]][,2]
dat = data.frame(dens = c(ub.itt.fullcomply.nc, ub.itt.fullcomply.c), lines = rep(c("itt.nc", "itt.c"), each = nIterations) )
b=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('risk difference estimates')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.ipw.fullcomply.nc = distributions.nonconfounding[[1]][[length(interval)]][,3] #distribution of the upper bounds when compliance is 100%, for IPW
ub.ipw.fullcomply.c = distributions.confounding[[1]][[length(interval)]][,3]
dat = data.frame(dens = c(ub.ipw.fullcomply.nc, ub.ipw.fullcomply.c), lines = rep(c("ipw.nc", "ipw.c"), each = nIterations) )
c=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + 
  xlab('risk difference estimates')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)


ub.pp.fullcomply.nc = distributions.nonconfounding[[1]][[length(interval)]][,4] #distribution of the upper bounds when compliance is 100%, for PP
ub.pp.fullcomply.c = distributions.confounding[[1]][[length(interval)]][,4]
dat = data.frame(dens = c(ub.pp.fullcomply.nc, ub.pp.fullcomply.c), lines = rep(c("pp.nc", "pp.c"), each = nIterations) )
d=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('risk difference estimates')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ggarrange(a, b, c, d, ncol=2, nrow = 2)
###C has a narrower distribution than NC upper bounds 

####broadening z value so the 2 distributions are similar 
type1.nonconfounding <- function(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop, cross.over){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate = c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, nonadhere.pop=nonadhere.pop, cross.over = cross.over, i=i)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'nonconfounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    # mean of estimates from iterated data 
    estimate[[i]] = matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE)
    # mean of upper bounds from iterated data 
    ub[[i]] = matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE)
  }
  return(list(estimate, ub))
  
}
type1.confounding <- function(n, p.experiment, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.97) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate = ub = c() #output from each simulation
  estimate = ..estimate= c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop, i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        #estimates
        .estimate[[l]] = analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        #upper bounds 
        se = analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin, confounding.type = 'confounding')
        ..estimate[[l]] = .estimate[[l]] + z * se
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of estimates from iterated data 
    estimate[[i]] = matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE)
    # mean of upper bounds from iterated data 
    ub[[i]] = matrix(as.numeric(unlist(..estimate)), ncol = length(analysis.method), byrow = TRUE)
  }
  
  return(list(estimate, ub))
  
}

distributions.nonconfounding = type1.nonconfounding(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop='both', cross.over=T)
distributions.confounding = type1.confounding(n, p.experiment, p.alt, NImargin, confounder.intervention='Increase likelihood', confounder.outcome='Decrease likelihood', interval, nIterations, cross.over=T, nonadhere.pop='both')

ub.iv.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,1] #distribution of the upper bounds when compliance is 100%, for IV 
ub.iv.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,1]
dat = data.frame(dens = c(ub.iv.fullcomply.nc, ub.iv.fullcomply.c), lines = rep(c("iv.nc", "iv.c"), each = nIterations) )
a=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.itt.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,2] #distribution of the upper bounds when compliance is 100%, for ITT
ub.itt.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,2]
dat = data.frame(dens = c(ub.itt.fullcomply.nc, ub.itt.fullcomply.c), lines = rep(c("itt.nc", "itt.c"), each = nIterations) )
b=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ub.ipw.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,3] #distribution of the upper bounds when compliance is 100%, for IPW
ub.ipw.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,3]
dat = data.frame(dens = c(ub.ipw.fullcomply.nc, ub.ipw.fullcomply.c), lines = rep(c("ipw.nc", "ipw.c"), each = nIterations) )
c=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)


ub.pp.fullcomply.nc = distributions.nonconfounding[[2]][[length(interval)]][,4] #distribution of the upper bounds when compliance is 100%, for PP
ub.pp.fullcomply.c = distributions.confounding[[2]][[length(interval)]][,4]
dat = data.frame(dens = c(ub.pp.fullcomply.nc, ub.pp.fullcomply.c), lines = rep(c("pp.nc", "pp.c"), each = nIterations) )
d=ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

ggarrange(a, b, c, d, ncol=2, nrow = 2)
##0.970 instead of 0.975 still too small 

########################################################################
############using SD from the simulated disgribution instead############
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/") #set working directory 
source(file = 'simulationcode_200120.R')

interval=1

bias.nonconfounding <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over){  
  z = qnorm(0.95)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, cross.over = cross.over, i=i)
        .estimate[[l]] = analysis.estimate(simdata=simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp 
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
 df=as.data.frame(matrix(unlist(estimate), byrow = T, ncol = 4))
 
 # for estimates 
 itt = (unlist(df[2]))
 iv =(unlist(df[1] ))
 mpp = (unlist(df[3] ))
 pp =(unlist(df[4] ))
 
  # for upper bounds
 
 
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv =(unlist(df[1] ))+ z * sd(unlist(df[1]))
  # mpp = (unlist(df[3] ))+ z * sd(unlist(df[3]))
  # pp =(unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  
  return(list(itt,pp,mpp,iv))
} 

bias.confounding <- function(n, p.experiment, p.stdcare, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  

  z= qnorm(0.95)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]] = analysis.estimate(simdata=simdata)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  
  # for estimates 
  itt = (unlist(df[2])) 
  iv = (unlist(df[1] )) 
  mpp =(unlist(df[3] )) 
  pp = (unlist(df[4] ))
  
  # for upper bounds
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv = (unlist(df[1] )) + z * sd(unlist(df[1]))
  # mpp =(unlist(df[3] )) + z * sd(unlist(df[3]))
  # pp = (unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  return(list(itt, iv, mpp, pp))
  
}

nIterations=100000
nc.ub.df=bias.nonconfounding(n=505, p.experiment=0.4, p.stdcare=0.4, p.alt=0.5, nIterations=nIterations, interval=interval, nonadhere.pop='both', cross.over=T)
c.ub.df=bias.confounding(n=505, p.experiment=0.4, p.stdcare=0.4, p.alt=0.45, NImargin=0.1, confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', interval=interval, nIterations=nIterations, cross.over=T, nonadhere.pop='both')

plot(density(c.ub.df[[1]]))
lines(density(nc.ub.df[[1]]))

dat = data.frame(dens = c(unlist(nc.ub.df[[4]]), unlist(c.ub.df[[4]])), lines = rep(c("nc", "c"), each =1000) )
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

#######################################
#######################################
#######################################
#######################################
# change of SD with increasing iterations 
interval=1

### when SD is calculated 
sd.nonconfounding.cal <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over){  
  
  z = qnorm(0.975)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, cross.over = cross.over, i=i)
        
        ## intention to treat 
        pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
        pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
        eff.itt = pz1.value-pz0.value
        ## intention to treat 
        var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
        itt =  sqrt(var.eff.itt)
        
        ## per protocol 
        pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
        p.experiment.vector= pp[which(pp[,2]==1), ][,5]
        p.experiment.value= mean(p.experiment.vector)
        p.stdcare.vector= pp[which(pp[,2]==0),][,5]
        p.stdcare.value= mean(p.stdcare.vector) 
        eff.pp = p.experiment.value-p.stdcare.value
        ## per protocol 
        var.eff.pp =   p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
        ppp =  sqrt(var.eff.pp)
        
        ## inverse probability weights on per protocol patients 
        pp=as.data.frame(pp)
        colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
        ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
        score=predict(ipwmodel, type="response")
        weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
        outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
        eff.mpp=coef(outcomemodel)[2]
        ## inverse probability weights
        ipw = summary(outcomemodel)$coefficients[2,2]
        
        # iv with 2 stage regression
        asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
        eff.iv = summary(asmm)$coefficients[2,1]
        # iv with 2 stage regression
        iv = summary(asmm)$coefficients [2,2]
        
        .estimate [[l]] = c(iv, itt, ipw, ppp)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df=as.data.frame(matrix(unlist(estimate), byrow = T, ncol = 4))
  
  return(df)
} 
sd.confounding.cal <- function(n, p.experiment, p.stdcare, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  z= qnorm(0.975)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        ## intention to treat 
        pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
        pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
        eff.itt = pz1.value-pz0.value
        ## intention to treat 
        var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
        itt =  sqrt(var.eff.itt)
        
        ## per protocol 
        pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
        p.experiment.vector= pp[which(pp[,2]==1), ][,5]
        p.experiment.value= mean(p.experiment.vector)
        p.stdcare.vector= pp[which(pp[,2]==0),][,5]
        p.stdcare.value= mean(p.stdcare.vector) 
        eff.pp = p.experiment.value-p.stdcare.value
        ## per protocol 
        var.eff.pp =   p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
        ppp =  sqrt(var.eff.pp)
        
        ## inverse probability weights on per protocol patients 
        pp=as.data.frame(pp)
        colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
        ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
        score=predict(ipwmodel, type="response")
        weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
        outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
        eff.mpp=coef(outcomemodel)[2]
        ## inverse probability weights
        ipw = summary(outcomemodel)$coefficients[2,2]
        
        # iv with 2 stage regression
        asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
        eff.iv = summary(asmm)$coefficients[2,1]
        # iv with 2 stage regression
        iv = summary(asmm)$coefficients [2,2]
        
        .estimate[[l]] = c(iv, itt, ipw, ppp)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))

  return(df)
  
}

### when SD is from the distribution 
sd.nonconfounding.sim <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over){  
  z = qnorm(0.975)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, cross.over = cross.over, i=i)
        .estimate[[l]] = analysis.estimate(simdata=simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp 
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df=as.data.frame(matrix(unlist(estimate), byrow = T, ncol = 4))
  
  # for estimates 
  itt = sd(unlist(df[2]))
  iv =sd(unlist(df[1] ))
  mpp = sd(unlist(df[3] ))
  pp =sd(unlist(df[4]))
  
  # for upper bounds
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv =(unlist(df[1] ))+ z * sd(unlist(df[1]))
  # mpp = (unlist(df[3] ))+ z * sd(unlist(df[3]))
  # pp =(unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  
  return(list(itt,pp,mpp,iv))
} 
sd.confounding.sim <- function(n, p.experiment, p.stdcare, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  z= qnorm(0.975)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]] = analysis.estimate(simdata=simdata)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  
  # for estimates 
  itt = sd(unlist(df[2]))
  iv = sd(unlist(df[1]))
  mpp = sd(unlist(df[3]))
  pp = sd(unlist(df[4]))
  
  # for upper bounds
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv = (unlist(df[1] )) + z * sd(unlist(df[1]))
  # mpp =(unlist(df[3] )) + z * sd(unlist(df[3]))
  # pp = (unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  return(list(itt, iv, mpp, pp))
  
}

sdnccal = sd.nonconfounding.cal(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, nIterations=5000, interval=interval, nonadhere.pop='both', cross.over=T)
sdccal = sd.confounding.cal(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, NImargin=0.1, confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', interval=interval, nIterations=5000, cross.over=T, nonadhere.pop='both')
sd(sdccal[[1]])
sd(sdnccal[[1]])
plot(density(unlist(sdccal[1])))
lines(density(unlist(sdnccal[1])))

sdncsim = sd.nonconfounding.sim(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, nIterations=10000, interval=interval, nonadhere.pop='both', cross.over=T)
sdcsim= sd.confounding.sim(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, NImargin=0.1, confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', interval=interval, nIterations=10000, cross.over=T, nonadhere.pop='both')
mean(unlist(sdcsim))
mean(unlist(sdncsim))

###EXPlore effect of iterations on estimates 
est.nonconfounding.sim <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over){  
  z = qnorm(0.975)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, cross.over = cross.over, i=i)
        .estimate[[l]] = analysis.estimate(simdata=simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp 
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df=as.data.frame(matrix(unlist(estimate), byrow = T, ncol = 4))
  
  # for estimates 
  itt = (unlist(df[2]))
  iv =(unlist(df[1] ))
  mpp = (unlist(df[3] ))
  pp =(unlist(df[4]))
  
  # for upper bounds
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv =(unlist(df[1] ))+ z * sd(unlist(df[1]))
  # mpp = (unlist(df[3] ))+ z * sd(unlist(df[3]))
  # pp =(unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  
  return(list(itt,pp,mpp,iv))
} 
est.confounding.sim <- function(n, p.experiment, p.stdcare, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  z= qnorm(0.975)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]] = analysis.estimate(simdata=simdata)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  
  # for estimates 
  itt = (unlist(df[2]))
  iv = (unlist(df[1]))
  mpp = (unlist(df[3]))
  pp = (unlist(df[4]))
  
  # for upper bounds
  # itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  # iv = (unlist(df[1] )) + z * sd(unlist(df[1]))
  # mpp =(unlist(df[3] )) + z * sd(unlist(df[3]))
  # pp = (unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  return(list(itt, iv, mpp, pp))
  
}
estncsim = est.nonconfounding.sim(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, nIterations=20000, interval=interval, nonadhere.pop='both', cross.over=T)
estcsim= est.confounding.sim(n=505, p.experiment=0.4, p.stdcare=0.3, p.alt=0.5, NImargin=0.1, confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', interval=interval, nIterations=20000, cross.over=T, nonadhere.pop='both')
mean(unlist(estncsim[[1]]))
mean(unlist(estcsim[[1]]))
sd(unlist(estncsim[[1]]))
sd(unlist(estcsim[[1]]))
plot(density(unlist(estcsim[1])))
lines(density(unlist(estncsim[1])))

###explore effect on upper bounds 
source(file = 'simulationcode_200120.R')
interval=1
ub.nonconfounding.sim <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over){  
  z = qnorm(0.975)
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.nonconfounding(n=n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, cross.over = cross.over, i=i)
        .estimate[[l]] = analysis.estimate(simdata=simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp 
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df=as.data.frame(matrix(unlist(estimate), byrow = T, ncol = 4))
  
  # for estimates 
  # itt = (unlist(df[2]))
  # iv =(unlist(df[1] ))
  # mpp = (unlist(df[3] ))
  # pp =(unlist(df[4]))
  
  # for upper bounds
  itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  iv =(unlist(df[1] ))+ z * sd(unlist(df[1]))
  mpp = (unlist(df[3] ))+ z * sd(unlist(df[3]))
  pp =(unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  
  return(list(itt,pp,mpp,iv))
} 
ub.confounding.sim <- function(n, p.experiment, p.stdcare, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  z= qnorm(0.975)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]] = analysis.estimate(simdata=simdata)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  
  # for estimates 
  # itt = (unlist(df[2]))
  # iv = (unlist(df[1]))
  # mpp = (unlist(df[3]))
  # pp = (unlist(df[4]))
  # 
  # for upper bounds
  itt = (unlist(df[2])) + z * sd(unlist(df[2]))
  iv = (unlist(df[1] )) + z * sd(unlist(df[1]))
  mpp =(unlist(df[3] )) + z * sd(unlist(df[3]))
  pp = (unlist(df[4] ))+ z * sd(unlist(df[4]))
  
  return(list(itt, iv, mpp, pp))
  
}

ubncsim = ub.nonconfounding.sim(n=505, p.experiment=0.4, p.stdcare=0.4, p.alt=0.5, nIterations=3000, interval=interval, nonadhere.pop='both', cross.over=T)
ubcsim= ub.confounding.sim(n=505, p.experiment=0.4, p.stdcare=0.4, p.alt=0.5, NImargin=0.1, confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', interval=interval, nIterations=3000, cross.over=T, nonadhere.pop='both')
plot(density(unlist(ubncsim[1])))
summary(unlist(ubcsim[1]))
plot(density(unlist(ubcsim[1])))

dat = data.frame(dens = c(unlist(ubncsim[[1]]), unlist(ubcsim[[1]])), lines = rep(c("nc", "c"), each =3000) )
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+ 
  xlab('95%CI upper bounds')+
  geom_vline(xintercept = 0.1, linetype="dotted", 
             color = "blue", size=1.5)

#relationship of shape1 and shape2
p = seq(0.4, 1, by=0.01)
shape2 = seq(0, 1, length.out=length(p))
shape1 = shape2*p/(1-p)
plot(x=shape2, y=shape1, type='l') #shape 1 increases exponentially with shape 2 as p approaches 1

#probabilities for outcome: 
p=0.5
shape2 = runif(1, min=0.5, max=3)
shape1 = shape2*p/(1-p)
hist(rbeta(10000,shape1, shape2))
#changing max from 1 to 100 = 
#<1 is concave(more people with small or big probabilties = LESS effect of confounders on outcome), 
#1 is uniform 
# >2 is normal distribution(more people have probabilities about middle = MORE effect of confounders on outcome)

#probabilities for intervention 
## to preserve positivity cannot have too concaved
## (a lot people with very low and many others with high probabilities) - 
## then taking up intervention will not be dependent on confounders 
p=0.6
shape2 = runif(1, min=0.5, max=3)
shape1 = shape2*p/(1-p)
hist(rbeta(10000,shape1, shape2))
# 2-10: 0.6 most normal increasingly right skewed until 1
# 0-1: 0.6 convave to right slant until 1 (too concaved)
# 0-2: too concave
# 1-2: 0.6 mostly normal (weak effect of confounder on intervention)
# 0.5-2:0.6 good mixture of normal/ flat and some mild convave but too skewed at 0.8



shape2.1 = seq(0, 1, length.out=length(p))
shape1.1 = shape2.1*p/(1-p)
shape2.10 = seq(0, 10, length.out=length(p))
shape1.10 = shape2.10*p/(1-p)
plot(x=p, y=shape1.1/shape2.1, type='l')
lines(x=p, y=shape1.10/shape2.10, type='l', col='red')
## with increasing probability, shape1/shape2 increasing exponentially to infinity 
## with runif(1) and runif(min=2, max=10), no difference

plot(x=p, y=shape1.1 - shape2.1, type='l')
lines(x=p, y=shape1.10 - shape2.10, type='l', col='red')
## with increasing probability, shape 1 - shape 2 increases faster when shape 2 is
## drawn from 0 to 10 compared to 0 to 1 
## 
# in generating probability for outcome - does not matter shape or direction
# in generating probability for intervention - need to ensure positivity, hence 
# shape 2 is drawn 0-1 and shape 1 is from 2-10 
hist(rbeta(10000, 10, 2)) 
hist(rbeta(10000, 2, 10)) 
hist(rbeta(10000, 2, 2)) 
hist(rbeta(10000, 10, 10)) 

##when shape1 and shape2 are the same 
hist(rbeta(10000, 10, 10)) #normal
hist(rbeta(10000, 1, 1)) #flat
hist(rbeta(10000, 0.4, 0.4)) #concave
hist(rbeta(10000, 0.2, 0.2)) #more concave

#when shape1 >> shape2 
hist(rbeta(10000, 1, 0.5)) #ratio 2
hist(rbeta(10000, 2, 1)) #more gradual increase to probability of 1

hist(rbeta(10000, 1, 0.1)) #ratio 10
hist(rbeta(10000, 10, 1)) #more gradual increase to probability of 1

#when shape 1 < shape 2
hist(rbeta(10000, 0.5, 1)) #ratio 0.5
hist(rbeta(10000, 2, 4)) #more gradual increase to probability of 1

##ratio defines the direction, difference defines the shape

########################################################
########################################################
########################################################
########################################################
power.test <- function (n, p.experiment, p.stdcare, p.alt, cross.over, nonadhere.pop, confounder.intervention, confounder.outcome, nIterations, interval, NImargin){
  
  z= qnorm(0.975)
  #make up vectors for simulations 
  estimate = c()  #for saving output from each interval
  .estimate = c() #output from each simulation
  
  #generate critical value from null hypothesis 
  #for (i in 1:length(interval)) {
    for (l in 1:nIterations){
      simdata = simdata.confounding(n=n, p.experiment = p.experiment, p.stdcare = p.experiment - NImargin, p.alt = p.alt, cross.over = cross.over, nonadhere.pop = nonadhere.pop, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention, i=1)
      
      ## intention to treat 
      pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
      pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
      eff.itt = pz1.value-pz0.value
      ## intention to treat 
      var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
      itt =  eff.itt + z * sqrt(var.eff.itt)
      
      ## per protocol 
      pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
      p.experiment.vector= pp[which(pp[,2]==1), ][,5]
      p.experiment.value= mean(p.experiment.vector)
      p.stdcare.vector= pp[which(pp[,2]==0),][,5]
      p.stdcare.value= mean(p.stdcare.vector) 
      eff.pp = p.experiment.value-p.stdcare.value
      ## per protocol 
      var.eff.pp = p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
      ppp =  eff.pp + z * sqrt(var.eff.pp)

      ## inverse probability weights on per protocol patients 
      pp=as.data.frame(pp)
      colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
      ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
      score=predict(ipwmodel, type="response")
      weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
      outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
      eff.mpp=coef(outcomemodel)[2]
      ## inverse probability weights
      ipw = eff.mpp + z * summary(outcomemodel)$coefficients[2,2]
      
      # iv with 2 stage regression
      asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
      eff.iv = summary(asmm)$coefficients[2,1]
      # iv with 2 stage regression
      iv = eff.iv + z * summary(asmm)$coefficients [2,2]
      
      .estimate[[l]] = c(iv, itt, ipw, ppp)
    }
    #save results of every interval
    estimate = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  #}
  df.z = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  z.iv = quantile(df.z[,1], probs = 0.025) 
  z.itt = quantile(df.z[,2], probs = 0.025)
  z.ipw = quantile(df.z[,3], probs = 0.025)
  z.pp = quantile(df.z[,4], probs = 0.025)

  ..estimate<-sd<-c()
  #simulate and derive treatment effect based on alternative hypothesis 
  #for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata = simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop=nonadhere.pop,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        ## intention to treat 
        pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
        pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
        eff.itt = pz1.value-pz0.value
        ## intention to treat 
        var.eff.itt =  pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
        itt =  eff.itt +  z * sqrt(var.eff.itt)
        sd.itt = sqrt(var.eff.itt)
        
        ## per protocol 
        pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
        p.experiment.vector= pp[which(pp[,2]==1), ][,5]
        p.experiment.value= mean(p.experiment.vector)
        p.stdcare.vector= pp[which(pp[,2]==0),][,5]
        p.stdcare.value= mean(p.stdcare.vector) 
        eff.pp = p.experiment.value-p.stdcare.value
        ## per protocol 
        var.eff.pp =   p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
        ppp =  eff.pp +  z * sqrt(var.eff.pp)
        sd.pp =sqrt(var.eff.pp)
        
        ## inverse probability weights on per protocol patients 
        pp=as.data.frame(pp)
        colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
        ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
        score=predict(ipwmodel, type="response")
        weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
        outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
        eff.mpp=coef(outcomemodel)[2]
        ## inverse probability weights
        ipw = eff.mpp +  z * summary(outcomemodel)$coefficients[2,2]
        sd.ipw = summary(outcomemodel)$coefficients[2,2]
        # iv with 2 stage regression
        asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
        eff.iv = summary(asmm)$coefficients[2,1]
        # iv with 2 stage regression
        iv = eff.iv +  z * summary(asmm)$coefficients [2,2]
        sd.iv=summary(asmm)$coefficients [2,2]
        
        .estimate[[l]] = c(iv, itt, ipw, ppp)
        ..estimate[[l]] = c(eff.iv, eff.itt, eff.mpp, eff.pp)
        sd[[l]]=c(sd.iv, sd.itt, sd.ipw, sd.pp)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  #}
  estimate.=matrix(unlist(..estimate), ncol = length(analysis.method), byrow = TRUE)
  sd.df=matrix(unlist(sd), ncol = length(analysis.method), byrow = TRUE)
  df = as.data.frame(matrix(unlist(estimate), byrow = T, ncol=4))
  
  iv.power = mean(df[1]< z.iv)
  itt.power = mean(df[2]< z.itt)
  mpp.power = mean(df[3] < z.ipw)
  pp.power = mean(df[4] < z.pp)
  
  return(list(df.z, df, estimate., sd.df,c(iv.power,itt.power, mpp.power, pp.power)))
}
df=power.test(n=505, p.experiment=0.4, p.stdcare=0.4, p.alt=0.5, cross.over=T, nonadhere.pop='both', confounder.intervention='Increase likelihood', confounder.outcome='Increase likelihood', nIterations=1000, interval=1, NImargin=0.1)
null=df[[1]][3]
alt = df[[2]][3]
alt.est=df[[3]][,3]
sd=df[[4]][,3]
plot(density(unlist(null)), xlim=c(-0.05,0.25), ylim=c(0,30))
lines(density(unlist(alt)))  
lines(density(unlist(alt.est))) 
lines(density(unlist(sd))) 
df[[5]]
  