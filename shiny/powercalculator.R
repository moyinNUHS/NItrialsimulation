rm(list=ls())
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/") #set working directory 
source(file = 'simulationcode_070220.R')
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/shiny/samplesize_nonadherence/")
sourceCpp(file='rcpp.cpp')

####Deploy app 
# setwd('/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/shiny/samplesize_nonadherence/')
# require('rsconnect')
# rsconnect::setAccountInfo(name='moru',
#                           token='30819BAEDD492333CE7CD293F3B08D42',
#                           secret='j24TwHndsiybKqQZERFOt0s3L3ro8IG47/OgpLk/')
# deployApp(account="moru",appName="samplesize_nonadherence")

####parameters####
p.experiment = 0.4
p.stdcare = 0.4 
NImargin=0.1
adhere.stdcare = 1
adhere.experiment = 1
n = 505
confounder.intervention = 'Decrease likelihood'
confounder.outcome = 'Decrease likelihood'

###functions###
analysis.ub <- function(simdata) {
  
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
  eff.mpp= coef(outcomemodel)[2]
  sd.mpp = sqrt(diag(vcovHC(outcomemodel)))[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
  eff.iv = summary(asmm)$coefficients[2,1]
  sd.iv = summary(asmm)$coefficients [2,2]
  
  return(list (c(eff.itt, eff.pp, eff.mpp, eff.iv), 
               c(sd.itt, sd.pp, eff.mpp, eff.iv)))
}

getoutcome <- function(outcome0, outcome1, intervention){
  outcome = c()
  outcome.matrix = matrix(c(outcome0, outcome1), ncol = 2)
  for (i in 1:length(intervention)){
    outcome[i] = outcome.matrix[i, (intervention[i]+1)]
  }
  return(unlist(outcome))
}

##########################################
##Using null distribution to get Z value##
##########################################

########Non-confounding##########
nc.power <- function (n, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, NImargin){
  
  nIterations = 1000
  
  ###NULL HYPOTHESIS to get z value 
  id = seq(1,(2*n), by = 1) #create participant id  
  randomisation = sample(rep(0:1, n)) #randomisation
  confounder = rbeta(n = (2*n), shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1 
  
  means.null.nc = sds.null.nc = c()
  p.stdcare.null.nc = p.experiment - NImargin
  for (i in 1:nIterations) {
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome0 = rbinom(2*n, prob=p.stdcare.null.nc, size = 1)
    outcome1 = rbinom(2*n, prob=p.experiment, size = 1) 
    
    #INTERVENTION dependent on adherence
    intervention = rep(NA, 2*n)
    intervention[sample(which(randomisation==1), size = adhere.experiment * n)] = 1
    intervention[which(randomisation==1) %in% which(is.na(intervention))] = 0
    intervention[sample(which(randomisation==0), size = adhere.stdcare * n)] = 0
    intervention[which(randomisation==0) %in% which(is.na(intervention))] = 1
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=outcome0, outcome1=outcome1, intervention=intervention)
    simdata = matrix(data=c(id, randomisation, confounder, intervention, outcome), nrow=(2*n))
    
    ub = analysis.ub(simdata = simdata)
    means.null.nc[[i]] = ub[[1]]
  }
  means.df.null.nc = matrix(unlist(means.null.nc), ncol = 4, byrow = T)
  sds.null.nc = unlist(apply(means.df.null.nc, 2, sd))
  ub.df.null.nc = means.df.null.nc + rep(sds.null.nc, each = nIterations)
  z.vector.null.nc = c(quantile(ub.df.null.nc[,1], probs = 0.025), 
                       quantile(ub.df.null.nc[,2], probs = 0.025), 
                       quantile(ub.df.null.nc[,3], probs = 0.025), 
                       quantile(ub.df.null.nc[,4], probs = 0.025))
 
  ###ALTERNATIVE HYPOTHESIS for power 
  means.alt.nc = sds.alt.nc = c()
  for (i in 1:nIterations) {
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome0 = rbinom(2*n, prob=p.stdcare, size = 1)
    outcome1 = rbinom(2*n, prob=p.experiment, size = 1) 
    
    #INTERVENTION dependent on adherence
    intervention = rep(NA,2*n)
    intervention[sample(which(randomisation==1), size = adhere.experiment * n)] = 1
    intervention[which(randomisation==1) %in% which(is.na(intervention))] = 0
    intervention[sample(which(randomisation==0), size = adhere.stdcare * n)] = 0
    intervention[which(randomisation==0) %in% which(is.na(intervention))] = 1
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=outcome0, outcome1=outcome1, intervention=intervention)
    simdata = matrix(data=c(id, randomisation, confounder, intervention, outcome), nrow=(2*n))
    
    ub = analysis.ub(simdata = simdata)
    means.alt.nc[[i]] = ub[[1]]
  }
  means.df.alt.nc = matrix(unlist(means.alt.nc), ncol = 4, byrow = T)
  sds.alt.nc = unlist(apply(means.df.alt.nc, 2, sd))
  ub.df.alt.nc = means.df.alt.nc + rep(sds.alt.nc, each = nIterations)
  
  p.itt.nc = mean(ub.df.alt.nc[,1] < z.vector.null.nc[1])
  p.pp.nc = mean(ub.df.alt.nc[,2] < z.vector.null.nc[2])
  p.mpp.nc = mean(ub.df.alt.nc[,3] < z.vector.null.nc[3])
  p.iv.nc = mean(ub.df.alt.nc[,4] < z.vector.null.nc[4])
  
  return (list(c(p.itt.nc, p.pp.nc, p.mpp.nc, p.iv.nc), 
               ub.df.null.nc,
               ub.df.alt.nc, 
               z.vector.null.nc ))
}

p.505.nc = nc.power (n=505, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, NImargin)

###########Confounding############
c.power <- function (n, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, confounder.intervention, confounder.outcome, NImargin){
  
  nIterations = 1000 
  
  ###NULL HYPOTHESIS to get z value 
  means.null.c = sds.null.c = c()
  p.stdcare.null.c = p.experiment - NImargin
  
  id = seq(1,(2*n), by = 1) #create participant id  
  randomisation = sample(rep(0:1, n)) #randomisation
  
  shape2.outcome = runif(1)
  shape2.intervention = runif(1, min=2, max=10)
  
  for (i in 1:nIterations) {
    
    confounder = rbeta(n = (2*n), shape1 = 2,shape2 = 2) #confounder beta distribution ranging 0-1 
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    if (confounder.outcome=="Increase likelihood") {
      #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
      shape1 =  shape2.outcome*p.stdcare.null.c/(1-p.stdcare.null.c) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.stdcare, in increasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment) 
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.experiment, in increasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome
    } else {
      shape1 =  shape2.outcome*p.stdcare.null.c/(1-p.stdcare.null.c) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.stdcare, in decreasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment)
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    }
    
    d = matrix(data = c(id, randomisation, confounder), ncol = 3)
    d.ordered = matrix(cbind(d[order(d[,3]), ], outcome0, outcome1), ncol = 5) #order confounder in ascending order 
    d.grouped = rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),]) #group according to randomisation 
    
    #INTERVENTION dependent on randomisation and confounders
    if (confounder.intervention=="Increase likelihood") {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.experiment, in increasing order
      int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
      
      shape1 =  shape2.intervention*(1-adhere.stdcare)/(1-(1-adhere.stdcare)) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in increasing order
      cont.intervention =  rbinom(n, 1, prob=adhere.stdcare.ind)
    } else {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
      
      shape1 =  shape2.intervention*adhere.stdcare/(1-adhere.stdcare) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in decreasing order
      cont.intervention =  rbinom(n, 1, prob = 1-adhere.stdcare.ind)
    }
    
    intervention = c(int.intervention,cont.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=d.grouped[,4], outcome1=d.grouped[,5], intervention)
    simdata = matrix(cbind(d.grouped[,-5:-4], intervention, outcome), ncol=5)
    
    ub = analysis.ub(simdata = simdata)
    means.null.c[[i]] = ub[[1]]
  }
  means.df.null.c = matrix(unlist(means.null.c), ncol = 4, byrow = T)
  sds.null.c = unlist(apply(means.df.null.c, 2, sd))
  ub.df.null.c = means.df.null.c + rep(sds.null.c, each = nIterations)
  z.vector.null.c = c(quantile(ub.df.null.c[,1], probs = 0.025), 
                      quantile(ub.df.null.c[,2], probs = 0.025), 
                      quantile(ub.df.null.c[,3], probs = 0.025), 
                      quantile(ub.df.null.c[,4], probs = 0.025))
  
  ###ALTERNATIVE HYPOTHESIS for power 
  means.alt.c = sds.alt.c = c()
  
  shape2.outcome = runif(1)
  shape2.intervention = runif(1, min=2, max=10)
  
  for (i in 1:nIterations) {
    confounder = rbeta(n = (2*n), shape1 = 2,shape2 = 2) #confounder beta distribution ranging 0-1 
    
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    if (confounder.outcome=="Increase likelihood") {
      #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
      shape1 =  shape2.outcome*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.stdcare, in increasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment) 
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.experiment, in increasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome
    } else {
      shape1 =  shape2.outcome*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.stdcare, in decreasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment)
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    }
    
    d = matrix(data = c(id, randomisation, confounder), ncol = 3)
    d.ordered = matrix(cbind(d[order(d[,3]), ], outcome0, outcome1), ncol = 5) #order confounder in ascending order 
    d.grouped = rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),]) #group according to randomisation 
    
    #INTERVENTION dependent on randomisation and confounders
    if (confounder.intervention=="Increase likelihood") {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.experiment, in increasing order
      int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
      
      shape1 =  shape2.intervention*(1-adhere.stdcare)/(1-(1-adhere.stdcare)) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in increasing order
      cont.intervention =  rbinom(n, 1, prob=adhere.stdcare.ind)
    } else {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
      
      shape1 =  shape2.intervention*adhere.stdcare/(1-adhere.stdcare) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in decreasing order
      cont.intervention =  rbinom(n, 1, prob = 1-adhere.stdcare.ind)
    }
    
    intervention = c(int.intervention,cont.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=d.grouped[,4], outcome1=d.grouped[,5], intervention)
    simdata = matrix(cbind(d.grouped[,-5:-4], intervention, outcome), ncol=5)
    
    ub = analysis.ub(simdata = simdata)
    means.alt.c[[i]] = ub[[1]]
  }
  means.df.alt.c = matrix(unlist(means.alt.c), ncol = 4, byrow = T)
  sds.alt.c = unlist(apply(means.df.alt.c, 2, sd))
  ub.df.alt.c = means.df.alt.c + rep(sds.alt.c, each = nIterations)
  
  p.itt.alt.c = mean(ub.df.alt.c[,1] < z.vector.null.c[1])
  p.pp.alt.c = mean(ub.df.alt.c[,2] < z.vector.null.c[2])
  p.mpp.alt.c = mean(ub.df.alt.c[,3] < z.vector.null.c[3])
  p.iv.alt.c = mean(ub.df.alt.c[,4] < z.vector.null.c[4])
  
  return (list(c(p.itt.alt.c, p.pp.alt.c, p.mpp.alt.c, p.iv.alt.c), 
               ub.df.null.c, 
               ub.df.alt.c,
               z.vector.null.c))
}

p.505.c = c.power (n=600, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, confounder.intervention, confounder.outcome, NImargin)

##compare null distribution 
p.505.nc[[1]]  
p.505.c[[1]]

ub.df.null.nc = p.505.nc[[2]]  
ub.df.alt.nc = p.505.nc [[3]]
ub.df.null.c = p.505.c[[2]]
ub.df.alt.c = p.505.c[[3]]
plot(density(ub.df.null.nc[,1]), xlim = c(-0.07, 0.25))
lines(density(ub.df.null.c[,1]), col ='red')
abline(v= mean(ub.df.null.nc[,1]), lwd=1, lty=2)
abline(v= mean(ub.df.null.c[,1]), col='red', lwd=1, lty=2)
lines(x=c(mean(ub.df.null.nc[,1]), mean(ub.df.null.nc[,1])+sd(ub.df.null.nc[,1])), col='grey', y=c(3,3), lwd=1)
lines(x=c(mean(ub.df.null.c[,1]), mean(ub.df.null.c[,1])+sd(ub.df.null.c[,1])), y=c(3.1,3.1), col='pink',lwd=1)
text(round(sd(ub.df.null.c[,1]), digits = 4), x=mean(ub.df.null.c[,1])+0.02, y=3.4)
text(round(sd(ub.df.null.nc[,1]), digits = 4), x=mean(ub.df.null.nc[,1])+0.02, y=2.6)
##compare alt distributions 
lines(density(ub.df.alt.nc[,1]))
lines(density(ub.df.alt.c[,1]), col ='red')
abline(v= mean(ub.df.alt.nc[,1]), lwd=1, lty=2)
abline(v= mean(ub.df.alt.c[,1]), col='red', lwd=1, lty=2)
lines(x=c(mean(ub.df.alt.nc[,1]), mean(ub.df.alt.nc[,1])+sd(ub.df.alt.nc[,1])), col='grey', y=c(3,3), lwd=1)
lines(x=c(mean(ub.df.alt.c[,1]), mean(ub.df.alt.c[,1])+sd(ub.df.alt.c[,1])), y=c(3.1,3.1), col='pink',lwd=1)
text(round(sd(ub.df.alt.c[,1]), digits = 4), x=mean(ub.df.alt.c[,1])+0.02, y=3.4)
text(round(sd(ub.df.alt.nc[,1]), digits = 4), x=mean(ub.df.alt.nc[,1])+0.02, y=2.6)
##compare z values 
abline(v= p.505.nc[[4]][1], col='grey', lwd=0.5)
abline(v= p.505.c[[4]][1], col='pink', lwd=0.5)

####################################################
########Using NI margin as threshold ###############
####################################################

####Non confounding ####
nc.power <- function (n, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, NImargin){
  
  nIterations = 3000
  
  ###NULL HYPOTHESIS to get z value 
  id = 1 : (2*n) #create participant id  
  
  ###ALTERNATIVE HYPOTHESIS for power 
  means.alt.nc = sds.alt.nc = c()
  
  for (i in 1:nIterations) {
    confounder = rbeta(n = 2 * n, shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1 
    randomisation = sample(rep(0:1, n)) #randomisation
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome0 = rbinom(2*n, prob=p.stdcare, size = 1)
    outcome1 = rbinom(2*n, prob=p.experiment, size = 1) 
    
    #INTERVENTION dependent on adherence
    intervention = rep(NA,2*n)
    intervention[sample(which(randomisation==1), size = adhere.experiment * n)] = 1
    intervention[which(randomisation==1) %in% which(is.na(intervention))] = 0
    intervention[sample(which(randomisation==0), size = adhere.stdcare * n)] = 0
    intervention[which(randomisation==0) %in% which(is.na(intervention))] = 1
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=outcome0, outcome1=outcome1, intervention=intervention)
    simdata = matrix(data=c(id, randomisation, confounder, intervention, outcome), nrow=(2*n))
    
    ub = analysis.ub(simdata = simdata)
    means.alt.nc[[i]] = ub[[1]]
    sds.alt.nc[[i]] = ub[[2]]
  }
  means.alt.df.nc = matrix(unlist(means.alt.nc), ncol = 4, byrow = T)
  sds.df.alt.nc = matrix(unlist(sds.alt.nc), ncol = 4, byrow = T)
  ub.df.alt.nc = means.alt.df.nc + qnorm(0.975) * sds.df.alt.nc
  
  p.itt.nc = mean(ub.df.alt.nc[,1] < NImargin)
  p.pp.nc = mean(ub.df.alt.nc[,2] < NImargin)
  p.mpp.nc = mean(ub.df.alt.nc[,3] < NImargin)
  p.iv.nc = mean(ub.df.alt.nc[,4] < NImargin)
  
  return (list(c(p.itt.nc, p.pp.nc, p.mpp.nc, p.iv.nc),
               ub.df.alt.nc))
}
c.power <- function (n, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, confounder.intervention, confounder.outcome, NImargin){
  
  nIterations = 3000 
  
  id = 1: (2*n) #create participant id  
  
  ###ALTERNATIVE HYPOTHESIS for power 
  means.alt.c = sds.alt.c = c()
  
  shape2.outcome = runif(1)
  shape2.intervention = runif(1, min=2, max=10)
  
  for (i in 1:nIterations) {
    confounder = rbeta(n = 2 * n, shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1 
    randomisation = sample(rep(0:1, n)) #randomisation
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    if (confounder.outcome=="Increase likelihood") {
      #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
      shape1 =  shape2.outcome*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.stdcare, in increasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment) 
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome)) #individual probability with mean of p.experiment, in increasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome
    } else {
      shape1 =  shape2.outcome*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.stdcare, in decreasing order
      outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
      
      shape1 =  shape2.outcome*p.experiment/(1-p.experiment)
      p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2.outcome), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    }
    
    d = matrix(data = c(id, randomisation, confounder), ncol = 3)
    d.ordered = matrix(cbind(d[order(d[,3]), ], outcome0, outcome1), ncol = 5) #order confounder in ascending order 
    d.grouped = rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),]) #group according to randomisation 
    
    #INTERVENTION dependent on randomisation and confounders
    if (confounder.intervention=="Increase likelihood") {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.experiment, in increasing order
      int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
      
      shape1 =  shape2.intervention*(1-adhere.stdcare)/(1-(1-adhere.stdcare)) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in increasing order
      cont.intervention =  rbinom(n, 1, prob=adhere.stdcare.ind)
    } else {
      shape1 =  shape2.intervention*adhere.experiment/(1-adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
      int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
      
      shape1 =  shape2.intervention*adhere.stdcare/(1-adhere.stdcare) 
      adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2.intervention)) #individual probability with mean of p.stdcare, in decreasing order
      cont.intervention =  rbinom(n, 1, prob = 1-adhere.stdcare.ind)
    }
    
    intervention = c(int.intervention,cont.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0=d.grouped[,4], outcome1=d.grouped[,5], intervention)
    simdata = matrix(cbind(d.grouped[,-5:-4], intervention, outcome), ncol=5)
    
    ub = analysis.ub(simdata = simdata)
    means.alt.c[[i]] = ub[[1]]
    sds.alt.c[[i]] = ub[[2]]
  }
  means.df.alt.c = matrix(unlist(means.alt.c), ncol = 4, byrow = T)
  sds.df.alt.c = matrix(unlist(sds.alt.c), ncol = 4, byrow = T)
  ub.df.alt.c = means.df.alt.c + qnorm(0.975) * sds.df.alt.c
  
  p.itt.alt.c = mean(ub.df.alt.c[,1] < NImargin)
  p.pp.alt.c = mean(ub.df.alt.c[,2] < NImargin)
  p.mpp.alt.c = mean(ub.df.alt.c[,3] < NImargin)
  p.iv.alt.c = mean(ub.df.alt.c[,4] < NImargin)
  
  return (list(c(p.itt.alt.c, p.pp.alt.c, p.mpp.alt.c, p.iv.alt.c),  
               ub.df.alt.c))
}

p.505.nc = nc.power (n=505, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, NImargin)
p.505.c = c.power (n=505, p.experiment, p.stdcare,adehere.experiment, adhere.stdcare, confounder.intervention, confounder.outcome, NImargin)

##compare null distribution 
p.505.nc[[1]]  
p.505.c[[1]]

ub.df.alt.nc = p.505.nc [[2]]
ub.df.alt.c = p.505.c[[2]]
##compare alt distributions 
plot(density(ub.df.alt.nc[,1]))
lines(density(ub.df.alt.c[,1]), col ='red')
abline(v= mean(ub.df.alt.nc[,1]), lwd=1, lty=2)
abline(v= mean(ub.df.alt.c[,1]), col='red', lwd=1, lty=2)
lines(x=c(mean(ub.df.alt.nc[,1]), mean(ub.df.alt.nc[,1])+sd(ub.df.alt.nc[,1])), col='grey', y=c(3,3), lwd=1)
lines(x=c(mean(ub.df.alt.c[,1]), mean(ub.df.alt.c[,1])+sd(ub.df.alt.c[,1])), y=c(3.1,3.1), col='pink',lwd=1)
text(round(sd(ub.df.alt.c[,1]), digits = 4), x=mean(ub.df.alt.c[,1])+0.02, y=3.4)
text(round(sd(ub.df.alt.nc[,1]), digits = 4), x=mean(ub.df.alt.nc[,1])+0.02, y=2.6)
abline(v= NImargin, col='blue', lwd=1, lty=2)
