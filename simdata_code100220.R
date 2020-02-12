################################################################################################################
###################Using causal inference to address non-adherence in non inferiority trials####################
#############################################Simulating data functions##########################################
################################################################################################################

simdata.nonconfounding <- function(n, p.experiment, p.stdcare, p.alt, nonadhere.pop, cross.over, adhere.experiment, adhere.stdcare, i){
  
  pt.id = 1 : (2*n) #create participant id
  randomisation = sample(rep(0:1, n)) #randomisation
  confounder = rbeta(n = 2*n, shape1 = 2,shape2 = 2) #confounder beta distribution ranging 0-1
  
  #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
  outcome0 = rbinom(2*n, prob = p.stdcare, size = 1)#probability of outcome if intervention = 0
  outcome1 = rbinom(2*n, prob = p.experiment, size = 1) #probability of outcome if intervention = 1
  outcome2 = rbinom(2*n, prob = p.alt, size = 1) #probability of outcome if intervention = alternate to both experimental and stdcare
  
  #INTERVENTION dependent on adherence
  #INTERVENTION dependent on adherence
  intervention = rep(NA, 2*n)
  intervention[sample(which(randomisation == 1), size = adhere.experiment[i] * n)] = 1
  intervention[sample(which(randomisation == 0), size = adhere.stdcare[i] * n)] = 0
  
  if (cross.over == TRUE) { # nonadherent participants take up inferior alternative treatments
    intervention[intersect(which(randomisation == 1), which(is.na(intervention)))] = 0
    intervention[intersect(which(randomisation == 0), which(is.na(intervention)))] = 1
  } else {
    intervention[intersect(which(randomisation == 1), which(is.na(intervention)))] = 2
    intervention[intersect(which(randomisation == 0), which(is.na(intervention)))] = 2
  }
  
  #ACTUAL OUTCOMES depend on intervention
  outcome = getoutcome(outcome0 = outcome0, outcome1=outcome1, outcome2 = outcome2, intervention = intervention)
  simdata = matrix(data = c(pt.id, randomisation, confounder, intervention, outcome), nrow = 2*n)

  return(simdata)
}

simdata.confounding <- function(n, p.experiment, p.stdcare, p.alt, nonadhere.pop, confounder.outcome, confounder.intervention, adhere.experiment, adhere.stdcare, cross.over, i){
  
  pt.id = 1 : (2*n) #create participant id
  randomisation = sample(rep(0:1, n)) #randomisation
  confounder = rbeta(n = 2*n, shape1 = 2,shape2 = 2) #confounder beta distribution ranging 0-1 
  
  #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
  if (confounder.outcome=="Increase likelihood") {
    #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
    shape2 = runif(1)
    
    shape1 =  shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.stdcare, in increasing order
    outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
    
    shape1 =  shape2*p.experiment/(1-p.experiment) 
    p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome
    
    shape1 =  shape2*p.alt/(1-p.alt)
    p.alt.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.alt, in increasing order
    outcome2 =  rbinom(2*n, 1, prob=p.alt.ind) #increasing confounder value will have increasing probability for outcome
    
  } else {
    shape2 = runif(1)
    
    shape1 =  shape2*p.stdcare / (1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.stdcare.ind = sort(rbeta(n = 2*n, shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.stdcare, in decreasing order
    outcome0 =  rbinom(2*n, 1, prob = p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
    
    shape1 =  shape2*p.experiment / (1-p.experiment)
    p.experiment.ind = sort(rbeta(n = 2*n, shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
    outcome1 =  rbinom(2*n, 1, prob = p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1 =  shape2*p.alt/(1-p.alt) 
    p.alt.ind = sort(rbeta(n = 2*n, shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.alt, in decreasing order
    outcome2 =  rbinom(2*n, 1, prob = p.alt.ind) #increasing confounder value will have decreasing probability for outcome
  }
  
  d = matrix(data = c(pt.id, randomisation, confounder), ncol = 3)
  d.ordered = matrix(cbind(d[order(d[,3]), ], outcome0, outcome1, outcome2), ncol = 6) #order confounder in ascending order 
  d.grouped = rbind(d.ordered[which(d.ordered[,2] == 1),], d.ordered[which(d.ordered[,2] == 0),]) #group according to randomisation 
  
  #INTERVENTION dependent on randomisation and confounders
  if (confounder.intervention=="Increase likelihood") {
    shape2 = runif(1, min = 2, max = 10)
    shape1 =  shape2*adhere.experiment[i]/(1-adhere.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1 =  shape2*(1-adhere.stdcare[i])/(1-(1-adhere.stdcare[i])) 
    adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.stdcare, in increasing order
    cont.intervention =  rbinom(n, 1, prob=adhere.stdcare.ind)
    
    if (cross.over == FALSE) { #if not cross over, non adherent patients take up alternative treatment
      int.intervention[which(int.intervention == 0)] = 2
      cont.intervention[which(cont.intervention == 1)] = 2
    }
  } else {
    shape2 = runif(1, min = 2, max = 10)
    shape1 =  shape2*adhere.experiment[i]/(1-adhere.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.experiment.ind = sort(rbeta(n = n, shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
    int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1 =  shape2*adhere.stdcare[i]/(1-adhere.stdcare[i]) 
    adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.stdcare, in decreasing order
    cont.intervention =  rbinom(n, 1, prob = 1 - adhere.stdcare.ind)
    
    if (cross.over == FALSE) { #if not cross over, non adherent patients take up alternative treatment
      int.intervention[which(int.intervention == 0)] = 2
      cont.intervention[which(cont.intervention == 1)] = 2
    }
  }
  
  intervention = c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome = getoutcome(outcome0 = d.grouped[,4], outcome1 = d.grouped[,5], outcome2 = d.grouped[,6], intervention)
  simdata = matrix(cbind(d.grouped[,-6:-4], intervention, outcome), ncol = 5)
  
  return(simdata)
}

simdata.unknownconfounding <- function(n, p.experiment, p.stdcare, p.alt, nonadhere.pop,confounder.outcome, confounder.intervention, adhere.experiment, adhere.stdcare, cross.over, i){
  
  pt.id = 1 : (2*n) #create participant id
  randomisation = sample(rep(0:1, n)) #randomisation
  known = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
  unknown1 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  unknown2 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  unknown3 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  confounder = known + unknown1 + unknown2 + unknown3
  
  #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
  if (confounder.outcome=="Increase likelihood") {
    #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
    shape2 = runif(1)
    
    shape1 =  shape2*p.experiment/(1-p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome 
    
    shape1 =  shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
    
  } else {
    shape2 = runif(1)
    shape1 =  shape2*p.experiment/(1-p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.experiment.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
    outcome1 =  rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1 =  shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    p.stdcare.ind = sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    outcome0 =  rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome 
  }
  
  d = matrix(data=c(pt.id, randomisation, confounder, known, unknown1, unknown2, unknown3), nrow=(2*n))
  d.ordered = matrix(cbind(d[order(d[,3]),], outcome0, outcome1),ncol=9)#order confounder in ascending order 
  d.grouped = rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),])
  
  #INTERVENTION dependent on randomisation and confounders
  if (confounder.intervention=="Increase likelihood") {
    shape2 = runif(1, min=2, max=10)
    shape1 =  shape2*adhere.experiment[i]/(1-adhere.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1 =  shape2*(1-adhere.stdcare[i])/(1-(1-adhere.stdcare[i])) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention =  rbinom(n, 1, prob=adhere.stdcare.ind)
    
  } else {
    shape2 = runif(1,min=2, max=10)
    shape1 =  shape2*adhere.experiment[i]/(1-adhere.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.experiment.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    int.intervention =  rbinom(n, 1, prob=adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1 =  shape2*adhere.stdcare[i]/(1-adhere.stdcare[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    adhere.stdcare.ind = sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention =  rbinom(n, 1, prob=1-adhere.stdcare.ind)
  }
  
  intervention = c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome = getoutcome.unknownconfounding.multi(vector.outcome1=d.grouped[,9], vector.outcome0=d.grouped[,8], intervention)
  simdata = matrix(cbind(d.grouped[,c(-8,-9)],intervention,outcome), ncol = 9)
  
  return(simdata)
}
