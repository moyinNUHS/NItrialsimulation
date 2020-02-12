################################################################################################################
###################Using causal inference to address non-adherence in non inferiority trials####################
#############################################Analysis functions#################################################
################################################################################################################

analysis.estimate <- function(simdata){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,4]),] # per protocol population
  p.experiment.vector = pp[which(pp[,2]==1), ][,5]
  p.experiment.value = mean(p.experiment.vector)
  p.stdcare.vector = pp[which(pp[,2]==0),][,5]
  p.stdcare.value = mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## inverse probability weights on per protocol patients 
  pp = as.data.frame(pp)
  colnames(pp) = c('id','randomisation','confounder','intervention','outcome')
  ipwmodel = glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score = predict(ipwmodel, type="response")
  weight = pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel = glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp = coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid") 
  eff.iv = summary(asmm)$coefficients[2,1] 
  
  return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
}

analysis.unknown <- function(simdata){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,7])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,7])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,6]),] # perprotocol population
  p.experiment.vector = pp[which(pp[,2]==1),][,7]
  p.experiment.value = mean(p.experiment.vector)
  p.stdcare.vector = pp[which(pp[,2]==0),][,7]
  p.stdcare.value = mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder', 'known', 'unknown1','unknown2','unknown3','intervention','outcome')
  
  ipwmodel.known=glm(intervention~known,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score.known=predict(ipwmodel.known, type="response")
  weight.known= pp$intervention*mean(pp$intervention)/score.known+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score.known)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel.known=glm(outcome~intervention, family=binomial(link="identity"), weights=weight.known, data=pp) #identity link for risk difference
  eff.mpp.known=coef(outcomemodel.known)[2]
  
  ipwmodel.unknown1=glm(intervention~known+unknown1,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score.unknown1=predict(ipwmodel.unknown1, type="response")
  weight.unknown1= pp$intervention*mean(pp$intervention)/score.unknown1+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score.unknown1)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel.unknown1=glm(outcome~intervention, family=binomial(link="identity"), weights=weight.unknown1, data=pp) #identity link for risk difference
  eff.mpp.unknown1=coef(outcomemodel.unknown1)[2]
  
  ipwmodel.unknown2=glm(intervention~known+unknown1+unknown2,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score.unknown2=predict(ipwmodel.unknown2, type="response")
  weight.unknown2= pp$intervention*mean(pp$intervention)/score.unknown2+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score.unknown2)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel.unknown2=glm(outcome~intervention, family=binomial(link="identity"), weights=weight.unknown2, data=pp) #identity link for risk difference
  eff.mpp.unknown2=coef(outcomemodel.unknown2)[2]
  
  ipwmodel.all=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score.all=predict(ipwmodel.all, type="response")
  weight.all= pp$intervention*mean(pp$intervention)/score.all+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score.all)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel.all=glm(outcome~intervention, family=binomial(link="identity"), weights=weight.all, data=pp) #identity link for risk difference
  eff.mpp.all=coef(outcomemodel.all)[2]
  
  # iv with 2 stage regression
  asmm  =  gmm(simdata[,7] ~ simdata[,6], x=simdata[,2], vcov="iid")
  eff.iv = summary(asmm)$coefficients [2,1]
  
  return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
}

sim.analysis <- function(nonconfounding, bias, nonadhere.pop, interval, cross.over, 
                         confounder.intervention, confounder.outcome, NImargin){
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  #non-adherent population 
  if (nonadhere.pop == "both") {
    adhere.experiment = interval
    adhere.stdcare = interval 
  } else if (nonadhere.pop == "experimental") {
    adhere.experiment =  interval
    adhere.stdcare =  rep(1, length(interval))
  } else { #nonadhere.pop=="stdcare"
    adhere.experiment = rep(1, length(interval))
    adhere.stdcare = interval
  }
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste("interval value", i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        if (nonconfounding == 'nonconfounding') { #simulate data 
          simdata = simdata.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, nonadhere.pop = nonadhere.pop, adhere.experiment = adhere.experiment, adhere.stdcare = adhere.stdcare, cross.over = cross.over, i=i)
        } else  if (nonconfounding == 'confounding' ) {
          simdata = simdata.confounding(n = n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop = nonadhere.pop, adhere.experiment = adhere.experiment, adhere.stdcare = adhere.stdcare, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention, i=i)
        } else {
          simdata = simdata.unknownconfounding(n = n, p.experiment=p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop = nonadhere.pop, i=i, adhere.experiment = adhere.experiment, adhere.stdcare = adhere.stdcare, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        }
        .estimate[[l]] = analysis.estimate(simdata = simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp for each iteration
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  if (bias == T) { #return the estimate matrix 
    
    return (estimate)
    
  } else { #go on to get upper bounds for type 1 error 
    
    t1 = c()
    #get type 1 error for each interval 
    for (i in 1:length(interval)) {
      #get sds for each interval 
      sd.vector = unlist(apply(estimate[[i]], 2, sd))
      #get upper bounds of 95%CI for each iteration of every interval
      ub = estimate[[i]] + qnorm(0.975) * rep(sd.vector, each = nIterations)
      #get critical values for each interval
      z.vector = unlist(apply(ub, 2, function (x) quantile(x, probs = 0.025)))
      
      t1[[i]] = colMeans(ub < NImargin)
    }
    
    t1.df = as.data.frame(matrix(unlist(t1), ncol=4, byrow=TRUE))
    
    return (t1.df)
  }
  
}