rm(list=ls())
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/") #set working directory 
source(file = 'simulationcode_110319.R')
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/shiny/samplesize_nonadherence/")
sourceCpp(file='rcpp.cpp')

ncpower<- function(n, significance, p.experiment, p.stdcare, comply.experiment, comply.stdcare, NImargin){  
  options(digits=2)
  
  #make up vectors for simulations 
  power.iter<-c() 
  
  #number of iterations 
  nIterations=1000 
  
  #significance 
  ifelse (significance=="1 sided 97.5%", z <- qnorm(0.975), z <- qnorm(0.95))
  
  if ((p.experiment-p.stdcare) >= NImargin) stop ("Error: NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation") #built with alternative hypothesis: true effect < NI 
  
  #simulate and derive treatment effect 
  for(l in 1:nIterations) { 
    tryCatch({
      
      id=seq(1,(2*n), by=1) #create participant id  
      randomisation=c(rep(1,n), rep(0,n)) #randomisation
      confounder=rbeta(n=(2*n),shape1=2,shape2=2) #confounder beta distribution ranging 0-1 
      
      #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
      outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment,1-p.experiment))  #probability of outcome if intervention = 1
      outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))       #probability of outcome if intervention = 0
      
      #INTERVENTION dependent on compliance 
      experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1-comply.experiment,comply.experiment))
      stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(comply.stdcare,1-comply.stdcare))
      intervention = c(experiment.intervention,stdcare.intervention)
      
      #ACTUAL OUTCOMES depend on intervention
      outcome<-getoutcome(outcome1, outcome0, intervention)
      
      simdata<-matrix(data=c(id,randomisation,confounder,intervention,outcome), nrow=(2*n))
      
      ## intention to treat 
      pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
      pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
      eff.itt = pz1.value-pz0.value
      var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
      CI.itt<- eff.itt + z*sqrt(var.eff.itt)
      itt<-CI.itt<NImargin
      
      ## per protocol 
      pp = simdata[which(simdata[,2]==simdata[,4]),] # perprotocol population
      p.experiment.vector= pp[which(pp[,2]==1),][,5]
      p.experiment.value= mean(p.experiment.vector)
      p.stdcare.vector= pp[which(pp[,2]==0),][,5]
      p.stdcare.value= mean(p.stdcare.vector) 
      eff.pp = p.experiment.value-p.stdcare.value
      var.eff.pp<-  p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
      CI.pp<- eff.pp + z*sqrt(var.eff.pp)
      ppp<- CI.pp<NImargin
      
      ## inverse probability weights on per protocol patients 
      pp=as.data.frame(pp)
      colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
      ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
      score=predict(ipwmodel, type="response")
      weight= ifelse(pp$intervention==1,1/score, 1/(1-score)) #create weights
      outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
      eff.mpp=coef(outcomemodel)[2]
      se<-squrt(diag(vcovHC(outcomemodel)))[2]
      CI.mpp<-eff.mpp+z*se
      mpp<- CI.mpp<NImargin

      # iv with 2 stage regression
      asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
      eff.iv=(summary(asmm))$ coefficients [2,1]
      se<-(summary(asmm))$coefficients [2,2]
      CI.iv<-eff.iv + z*se
      iv<- CI.iv<NImargin
      
      power.iter[[l]]<-c(iv, itt, mpp, ppp)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
  }
  
  # mean of power from iterated data 
  power.matrix=matrix(as.numeric(unlist(power.iter)), nrow=nIterations, ncol=4,byrow = TRUE)
  power.cal= colMeans(power.matrix, na.rm = TRUE)
  return(power.cal)
}
cpower<- function(n, p.experiment, p.stdcare, comply.experiment, comply.stdcare, NImargin, significance, confounder.intervention, confounder.outcome){  
  options(digits=2)
  #make up vectors for simulations 
  power.iter<-c() 
  
  #number of iterations 
  nIterations=1000
  
  #significance 
  ifelse (significance=="1 sided 97.5%", z <- qnorm(0.975), z <- qnorm(0.95))
  
  if ((p.experiment-p.stdcare) > NImargin) stop ("Error: NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation") #built with alternative hypothesis: true effect < NI 
  
  #simulate and derive treatment effect 
  for(l in 1:nIterations) { 
    tryCatch({
      id = seq(1,(2*n), by=1) #create participant id  
      randomisation = c(rep(1,n), rep(0,n)) #randomisation
      confounder = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
      
      #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
      if (confounder.outcome=="Increase likelihood") {
        #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome 
        shape2<-runif(1)
        shape1<- shape2*p.experiment/(1-p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        p.experiment.ind<-sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
        outcome1<- rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome 
        
        shape1<- shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        p.stdcare.ind<-sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
        outcome0<- rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have increasing probability for outcome 
      } else {
        shape2<-runif(1)
        shape1<- shape2*p.experiment/(1-p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        p.experiment.ind<-sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in decreasing order
        outcome1<- rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
        
        shape1<- shape2*p.stdcare/(1-p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        p.stdcare.ind<-sort(rbeta(n=(2*n), shape1 = shape1, shape2 = shape2), decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
        outcome0<- rbinom(2*n, 1, prob=p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome 
      }
      
      d<-matrix(data=c(id, randomisation, confounder), nrow=(2*n))
      d.ordered<-matrix(cbind(d[order(d[,3]),], outcome0, outcome1),ncol=5)#order confounder in ascending order 
      d.grouped<-rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),])
      
      #INTERVENTION dependent on randomisation and confounders
      if (confounder.intervention=="Increase likelihood") {
        shape2<-runif(1)
        shape1<- shape2*comply.experiment/(1-comply.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
        int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
        
        shape1<- shape2*(1-comply.stdcare)/(1-(1-comply.stdcare)) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
        cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)
        
      } else {
        shape2<-runif(1)
        shape1<- shape2*comply.experiment/(1-comply.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
        int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
        
        shape1<- shape2*comply.stdcare/(1-comply.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
        comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
        cont.intervention<- rbinom(n, 1, prob=1-comply.stdcare.ind)
      }
      
      intervention<-c(int.intervention,cont.intervention)
      
      #ACTUAL OUTCOMES depend on intervention
      outcome<-getoutcome(d.grouped[,5], d.grouped[,4], intervention)
      
      simdata<-matrix(cbind(d.grouped[,c(-4,-5)],intervention,outcome), ncol=5)
      
      ## intention to treat 
      pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
      pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
      eff.itt = pz1.value-pz0.value
      var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
      CI.itt<- eff.itt + z*sqrt(var.eff.itt)
      itt<-CI.itt<NImargin
      
      ## per protocol 
      pp = simdata[which(simdata[,2]==simdata[,4]),] # perprotocol population
      p.experiment.vector= pp[which(pp[,2]==1),][,5]
      p.experiment.value= mean(p.experiment.vector)
      p.stdcare.vector= pp[which(pp[,2]==0),][,5]
      p.stdcare.value= mean(p.stdcare.vector) 
      eff.pp = p.experiment.value-p.stdcare.value
      var.eff.pp<-  p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
      CI.pp<- eff.pp + z*sqrt(var.eff.pp)
      ppp<- CI.pp<NImargin

      ## inverse probability weights on per protocol patients 
      pp=as.data.frame(pp)
      colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
      ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
      score=predict(ipwmodel, type="response")
      weight= ifelse(pp$intervention==1,1/score, 1/(1-score)) #create weights
      outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
      eff.mpp=coef(outcomemodel)[2]
      se<-squrt(diag(vcovHC(outcomemodel)))[2]
      CI.mpp<-eff.mpp+z*se
      mpp<- CI.mpp<NImargin
      
      # iv with 2 stage regression
      asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
      eff.iv=(summary(asmm))$ coefficients [2,1]
      se<-(summary(asmm))$coefficients [2,2]
      CI.iv<-eff.iv + z*se
      iv<- CI.iv<NImargin
      
      power.iter[[l]]<-c(iv, itt, mpp,ppp)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
  }
  
  # mean of power from iterated data 
  power.matrix=matrix(as.numeric(unlist(power.iter)), nrow=nIterations, ncol=4,byrow = TRUE)
  power.cal= colMeans(power.matrix, na.rm = TRUE)
  return(power.cal)
  
}

ncpower(n=330, significance = "1 sided 97.5%", p.experiment = 0.3, p.stdcare = 0.3, comply.experiment = 1, comply.stdcare = 1, NImargin = 0.1)
cpower(n=330, significance = "1 sided 97.5%", p.experiment = 0.3, p.stdcare = 0.3, comply.experiment = 1, comply.stdcare = 1, NImargin = 0.1, 
       confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood")

setwd('/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/shiny/samplesize_nonadherence/')
require('rsconnect')
rsconnect::setAccountInfo(name='moru', 
                          token='30819BAEDD492333CE7CD293F3B08D42', 
                          secret='j24TwHndsiybKqQZERFOt0s3L3ro8IG47/OgpLk/')
deployApp(account="moru",appName="samplesize_nonadherence")
