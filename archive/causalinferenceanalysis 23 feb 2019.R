################################################################################################################
###################Using causal inference to address non-compliance in non inferiority trials###################
################################################################################################################

######## Set up #########
setwd("/Users/moyin/Desktop/angelsfly/NItrialsimulation/codes/") #set working directory 

# Required libraries 
library(Hmisc); library(rms); library(gsDesign)
library(data.table) #for melt 
library(dplyr) 
library(tableone); library(scales);
library(gtools) #for logit 
library(Matching); library(survey);library(gmm) #for analysis 
library(ggpubr); library(ggplot2); library(gridExtra);library(plotly) #for plots 

# The colour blind friendly palette 
cbPalette <- c("#009E73", "#CC79A7", "#E69F00", "#0072B2", "#D55E00", "#56B4E9")

#Analysis methods 
analysis.method<- c("As treated","Instrumental variable","Intention to treat","Inverse probability weighting","Matching","Per protocol")

#set number of data points at which simulated data is analysed  
interval <- seq(from=0.6,to=1,by=0.05)

################Source functions#################

simdata.nonconfounding<- function(n, p.experiment, p.stdcare, noncomply,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  }
  
  id=seq(1,(2*n), by=1) #create participant id  
  randomisation=c(rep(1,n), rep(0,n)) #randomisation
  confounder=rbeta(n=(2*n),shape1=2,shape2=2) #confounder beta distribution ranging 0-1 
  
  #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
  outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment,1-p.experiment))  #probability of outcome if intervention = 1
  outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))       #probability of outcome if intervention = 0
  
  #INTERVENTION dependent on compliance 
  experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1-comply.experiment[i],comply.experiment[i]))
  stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(comply.stdcare[i],1-comply.stdcare[i]))
  intervention = c(experiment.intervention,stdcare.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome<- rep(NA,2*n)
  for (y in 1:(2*n)) {
    ifelse (intervention[y]==1, outcome[y]<-outcome1[y], outcome[y]<-outcome0[y])
  }
  
  return(matrix(data=c(id,randomisation,confounder,intervention,outcome), nrow=(2*n)))
  
}

simdata.confounding<- function(n, p.experiment, p.stdcare, noncomply,confounder.outcome, confounder.intervention,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  }
  
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
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1<- shape2*(1-comply.stdcare[i])/(1-(1-comply.stdcare[i])) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)

  } else {
    shape2<-runif(1)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 

    shape1<- shape2*comply.stdcare[i]/(1-comply.stdcare[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=1-comply.stdcare.ind)
  }
  
  intervention<-c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome<-c()
  for (y in 1:(2*n)) {
    if (intervention[y]==1) {outcome[y]<-d.grouped[,5][y]} else { 
      outcome[y]<-d.grouped[,4][y]}
  }

  return(matrix(cbind(d.grouped[,c(-4,-5)],intervention,outcome), ncol=5))
}

analysis<- function(simdata, type,n,z,NImargin){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,5])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,5])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,4]),] # perprotocol population
  p.experiment.vector= pp[which(pp[,2]==1),][,5]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,5]
  p.stdcare.value= mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## as treated
  pd1.vector= simdata[which(simdata[,4]==1),][,5]
  pd1.value=mean(pd1.vector)
  pd0.vector= simdata[which(simdata[,4]==0),][,5]
  pd0.value= mean(pd0.vector)
  eff.at = pd1.value-pd0.value  
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
  ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= ifelse(pp$intervention==1,1/score, 1/(1-score)) #create weights
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  
  ## Matching 
  m<- glm(V4~V3, family = binomial(link='logit'), data=as.data.frame(simdata)) 
  mscore<- m$fitted.values
  match<-Match(Tr=simdata[,4], M=1, X= logit(mscore), replace=FALSE, caliper=0.2)
  matched<-simdata[unlist(match[c('index.treated','index.control')]),]
  m1.vector<-matched[,5][matched[,4]==1]
  m0.vector<-matched[,5][matched[,4]==0]
  eff.ps=mean(m1.vector)-mean(m0.vector)
  
  # iv with 2 stage regression
  asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
  eff.iv=(summary(asmm))$ coefficients [2,1]
  
  if (type=="bias"){
    
    return(c(eff.at, eff.iv, eff.itt, eff.mpp, eff.ps, eff.pp))
    
  } else { #type 1 error and power 
    
    ## intention to treat 
    var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
    CI.itt<- eff.itt + z*sqrt(var.eff.itt)
    itt<-CI.itt<NImargin
    
    ## per protocol 
    var.eff.pp<-  p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
    CI.pp<- eff.pp + z*sqrt(var.eff.pp)
    ppp<- CI.pp<NImargin
    
    ## as treated
    var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
    CI.at<- eff.at + z*sqrt(var.eff.at)
    at<- CI.at<NImargin
    
    ## inverse probability weights
    se<-sqrt(diag(vcovHC(outcomemodel,type="HC0")))[2]
    CI.mpp<-eff.mpp+z*se
    mpp<- CI.mpp<NImargin
    
    ## Matching 
    var.eff.ps<-mean(m1.vector)*(1-mean(m1.vector))/length(m1.vector) + mean(m0.vector)*(1-mean(m0.vector))/length(m0.vector)
    CI.ps<- eff.ps + z*sqrt(var.eff.ps)
    ps<- CI.ps<NImargin
    
    # iv with 2 stage regression
    se<-(summary(asmm))$coefficients [2,2]
    CI.iv<-eff.iv + z*se
    iv<- CI.iv<NImargin
    
    return(c(at, iv, itt, mpp, ps, ppp))
  }
}

plot.eff<- function(df,method,nIterations){
  
  k<-which(method==analysis.method)
  x<-c()
  y<-c()
  
  for (i in 1:length(interval)) {
    x[[i]]<-df[[i]][,k]
    y[[i]]<-rep(interval[i],dim(df[[i]])[1])
  }
  
  plotdata<-matrix(c(unlist(x),unlist(y)), ncol = 2)
  
  plot<-ggplot(as.data.frame(plotdata), aes(x=V2, y=V1))+ 
    geom_point(aes(y=V1, colour=method, alpha=0.01), size=0.2) + 
    geom_smooth(aes(y=V1, colour=method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    xlab("Proportion of compliant participants")+
    ylab("Effect")+
    theme_bw()+
    scale_colour_manual(values=cbPalette[k])+
    theme(legend.title=element_blank(), legend.position="none")+
    scale_x_continuous(limits=c(0.60, 1))+
    scale_y_continuous(limits=c(-0.4, 0.2))
  
  return(plot)
}

#################CASE 1 Non-compliance caused by non confounding process#########################
#################################################################################################

#BIAS
bias.nonconfounding<- function(n, p.experiment, p.stdcare, nIterations, interval, noncomply){  
  
  #make up vectors for simulations 
  .estimate<-c() #output from each simulation
  estimate<-c()  #for saving output from each interval
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata<-simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i)
        .estimate[[l]]<-analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]]<-matrix(unlist(.estimate), ncol = 6, byrow = TRUE)
  }
  
  #plot
  itt<- plot.eff(df=estimate,analysis.method[1],nIterations=nIterations)
  pp<- plot.eff(df=estimate,analysis.method[2],nIterations=nIterations)
  at<- plot.eff(df=estimate,analysis.method[3],nIterations=nIterations)
  mpp<- plot.eff(df=estimate,analysis.method[4],nIterations=nIterations)
  ps<- plot.eff(df=estimate,analysis.method[5],nIterations=nIterations)
  iv<- plot.eff(df=estimate,analysis.method[6],nIterations=nIterations)
  
  bias.plot<-ggarrange(itt, pp, at,
                       iv, mpp, ps,
                       ncol = 3, nrow = 2,
                       common.legend = FALSE )
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.nonconfounding<- function(n, p.experiment, nIterations, interval,NImargin, noncomply){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate<-c() #output from each simulation
  estimate<-c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata<-simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i)
        .estimate[[l]]<-analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    # mean of type 1 error from iterated data 
    estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = 6, byrow = TRUE))
  }
  
  df<-as.data.frame(matrix(unlist(estimate),ncol=6, byrow=TRUE))
  
  plot <- ggplot(df, aes(interval))+ 
    geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.01)) + 
    geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.01)) + 
    geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.01)) + 
    geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.01)) + 
    geom_point(aes(y=df$V5, colour=analysis.method[5], alpha=0.01)) +
    geom_point(aes(y=df$V6, colour=analysis.method[6], alpha=0.01)) + 
    geom_smooth(aes(y=df$V1, colour=analysis.method[1]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V2, colour=analysis.method[2]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V3, colour=analysis.method[3]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V4, colour=analysis.method[4]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V5, colour=analysis.method[5]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V6, colour=analysis.method[6]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    xlab("Proportion of compliant participants")+
    ylab("Type 1 error")+
    theme_bw()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    scale_x_continuous(limits=c(0.6, 1))+
    scale_y_continuous(limits=c(0, 0.9),breaks=seq(0,0.9,0.1))
  
  return(plot)
  
}

#POWER
power.nonconfounding<- function(n, p.experiment, p.stdcare, NImargin,interval,nIterations,noncomply){  
  
  #make up vectors for simulations 
  .estimate<-c() #output from each simulation
  estimate<-c()  #for saving output from each interval
  
  #alpha error and critical value 
  z <- qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare 
  
  if ((NImargin>0) & (true.eff < NImargin)) { #built with alternative hypothesis: true effect < NI 
    
    #simulate and derive treatment effect 
    for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
      for(l in 1:nIterations) { 
        tryCatch({
          simdata<-simdata.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i)
          .estimate[[l]]<-analysis(simdata=simdata, type='power',n=n,z=z,NImargin=NImargin)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
        
      }
      # mean of power from iterated data 
      estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = 6, byrow = TRUE))
    }
    
    df<-as.data.frame(matrix(unlist(estimate),ncol=6, byrow=TRUE))
    
    plot <- ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.01)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.01)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.01)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.01)) + 
      geom_point(aes(y=df$V5, colour=analysis.method[5], alpha=0.01)) +
      geom_point(aes(y=df$V6, colour=analysis.method[6], alpha=0.01)) + 
      geom_smooth(aes(y=df$V1, colour=analysis.method[1]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V2, colour=analysis.method[2]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V3, colour=analysis.method[3]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V4, colour=analysis.method[4]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V5, colour=analysis.method[5]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V6, colour=analysis.method[6]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      xlab("Proportion of compliant participants")+
      ylab("Power")+
      theme_bw()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom")+
      theme(legend.title=element_blank())+
      scale_x_continuous(limits=c(0.6, 1))+
      scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1))
    
    return(plot)
    
  }
  else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

############################CASE 2 Non-compliance caused by confounding process  ################################
################################################################################################################
#BIAS
bias.confounding<- function(n, p.experiment, p.stdcare, confounder.intervention, confounder.outcome,interval,nIterations,noncomply){  
  
  #make up vectors for simulations 
  estimate<-c()  #for saving output from each interval
  .estimate<-c() #output from each simulation
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata<-simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]]<-matrix(unlist(.estimate), ncol = 6, byrow = TRUE)
  }
  
  #plot
  itt<- plot.eff(df=estimate,analysis.method[1],nIterations=nIterations)
  pp<- plot.eff(df=estimate,analysis.method[2],nIterations=nIterations)
  at<- plot.eff(df=estimate,analysis.method[3],nIterations=nIterations)
  mpp<- plot.eff(df=estimate,analysis.method[4],nIterations=nIterations)
  ps<- plot.eff(df=estimate,analysis.method[5],nIterations=nIterations)
  iv<- plot.eff(df=estimate,analysis.method[6],nIterations=nIterations)
  
  bias.plot<-ggarrange(itt, pp, at,
                       iv, mpp, ps,
                       ncol = 3, nrow = 2,
                       common.legend = FALSE )
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.confounding<- function(n, p.experiment, NImargin, confounder.intervention, confounder.outcome,interval,nIterations,noncomply){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff=p.experiment-p.stdcare
  
  #make up vectors for simulations 
  .estimate<-c() #output from each simulation
  estimate<-c()  #for saving output from each interval
  
  #simulate data
  for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
    for(l in 1:nIterations) {  
      tryCatch({
        simdata<-simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i,confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysis(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of type 1 error from iterated data 
    estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = 6, byrow = TRUE),na.rm = TRUE)
  }
  
  df<-as.data.frame(matrix(unlist(estimate),ncol=6, byrow=TRUE))
  
  plot <- ggplot(df, aes(interval))+ 
    geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.01)) + 
    geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.01)) + 
    geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.01)) + 
    geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.01)) + 
    geom_point(aes(y=df$V5, colour=analysis.method[5], alpha=0.01)) +
    geom_point(aes(y=df$V6, colour=analysis.method[6], alpha=0.01)) + 
    geom_smooth(aes(y=df$V1, colour=analysis.method[1]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V2, colour=analysis.method[2]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V3, colour=analysis.method[3]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V4, colour=analysis.method[4]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V5, colour=analysis.method[5]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=df$V6, colour=analysis.method[6]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    xlab("Proportion of compliant participants")+
    ylab("Type 1 error")+
    theme_bw()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    scale_x_continuous(limits=c(0.6, 1))+
    scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1))
  
  return(plot)
  
}

#POWER
power.confounding<- function(n, p.experiment, p.stdcare, NImargin, confounder.intervention, confounder.outcome,interval,nIterations,noncomply){  
  
  #make up vectors for simulations 
  .estimate<-c() #output from each simulation
  estimate<-c()  #for saving output from each interval
  
  #alpha error and critical value 
  z <- qnorm(0.975) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff<-p.experiment-p.stdcare 
  
  if ((NImargin>0) & (true.eff < NImargin)) { #built with alternative hypothesis: true effect < NI 
    
    #simulate and derive treatment effect 
    for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
      for(l in 1:nIterations) { 
        tryCatch({
          simdata<-simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i,confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
          .estimate[[l]]<-analysis(simdata=simdata, type='power',n=n,z=z,NImargin=NImargin)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
        
      }
      # mean of power from iterated data 
      estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = 6, byrow = TRUE),na.rm = TRUE)
    }
    
    df<-as.data.frame(matrix(unlist(estimate),ncol=6, byrow=TRUE))
    
    plot <- ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.01)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.01)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.01)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.01)) + 
      geom_point(aes(y=df$V5, colour=analysis.method[5], alpha=0.01)) +
      geom_point(aes(y=df$V6, colour=analysis.method[6], alpha=0.01)) + 
      geom_smooth(aes(y=df$V1, colour=analysis.method[1]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V2, colour=analysis.method[2]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V3, colour=analysis.method[3]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V4, colour=analysis.method[4]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V5, colour=analysis.method[5]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=df$V6, colour=analysis.method[6]), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      xlab("Proportion of compliant participants")+
      ylab("Power")+
      theme_bw()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom")+
      theme(legend.title=element_blank())+
      scale_x_continuous(limits=c(0.6, 1))+
      scale_y_continuous(limits=c(0.2, 1),breaks=seq(0.2,1,0.1))
    
    return(plot)
    
  }
  else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

