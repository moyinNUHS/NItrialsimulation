################################################################################################################
###################Using causal inference to address non-compliance in non inferiority trials###################
################################################################################################################

######## Set up #########
# Required libraries 
library(Rcpp)
library(Hmisc); library(rms); library(gsDesign)
library(dplyr); library(tableone); library(scales);
library(survey);library(gmm) #for analysis 
library(ggpubr); library(ggplot2); library(gridExtra);library(plotly) #for plots 

# The colour blind friendly palette 
cbPalette <- c("#009E73", "#CC79A7", "#E69F00", "#0072B2")
analysis.method<- c("Instrumental variable","Intention to treat","Inverse probability weighting","Per protocol")

#set number of data points at which simulated data is analysed  
start.interval= 0.6
interval <- seq(from=start.interval, to=1,by=0.025)

################Source functions#################

getoutcome<-function(vector.outcome1, vector.outcome0, intervention){
  outcome<-c()
  for (i in 1:length(vector.outcome0)){
    if (intervention[i]==1) {outcome[[i]]=vector.outcome1[i]} else {outcome[[i]]=vector.outcome0[i]}
  }
  return(unlist(outcome))
}

simdata.nonconfounding<- function(n, p.experiment, p.stdcare, noncomply,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- interval
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-interval
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
  outcome<-getoutcome(vector.outcome1=outcome1, vector.outcome0=outcome0, intervention=intervention)
  
  return(matrix(data=c(id,randomisation,confounder,intervention,outcome), nrow=(2*n)))
}

simdata.confounding<- function(n, p.experiment, p.stdcare, noncomply,confounder.outcome, confounder.intervention,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- interval
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-interval
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
    shape2<-runif(1, min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1<- shape2*(1-comply.stdcare[i])/(1-(1-comply.stdcare[i])) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)
    
  } else {
    shape2<-runif(1,min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1<- shape2*comply.stdcare[i]/(1-comply.stdcare[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=1-comply.stdcare.ind)
  }
  
  intervention<-c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome<-getoutcome(vector.outcome1=d.grouped[,5], vector.outcome0=d.grouped[,4], intervention)
  
  return(matrix(cbind(d.grouped[,c(-4,-5)],intervention,outcome), ncol=5))
}

simdata.unknownconfounding<- function(n, p.experiment, p.stdcare, noncomply,confounder.outcome, confounder.intervention,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- interval
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-interval
  }
  
  id = seq(1,(2*n), by=1) #create participant id  
  randomisation = c(rep(1,n), rep(0,n)) #randomisation
  known = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
  unknown = rep(rbeta(n=n,shape1=2,shape2=2),2)
  confounder=known*unknown
  
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
  
  d<-matrix(data=c(id, randomisation, confounder, known, unknown), nrow=(2*n))
  d.ordered<-matrix(cbind(d[order(d[,3]),], outcome0, outcome1),ncol=7)#order confounder in ascending order 
  d.grouped<-rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),])
  
  #INTERVENTION dependent on randomisation and confounders
  if (confounder.intervention=="Increase likelihood") {
    shape2<-runif(1, min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1<- shape2*(1-comply.stdcare[i])/(1-(1-comply.stdcare[i])) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)
    
  } else {
    shape2<-runif(1,min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1<- shape2*comply.stdcare[i]/(1-comply.stdcare[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=1-comply.stdcare.ind)
  }
  
  intervention<-c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome<-getoutcome(vector.outcome1=d.grouped[,7], vector.outcome0=d.grouped[,6], intervention)
  
  return(matrix(cbind(d.grouped[,c(-6,-7)],intervention,outcome), ncol=7))
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
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder','intervention','outcome')
  ipwmodel=glm(intervention~confounder,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
  eff.iv=summary(asmm)$coefficients [2,1]
  
  if (type=="bias"){
    
    return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
    
  } else { #type 1 error and power 
    
    ## intention to treat 
    var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
    CI.itt<- eff.itt + z*sqrt(var.eff.itt)
    itt<-CI.itt<NImargin
    
    ## per protocol 
    var.eff.pp<-  p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
    CI.pp<- eff.pp + z*sqrt(var.eff.pp)
    ppp<- CI.pp<NImargin
    
    ## inverse probability weights
    se<-summary(outcomemodel)$coefficients[2,2]
    CI.mpp<-eff.mpp+z*se
    mpp<- CI.mpp<NImargin
    
    # iv with 2 stage regression
    se<-summary(asmm)$coefficients [2,2]
    CI.iv<-eff.iv + z*se
    iv<- CI.iv<NImargin
    
    return(c(iv, itt, mpp, ppp))
  }
}

analysisunknown<- function(simdata, type,n,z,NImargin){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,7])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,7])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,6]),] # perprotocol population
  p.experiment.vector= pp[which(pp[,2]==1),][,7]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,7]
  p.stdcare.value= mean(p.stdcare.vector) 
  eff.pp = p.experiment.value-p.stdcare.value
  
  ## inverse probability weights on per protocol patients 
  pp=as.data.frame(pp)
  colnames(pp)=c('id','randomisation','confounder', 'known', 'unknown','intervention','outcome')
  ipwmodel=glm(intervention~known,family=binomial(link="logit"), data=pp) #calculate denominators used in inverse probability weights
  score=predict(ipwmodel, type="response")
  weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
  eff.mpp=coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm <- gmm(simdata[,7] ~ simdata[,6], x=simdata[,2], vcov="iid")
  eff.iv=summary(asmm)$coefficients [2,1]
  
  if (type=="bias"){
    
    return(c(eff.iv, eff.itt, eff.mpp, eff.pp))
    
  } else { #type 1 error and power 
    
    ## intention to treat 
    var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
    CI.itt<- eff.itt + z*sqrt(var.eff.itt)
    itt<-CI.itt<NImargin
    
    ## per protocol 
    var.eff.pp<-  p.experiment.value*(1-p.experiment.value)/length(p.experiment.vector) + p.stdcare.value*(1-p.stdcare.value)/length(p.stdcare.vector)
    CI.pp<- eff.pp + z*sqrt(var.eff.pp)
    ppp<- CI.pp<NImargin
    
    ## inverse probability weights
    se<-summary(outcomemodel)$coefficients[2,2]
    CI.mpp<-eff.mpp+z*se
    mpp<- CI.mpp<NImargin
    
    # iv with 2 stage regression
    se<-summary(asmm)$coefficients [2,2]
    CI.iv<-eff.iv + z*se
    iv<- CI.iv<NImargin
    
    return(c(iv, itt, mpp, ppp))
  }
}

plot.eff<- function(df,method,nIterations, true.effect, ymin, ymax){
  
  k<-which(method==analysis.method)
  x<-c()
  y<-c()
  
  for (i in 1:length(interval)) {
    x[[i]]<-df[[i]][,k]
    y[[i]]<-rep(interval[i],dim(df[[i]])[1])
  }
  
  plotdata<-matrix(c(unlist(x),unlist(y)), ncol = 2)
  
  plot<-ggplot(as.data.frame(plotdata), aes(x=V2, y=V1))+ 
    geom_point(aes(y=V1, colour=method, alpha=0.0001), size=0.2) + 
    geom_smooth(aes(y=V1, colour=method), method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants")+
    ylab("Effect estimate")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette[k])+
    theme(legend.title=element_blank(), legend.position="none", legend.text=element_text(size=legendfontsize))+
    scale_x_continuous(limits=c(0.58, 1), breaks = seq(0.6, 1, by=0.1))+
    scale_y_continuous(limits=c(ymin, ymax))+
    geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = 0.25),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
    geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = -0.05),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
    annotate(geom="text", x=0.58, y=0.11, angle = 90, label="Favour control", color="grey40", hjust=0) +
    annotate(geom="text", x=0.58, y=0.09, angle = 90, label="Favour experiment", color="grey40", hjust=1) +
    geom_hline(yintercept=true.effect, linetype='dashed', color='red', size=0.5)
  
  return(plot)
}

#################CASE 1 Non-compliance caused by non confounding process#########################
#################################################################################################

#BIAS
bias.nonconfounding<- function(n, p.experiment, p.stdcare, nIterations, interval, noncomply, true.effect,ymin, ymax){  
  
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
    estimate[[i]]<-matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  #plot
  iv<- plot.eff(df=estimate,analysis.method[1],nIterations=nIterations,  true.effect = true.effect ,  ymin=ymin, ymax=ymax)
  itt<- plot.eff(df=estimate,analysis.method[2],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  mpp<- plot.eff(df=estimate,analysis.method[3],nIterations=nIterations,  true.effect = true.effect,  ymin=ymin, ymax=ymax)
  pp<- plot.eff(df=estimate,analysis.method[4],nIterations=nIterations,  true.effect = true.effect,  ymin=ymin, ymax=ymax)
  
  bias.plot<-ggarrange(itt, pp,
                       mpp, iv,
                       ncol = 2, nrow = 2,
                       common.legend = FALSE )
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.nonconfounding<- function(n, p.experiment, nIterations, interval,NImargin, noncomply){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare=p.experiment-NImargin 
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
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
    estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE))
  }
  
  df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
  
  plot <- ggplot(df, aes(interval))+ 
    geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
    geom_point(aes(y=df$V2, colour=analysis.method[2], alpha=0.0001)) + 
    geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
    geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
    stat_smooth(aes(y=df$V1, colour=analysis.method[1]), se = FALSE)+
    stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
    stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
    stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants")+
    ylab("Type 1 error")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom", legend.text=element_text(size=legendfontsize))+
    theme(legend.title=element_blank())+
    scale_x_continuous(limits=c(start.interval, 1))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    geom_hline(yintercept=0.025, linetype='dashed', color='red', size=0.5)
  
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
      estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE))
    }
    
    df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
    
    plot <- ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.0001)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
      stat_smooth(aes(y=df$V1, colour=analysis.method[1]),  se = FALSE)+
      stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
      stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
      stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
 
      xlab("Proportion of adherent participants")+
      ylab("Power")+
      theme_minimal()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=legendfontsize))+
      scale_x_continuous(limits=c(start.interval, 1))+
      scale_y_continuous(breaks=seq(0,1,0.1))
    
    return(plot)
    
  }
  else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

############################CASE 2 Non-compliance caused by confounding process  ################################
################################################################################################################
#BIAS
bias.confounding<- function(n, p.experiment, p.stdcare, confounder.intervention, confounder.outcome,interval,nIterations,noncomply, true.effect, ymin, ymax){  
  
  #make up vectors for simulations 
  estimate<-c()  #for saving output from each interval
  .estimate<-c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata<-simdata.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysis(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]]<-matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  #plot
  iv<- plot.eff(df=estimate,analysis.method[1],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  itt<- plot.eff(df=estimate,analysis.method[2],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  mpp<- plot.eff(df=estimate,analysis.method[3],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  pp<- plot.eff(df=estimate,analysis.method[4],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  
  bias.plot<-ggarrange(itt, pp, 
                       mpp, iv,
                       ncol = 2, nrow = 2,
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
    estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE),na.rm = TRUE)
  }
  
  df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
  
  plot <- ggplot(df, aes(interval))+ 
    geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
    geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.0001)) + 
    geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
    geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
    stat_smooth(aes(y=df$V1, colour=analysis.method[1]),  se = FALSE)+
    stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
    stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
    stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants")+
    ylab("Type 1 error")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom", legend.text=element_text(size=legendfontsize), legend.title=element_blank())+
    scale_x_continuous(limits=c(start.interval, 1))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    geom_hline(yintercept=0.025, linetype='dashed', color='red', size=0.5)
  
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
      estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE),na.rm = TRUE)
    }
    
    df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
    
    plot <- ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.0001)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
      stat_smooth(aes(y=df$V1, colour=analysis.method[1]),  se = FALSE)+
      stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
      stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
      stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
 
      xlab("Proportion of adherent participants")+
      ylab("Power")+
      theme_minimal()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom", legend.text=element_text(size=legendfontsize), legend.title=element_blank())+
      scale_x_continuous(limits=c(start.interval, 1))+
      scale_y_continuous(breaks=seq(0,1,0.1))
    
    return(plot)
    
  }
  else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

############################CASE 3 Non-compliance caused by unknown confounding process  #######################
################################################################################################################
#BIAS
bias.unknownconfounding<- function(n, p.experiment, p.stdcare, confounder.intervention, confounder.outcome,interval,nIterations,noncomply, true.effect, ymin, ymax){  
  
  #make up vectors for simulations 
  estimate<-c()  #for saving output from each interval
  .estimate<-c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata<-simdata.unknownconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysisunknown(simdata=simdata, type='bias',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]]<-matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  #plot
  iv<- plot.eff(df=estimate,analysis.method[1],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  itt<- plot.eff(df=estimate,analysis.method[2],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  mpp<- plot.eff(df=estimate,analysis.method[3],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  pp<- plot.eff(df=estimate,analysis.method[4],nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  
  bias.plot<-ggarrange(itt, pp, 
                       mpp, iv,
                       ncol = 2, nrow = 2,
                       common.legend = FALSE )
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.unknownconfounding<- function(n, p.experiment, NImargin, confounder.intervention, confounder.outcome,interval,nIterations,noncomply){  
  
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
        simdata<-simdata.unknownconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i,confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysisunknown(simdata=simdata, type='type1',n=n,z=z,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    
    # mean of type 1 error from iterated data 
    estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE),na.rm = TRUE)
  }
  
  df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
  
  plot <- ggplot(df, aes(interval))+ 
    geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
    geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.0001)) + 
    geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
    geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
    stat_smooth(aes(y=df$V1, colour=analysis.method[1]),  se = FALSE)+
    stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
    stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
    stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants")+
    ylab("Type 1 error")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom", legend.text=element_text(size=legendfontsize), legend.title=element_blank())+
    scale_x_continuous(limits=c(start.interval, 1))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    geom_hline(yintercept=0.025, linetype='dashed', color='red', size=0.5)
  
  return(plot)
  
}

#POWER
power.unknownconfounding<- function(n, p.experiment, p.stdcare, NImargin, confounder.intervention, confounder.outcome,interval,nIterations,noncomply){  
  
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
          simdata<-simdata.unknownconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i,confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
          .estimate[[l]]<-analysisunknown(simdata=simdata, type='power',n=n,z=z,NImargin=NImargin)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
        
      }
      # mean of power from iterated data 
      estimate[[i]]<-colMeans(matrix(as.numeric(unlist(.estimate)), ncol = length(analysis.method), byrow = TRUE),na.rm = TRUE)
    }
    
    df<-as.data.frame(matrix(unlist(estimate),ncol=length(analysis.method), byrow=TRUE))
    
    plot <- ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2],alpha=0.0001)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
      stat_smooth(aes(y=df$V1, colour=analysis.method[1]),  se = FALSE)+
      stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
      stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
      stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
      xlab("Proportion of adherent participants")+
      ylab("Power")+
      theme_minimal()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom", legend.text=element_text(size=legendfontsize), legend.title=element_blank())+
      scale_x_continuous(limits=c(start.interval, 1))+
      scale_y_continuous(breaks=seq(0,1,0.1))
    
    return(plot)
    
  }
  else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#########analysis of multiple unknown covariates with IPW 

simdata.unknownconfounding.multi<- function(n, p.experiment, p.stdcare, noncomply,confounder.outcome, confounder.intervention,i){
  
  #compliance
  if (noncomply=="both") {
    comply.experiment<-interval
    comply.stdcare<-interval 
  } else if (noncomply=="experimental") {
    comply.experiment<- interval
    comply.stdcare<- rep(1, length(interval))
  } else { #noncomply=="stdcare"
    comply.experiment<-rep(1, length(interval))
    comply.stdcare<-interval
  }
  
  id = seq(1,(2*n), by=1) #create participant id  
  randomisation = c(rep(1,n), rep(0,n)) #randomisation
  known = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
  unknown1 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  unknown2 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  unknown3 = rep(rbeta(n=n,shape1=2,shape2=2),2)
  confounder=known+unknown1+unknown2+unknown3
  
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
  
  d<-matrix(data=c(id, randomisation, confounder, known, unknown1, unknown2, unknown3), nrow=(2*n))
  d.ordered<-matrix(cbind(d[order(d[,3]),], outcome0, outcome1),ncol=9)#order confounder in ascending order 
  d.grouped<-rbind(d.ordered[which(d.ordered[,2]==1),], d.ordered[which(d.ordered[,2]==0),])
  
  #INTERVENTION dependent on randomisation and confounders
  if (confounder.intervention=="Increase likelihood") {
    shape2<-runif(1, min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for intervention
    
    shape1<- shape2*(1-comply.stdcare[i])/(1-(1-comply.stdcare[i])) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)
    
  } else {
    shape2<-runif(1,min=2, max=10)
    shape1<- shape2*comply.experiment[i]/(1-comply.experiment[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2),decreasing = TRUE) #individual probability with mean of p.experiment, in increasing order
    int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
    
    shape1<- shape2*comply.stdcare[i]/(1-comply.stdcare[i]) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
    comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
    cont.intervention<- rbinom(n, 1, prob=1-comply.stdcare.ind)
  }
  
  intervention<-c(int.intervention,cont.intervention)
  
  #ACTUAL OUTCOMES depend on intervention
  outcome<-getoutcome(vector.outcome1=d.grouped[,9], vector.outcome0=d.grouped[,8], intervention)
  
  return(matrix(cbind(d.grouped[,c(-8,-9)],intervention,outcome), ncol=9))
}

analysisunknown.multi<- function(simdata,n, NImargin){
  
  ## intention to treat 
  pz1.value = mean(simdata[which(simdata[,2]==1),][,9])
  pz0.value = mean(simdata[which(simdata[,2]==0),][,9])
  eff.itt = pz1.value-pz0.value
  
  ## per protocol 
  pp = simdata[which(simdata[,2]==simdata[,8]),] # perprotocol population
  p.experiment.vector= pp[which(pp[,2]==1),][,9]
  p.experiment.value= mean(p.experiment.vector)
  p.stdcare.vector= pp[which(pp[,2]==0),][,9]
  p.stdcare.value= mean(p.stdcare.vector) 
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
  asmm =gmm(simdata[,9] ~ simdata[,8], x=simdata[,2], vcov="iid")
  eff.iv=summary(asmm)$coefficients [2,1]
  
  return(c(eff.iv, eff.itt, eff.mpp.all, eff.pp, eff.mpp.unknown2, eff.mpp.unknown1, eff.mpp.known))
  
}

bias.unknownconfounding.multi<- function(n, p.experiment, p.stdcare, confounder.intervention, confounder.outcome,interval,nIterations,noncomply, ymin, ymax){  
  
  #make up vectors for simulations 
  estimate<-c()  #for saving output from each interval
  .estimate<-c() #output from each simulation
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  #simulate and derive treatment effect 
  for(i in 1:length(interval)) { print(paste ("interval value",i,"out of", length(interval)))
    for(l in 1:nIterations)  {
      tryCatch({
        simdata<-simdata.unknownconfounding.multi(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, noncomply=noncomply,i=i, confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention)
        .estimate[[l]]<-analysisunknown.multi(simdata=simdata,n=n,NImargin=NImargin)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval
    estimate[[i]]<-matrix(unlist(.estimate), ncol = 7, byrow = TRUE)
  }
  
  #plot
  iv= plot.eff(df=estimate, method=analysis.method[1], nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  itt= plot.eff(df=estimate, analysis.method[2], nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  pp=plot.eff(df=estimate, analysis.method[4], nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  
  x<-c()
  
  for (h in 1:length(interval)) {
    x[[h]]=estimate[[h]][,c(3,5,6,7)]
  }
  
  mppdata=as.data.frame(matrix(c(unlist(x),
                                 rep(interval, each=length(x[[1]])), 
                                 rep(rep(1:4, each=nIterations),length(interval))), ncol = 3))
  colnames(mppdata)=c('y','x','type')
  mppdata$type=as.factor(mppdata$type)
  
  # The colour palette 
  cbPalette.multi <- c("#e69f00", "#E6AE00", "#E6BD00", "#E6CD00")
  
  mpp= ggplot(mppdata)+ 
    geom_point(aes(x=x, y=y, color=type, alpha=0.0001), size=0.2) +
    geom_smooth(aes(x=x,y=y, colour=type), method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants")+
    ylab("Treatment effect")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette.multi)+
    theme(legend.title=element_blank(), legend.text=element_text(size=legendfontsize), legend.position="none")+
    scale_x_continuous(limits=c(start.interval, 1))+
    scale_y_continuous(limits=c(ymin, ymax))+
    geom_hline(yintercept=true.effect, linetype='dashed', color='red', size=0.5)
  
  bias.plot<-ggarrange(itt, pp, mpp, iv,
                       ncol = 2, nrow = 2,
                       legend = "none")
  
  return(bias.plot)
} 


