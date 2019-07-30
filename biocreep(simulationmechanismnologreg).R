
################################################################################################################
###################Using causal inference to address non-Adherence in non inferiority trials###################
################################################################################################################

######## Set up #########
rm(list=ls()) # Clean working environment 

# Required libraries 
library(plotrix)
library(gtools) #for logit function
library(parallel) #use multiple cores
library(forcats) #ridge plot 
library(ggridges)
library(ggplot2)
library(gridExtra)

setwd("/Users/moyin/Desktop/angelsfly/NItrialsimulation/codes")

cl<-makeCluster(detectCores()-2)

upperCI<- function(n, comply.experiment, comply.stdcare){  
  
  #make up vectors for simulations 
  CI.itt<-.CI.itt<-c() 
  
  nIterations<-1000
  
  #alpha error and critical value 
  z = qnorm(0.975) #alpha=0.025 (one sided)
  
  #decreasing trial efficacy 
  p.experiment<- c(0.4,0.5,0.6,0.7)
  p.stdcare<- c(0.3,0.4,0.5,0.6)
  
  for (q in 1:4){
    #simulate and derive treatment effect 
    for(l in 1:nIterations) { 
      tryCatch({ #allow the function to run in case of errors
        #simulate trial data frame 
        
        #RANDOMISATION ratio 1:1
        simdata<- data.table(id=seq(1,(2*n), by=1), #create participant id 
                             randomisation=c(rep(1,n), rep(0,n))) #randomisation 
        
        #CONFOUNDER normal distribution ranging 0-1 
        simdata$confounder<- rbeta(n=(2*n),shape1=2,shape2=2)
        
        #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
        sd<-runif(1,min = 1, max=3) #sd of probability of outcome, indicates how strongly confounder affects outcome 
        logit.p.experiment<-logit(p.experiment[q])
        p.experiment.ind<- sort(inv.logit(rnorm(n=(2*n), mean=logit.p.experiment,sd=sd))) #arrange individual probabilities in ascending order
        simdata<-simdata[order(simdata$confounder),] #arrange confounder in ascending order 
        simdata$outcome1<- rbinom(2*n, 1, prob=p.experiment.ind) #increasing confounder value will have increasing probability for outcome 
        
        simdata$outcome0<- rbinom(2*n, 1, prob=p.experiment.ind-0.1)
        simdata$outcome0[is.na(simdata$outcome0)]<- 0
        
        #INTERVENTION dependent on randomisation and confounders
        simdata$intervention<- rep(NA,2*n)
        
        int<-simdata[which(simdata$randomisation==1),] #set of simulated data with randomisation =1
        cont<-simdata[which(simdata$randomisation==0),] #set of simulated data with randomisation =0
        
        nocomply<-which(int$confounder>quantile(int$confounder, comply.experiment)) #identify those who will not comply in intervention group (high confounder value will have control)
        int$intervention[nocomply]<-0
        int$intervention[which(is.na(int$intervention))]<-1
        
        comply<-which(cont$confounder>quantile(cont$confounder, 1-comply.stdcare)) #identify those who will comply in control group (high confounder value will have control)
        cont$intervention[comply]<-0
        cont$intervention[which(is.na(cont$intervention))]<-1
        
        simdata<- rbind(int,cont)
        
        #ACTUAL OUTCOMES depend on intervention
        simdata$outcome<- rep(NA,2*n)
        for (y in 1:(2*n)) {
          if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
            simdata$outcome[y]<-simdata$outcome0[y]}
        }
        
        #estimating treatment effect 
        ## intention to treat 
        pz1.value<- mean(simdata[which(simdata$randomisation==1),]$outcome)
        pz0.value<- mean(simdata[which(simdata$randomisation==0),]$outcome)
        eff.itt<- pz1.value-pz0.value    
        var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
        .CI.itt[l]<- eff.itt + z*sqrt(var.eff.itt)
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
      }
    CI.itt[[q]]<-.CI.itt
  }
  CImatrix<- matrix(unlist(CI.itt), nrow = 4, ncol=nIterations, byrow = TRUE)
  return(CImatrix)
}

trial1<-upperCI(n=500, comply.experiment = 1,comply.stdcare = 1)
trial2<-upperCI(n=500, comply.experiment = 0.9,comply.stdcare = 0.9)
trial3<-upperCI(n=500, comply.experiment = 0.8,comply.stdcare = 0.8)
trial4<-upperCI(n=500, comply.experiment = 0.7,comply.stdcare = 0.7)

stopCluster(cl)

fullset<-list(trial1, trial2, trial3, trial4)
save(fullset, file="fullset.Rdata")

#format data for plot 
load("fullset.RData")

formatforplot<- function(fullset,nIterations) {
  
  forplot<-c()
  
  for (i in 1:4) { #every i is a Adherence level 
    
    set1<-fullset[[i]][1,] #trial1
    set2<-fullset[[i]][2,] #trial2
    set3<-fullset[[i]][3,] #trial3
    set4<-fullset[[i]][4,] #trial4
    
    false1<- sum(set1<=0.1)/nIterations
    false2<- sum(set2<=0.1)/nIterations*false1
    false3<- sum(set3<=0.1)/nIterations*false2
    false4<- sum(set4<=0.1)/nIterations*false3
    
    set1.df<- as.data.frame(cbind(rep(1, nIterations), set1))
    colnames(set1.df)<-c('y','x')
    
    .set2.plot<- set2+0.1
    set2.diff<- 0.2-max(.set2.plot[{q<-rank(.set2.plot)/nIterations;q<false2}])
    if (set2.diff==Inf) {set2.diff=0.2-min(.set2.plot)}
    set2.plot<- .set2.plot+set2.diff
    set2.df<- as.data.frame(cbind(rep(2, nIterations), set2.plot))
    colnames(set2.df)<-c('y','x')
    
    .set3.plot<- set3+0.2
    set3.diff<- 0.3-max(.set3.plot[{q<-rank(.set3.plot)/nIterations;q<false3}])
    if (set3.diff==Inf) {set3.diff=0.3-min(.set3.plot)}
    set3.plot<- .set3.plot+set3.diff
    set3.df<- as.data.frame(cbind(rep(3, nIterations), set3.plot))
    colnames(set3.df)<-c('y','x')
    
    .set4.plot<- set4+0.3
    set4.diff<- 0.4-max(.set4.plot[{q<-rank(.set4.plot)/nIterations;q<false4}])
    if (set4.diff==Inf) {set4.diff=0.4-min(.set4.plot)}
    set4.plot<- .set4.plot+set4.diff
    set4.df<- as.data.frame(cbind(rep(4, nIterations), set4.plot))
    colnames(set4.df)<-c('y','x')
    
    set.df<- as.data.frame(rbind(set1.df, set2.df, set3.df,set4.df))
    set.df$y<-fct_rev(as.factor(set.df$y))
    
    forplot[[i]]<-list(set.df, c(false1,false2,false3,false4))
  }
  return(forplot)
}

plot<-function(set.df,Adherence){
  
  gg<- ggplot(set.df, aes(x = x, y = y, group = y),height = ..count..) +
    geom_density_ridges(aes(x = x, fill = paste(y)), 
                        alpha = .8, color = "grey", from = 0, to = 0.55) +
    labs(
      subtitle = paste(Adherence, "%","", "Adherence"),
      y= "", 
      x=""
    ) +
    scale_fill_cyclical(
      values = c("#0868ac", "#43a2ca","#7bccc4","#bae4bc")
    ) +
    scale_x_continuous(limits=c(-0.02,0.55), breaks = c(0,0.10, 0.20, 0.30,0.40), labels = c("0", "10", "20", "30", "40")) +
    scale_y_discrete(labels=c("C vs D","B vs C","A vs B", "Standard vs A"))+
    theme_ridges(grid = F) +theme(axis.text=element_text(size=10))
  
  d <- ggplot_build(gg)$data[[1]]
  .d1<- d[which(d$x <=0.1),]
  d1<- .d1[which(.d1$ymin ==4),]
  .d2<- d[which(d$x <=0.2),]
  d2<- .d2[which(.d2$ymin ==3),]
  .d3<- d[which(d$x <=0.3),]
  d3<- .d3[which(.d3$ymin ==2),]
  .d4<- d[which(d$x <=0.4),]
  d4<- .d4[which(.d4$ymin ==1),]
  
  # Add geom_ribbon for shaded area
  plot<- gg +
    geom_ribbon(
      data = transform(d1, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d2, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d3, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d4, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")
  
  return(plot)
}

plot1<-function(set.df,Adherence){
  
  gg<-ggplot(set.df, aes(x = x, y = y, group = y),height = ..count..) +
    geom_density_ridges(aes(x = x, fill = paste(y)), 
                        alpha = .8, color = "grey", from = 0, to = 0.55) +
    labs(
      subtitle = paste(Adherence, "%","", "Adherence"),
      y= "", 
      x="Decrease in treatment efficacy with respect to standard-of-care (%)"
    ) +
    scale_fill_cyclical(
      values = c("#0868ac", "#43a2ca","#7bccc4","#bae4bc")
    ) +
    scale_x_continuous(limits=c(-0.02,0.55), breaks = c(0,0.10, 0.20, 0.30,0.40),labels = c("0", "10", "20", "30", "40")) +
    scale_y_discrete(labels=c("C vs D","B vs C","A vs B", "Standard vs A"))+ 
    theme_ridges(grid=F) +theme(axis.text=element_text(size=10), axis.title.x = element_text(hjust =0.5, size=10))
  
  d <- ggplot_build(gg)$data[[1]]
  .d1<- d[which(d$x <=0.1),]
  d1<- .d1[which(.d1$ymin ==4),]
  .d2<- d[which(d$x <=0.2),]
  d2<- .d2[which(.d2$ymin ==3),]
  .d3<- d[which(d$x <=0.3),]
  d3<- .d3[which(.d3$ymin ==2),]
  .d4<- d[which(d$x <=0.4),]
  d4<- .d4[which(.d4$ymin ==1),]
  
  # Add geom_ribbon for shaded area
  plot<- gg +
    geom_ribbon(
      data = transform(d1, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d2, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d3, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d4, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")
  
  return(plot)
  
}


###########################################
###########################################
forplot<-formatforplot(fullset,1000)
plot100<- plot(forplot[[1]][[1]],100)
plot90<- plot(forplot[[2]][[1]],90)
plot80<- plot(forplot[[3]][[1]],80)
plot70<- plot1(forplot[[4]][[1]],70)

grid.arrange(plot100, plot90,plot80,plot70,nrow = 4)




