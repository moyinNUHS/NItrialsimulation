################################################################################################################
###################Using causal inference to address non-compliance in non inferiority trials###################
################################################################################################################

######## Set up #########
setwd("/Users/moyin/Desktop/VAP studd/Causal inference simulation/graphs") #set working directory 
rm(list=ls()) # Clean working environment 

# Required libraries 
library(Hmisc); library(rms); library(gsDesign);library(ivpack); library(gmm)
library(data.table); library(dplyr); library(plotly); library(ggpubr); library(ggplot2); library(gridExtra)
library(Matching); library(tableone); library(MatchIt); library(geepack); library(scales)

# The colour blind friendly palette with black:
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
col.itt<- "#009E73" 
col.pp <- "#D55E00"
col.at<-  "#E69F00" 
col.iv<-  "#56B4E9"
col.mpp<- "#0072B2" 
col.ps<-  "#CC79A7"

#set seed for simulation 
set.seed(1234)

#################CASE 1.1 Non-compliance caused by random process (affect both groups) #########################
################################################################################################################
#BIAS
simdata.bias.1.1<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() #output from each simulation
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-50 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER unknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p1+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100%,number of people who complied to protocol
                compliance<- (sample(1:(2*n), compliance[i]))       # IDs of participants who complied to protocol
                for (d in 1: (2*n)) {        
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                    if ((is.na(simdata$intervention))[d]) { 
                        if (simdata$randomisation[d] == 1) {
                            simdata$intervention[d] = 0
                        } else {
                            simdata$intervention [d]= 1
                        }
                    }
                }
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                y<- simdata$outcome
                x<- simdata$intervention
                z<-simdata$randomisation
                data<-data.frame(y,x,z)
                asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                .eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                
                #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                #.eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
            #save results of every iteration 
            eff.itt[[i]]<- .eff.itt
            eff.pp[[i]]<- .eff.pp
            eff.at[[i]]<- .eff.at
            eff.iv[[i]]<- .eff.iv
            eff.mpp[[i]]<- .eff.mpp
            eff.ps[[i]]<- .eff.ps
            
        }
        
        # make a dataframe of the outcomes of the iterations 
        eff.itt.df<- cbind(prop.c,t(data.frame(eff.itt))) 
        eff.pp.df<- cbind(prop.c,t(data.frame(eff.pp))) 
        eff.at.df<- cbind(prop.c,t(data.frame(eff.at))) 
        eff.iv.df<- cbind(prop.c,t(data.frame(eff.iv))) 
        eff.mpp.df<- cbind(prop.c,t(data.frame(eff.mpp))) 
        eff.ps.df<- cbind(prop.c,t(data.frame(eff.ps))) 
        
        return(list(eff.itt.df, eff.pp.df,eff.at.df,eff.iv.df,eff.mpp.df,eff.ps.df))
    }
    
    sim.df<-sim()
    
    itt.df<-as.data.frame(sim.df[1])
    colnames(itt.df)<- c('prop.c', 1:nIterations)
    row.names(itt.df)<- c(1:interval)
    itt.melt<- melt(itt.df, id.vars='prop.c')
    itt<-ggplot(itt.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=itt.melt$value, colour="Intention to treat", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=itt.melt$value, colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.itt)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    pp.df<-as.data.frame(sim.df[2])
    colnames(pp.df)<- c('prop.c', 1:nIterations)
    row.names(pp.df)<- c(1:interval)
    pp.melt<- melt(pp.df, id.vars='prop.c')
    pp<-ggplot(pp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=pp.melt$value, colour="Per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=pp.melt$value, colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.pp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    at.df<-as.data.frame(sim.df[3])
    colnames(at.df)<- c('prop.c', 1:nIterations)
    row.names(at.df)<- c(1:interval)
    at.melt<- melt(at.df, id.vars='prop.c')
    at<-ggplot(at.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=at.melt$value, colour="As treated", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=at.melt$value, colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.at)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    iv.df<-as.data.frame(sim.df[4])
    colnames(iv.df)<- c('prop.c', 1:nIterations)
    row.names(iv.df)<- c(1:interval)
    iv.melt<- melt(iv.df, id.vars='prop.c')
    iv<-ggplot(iv.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=iv.melt$value, colour="Instrumental variable", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=iv.melt$value, colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.iv)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    mpp.df<-as.data.frame(sim.df[5])
    colnames(mpp.df)<- c('prop.c', 1:nIterations)
    row.names(mpp.df)<- c(1:interval)
    mpp.melt<- melt(mpp.df, id.vars='prop.c')
    mpp<-ggplot(mpp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=mpp.melt$value, colour="Modified per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=mpp.melt$value, colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.mpp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    ps.df<-as.data.frame(sim.df[6])
    colnames(ps.df)<- c('prop.c', 1:nIterations)
    row.names(ps.df)<- c(1:interval)
    ps.melt<- melt(ps.df, id.vars='prop.c')
    ps<-ggplot(ps.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=ps.melt$value, colour="Propensity score matching", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=ps.melt$value, colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.ps)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    bias.plot<-ggarrange(itt, pp, at,
                         iv, mpp, ps,
                         ncol = 3, nrow = 2,
                         common.legend = FALSE )
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.1.1<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #alpha error and critical value 
    alpha=0.025
    z = qnorm(1-0.05/2)
    
    #number of data points in simulation 
    interval<-50 
    
    #make up vectors for simulations 
    .type1.error.at<-.type1.error.itt<-.type1.error.iv<-.type1.error.mpp<-.type1.error.pp<-.type1.error.ps<-c() 
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) {  
                    
                    tryCatch({ #allow the function to run in case of errors
                        
                        #simulate data frame 
                        
                        #RANDOMISATION ratio 1:1
                        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                        simdata$id<- seq(1,(2*n), by=1) #create participant id 
                        
                        #CONFOUNDER normal distribution ranging 0-1 
                        simdata$confounder<- rep(NA,2*n) 
                        simdata$confounder<- rnorm(2*n) 
                        simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #CONFOUNDER nknown normal distribution ranging 0-1 
                        simdata$confounder.u<- rep(NA,2*n) 
                        simdata$confounder.u<- rnorm(2*n) 
                        simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                        simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                        simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                        
                        prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                        simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                        
                        prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                        simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                        
                        #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                        simdata$intervention <- rep(NA,2*n)
                        
                        compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100% 
                        compliance<- (sample(1:(2*n), compliance[i]))
                        for (d in 1: (2*n)) {        
                            if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                                simdata$intervention[d] = simdata$randomisation[d]
                            } 
                        }
                        
                        for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                            if ((is.na(simdata$intervention))[d]) { 
                                if (simdata$randomisation[d] == 1) {
                                    simdata$intervention[d] = 0
                                } else {
                                    simdata$intervention [d]= 1
                                }
                            }
                        }
                        
                        #ACTUAL OUTCOMES depend on intervention
                        simdata$outcome<- rep(NA,2*n)
                        for (y in 1:(2*n)) {
                            if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                                simdata$outcome[y]<-simdata$outcome0[y]}
                        }
                        
                        #estimating treatment effect 
                        ## intention to treat 
                        pz1.vector<- (filter(simdata, randomisation==1))$outcome
                        pz1.value<- mean(pz1.vector)  
                        pz0.vector<- (filter(simdata, randomisation==0))$outcome
                        pz0.value<- mean(pz0.vector)
                        eff.itt<- pz1.value-pz0.value    
                        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
                        CI.itt<- eff.itt + c(-1,1)*z*sqrt(var.eff.itt)
                        .type1.error.itt[l]<-CI.itt[2]<NImargin
                        
                        ## per protocol 
                        pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                        .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome
                        p11.value<- mean(p11.vector)
                        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome
                        p00.value<- mean(p00.vector) 
                        eff.pp <- p11.value-p00.value   
                        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
                        CI.pp<- eff.pp + c(-1,1)*z*sqrt(var.eff.pp)
                        .type1.error.pp[l]<- CI.pp[2]<NImargin
                        
                        ## as treated
                        pd1.vector<- (filter(simdata, intervention==1))$outcome
                        pd1.value<- mean(pd1.vector)
                        pd0.vector<- (filter(simdata, intervention==0))$outcome
                        pd0.value<- mean(pd0.vector)
                        eff.at <- pd1.value-pd0.value        
                        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
                        CI.at<- eff.at + c(-1,1)*z*sqrt(var.eff.at)
                        .type1.error.at[l]<- CI.at[2]<NImargin
                        
                        ## linear binomial generalized linear models with inverse probability weights
                        E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                        ps<-predict(E.out, type="response")
                        sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                        mpp<- geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)
                        eff.mpp<-mpp$coefficients['intervention']# fit linear binomial model to the weighted data
                        se<-summary(mpp)$coefficients[2,2] 
                        CI.mpp<-eff.mpp + c(-1,1)*z*se
                        .type1.error.mpp[l]<- CI.mpp[2]<NImargin
                        
                        ## propensity score matching 
                        m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                        pscore<- m$fitted.values
                        mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                        match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                        matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                        y_trt<-matched$outcome[matched$intervention==1]
                        y_con<-matched$outcome[matched$intervention==0]
                        diffy<-y_trt-y_con
                        eff.ps<-(t.test(diffy))$estimate
                        CI.psm<-(t.test(diffy))$conf.int
                        .type1.error.ps[l]<- CI.psm[2]<NImargin
                        
                        # iv with 2 stage regression
                        y<- simdata$outcome
                        x<- simdata$intervention
                        z<-simdata$randomisation
                        data<-data.frame(y,x,z)
                        asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                        eff.iv<-(summary(asmm))$ coefficients [2,1]
                        se<-(summary(asmm))$ coefficients [2,2]
                        CI.iv<-eff.iv + c(-1,1)*z*se
                        .type1.error.iv[l]<- CI.iv[2]<NImargin
                        #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                        #eff.iv<-ivmodel$coef[2] 
                        #se<-robust.se(ivmodel)[2,2]
                        #CI.iv<-eff.iv + c(-1,1)*z*se
                       # .type1.error.iv[l]<- CI.iv[2]<NImargin
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(.type1.error.itt, na.rm=TRUE)
                type1.error.pp[i] <- mean(.type1.error.pp, na.rm=TRUE)
                type1.error.at[i] <- mean(.type1.error.at, na.rm=TRUE)
                type1.error.mpp[i]<- mean(.type1.error.mpp, na.rm=TRUE)
                type1.error.iv[i] <- mean(.type1.error.iv, na.rm=TRUE)
                type1.error.ps[i] <- mean(.type1.error.ps, na.rm=TRUE)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.1.1<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-50
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100% 
                    compliance<- (sample(1:(2*n), compliance[i]))
                    for (d in 1: (2*n)) {        
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                        if ((is.na(simdata$intervention))[d]) { 
                            if (simdata$randomisation[d] == 1) {
                                simdata$intervention[d] = 0
                            } else {
                                simdata$intervention [d]= 1
                            }
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    y<- simdata$outcome
                    x<- simdata$intervention
                    z<-simdata$randomisation
                    data<-data.frame(y,x,z)
                    asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                    eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                    
                    #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    #.eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.1
b1.1<-simdata.bias.1.1   (n=300, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=200)
t1.1<-simdata.type1.1.1  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.1<- simdata.power.1.1 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case1.1<-ggarrange(b1.1, 
                   ggarrange (t1.1, p1.1, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case1.1, top = text_grob("Case 1.1: Non-compliance caused by random processes affecting both groups", face = "bold", size = 12))
ggsave(filename = paste("case1.1", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

#################CASE 1.2 Non-compliance caused by random process (affect intervention group) ##################
################################################################################################################
#BIAS
simdata.bias.1.2<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() #output from each simulation
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-50 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                for (d in 1: (2*n)) {        
                    if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                compliance<- (sample(1:n, compliance[i]))
                
                for (d in 1:n) { 
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } else {
                        simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                    }
                }
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                y<- simdata$outcome
                x<- simdata$intervention
                z<-simdata$randomisation
                data<-data.frame(y,x,z)
                asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                .eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                
                #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                #.eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
            #save results of every iteration 
            eff.itt[[i]]<- .eff.itt
            eff.pp[[i]]<- .eff.pp
            eff.at[[i]]<- .eff.at
            eff.iv[[i]]<- .eff.iv
            eff.mpp[[i]]<- .eff.mpp
            eff.ps[[i]]<- .eff.ps
            
        }
        
        # make a dataframe of the outcomes of the iterations 
        eff.itt.df<- cbind(prop.c,t(data.frame(eff.itt))) 
        eff.pp.df<- cbind(prop.c,t(data.frame(eff.pp))) 
        eff.at.df<- cbind(prop.c,t(data.frame(eff.at))) 
        eff.iv.df<- cbind(prop.c,t(data.frame(eff.iv))) 
        eff.mpp.df<- cbind(prop.c,t(data.frame(eff.mpp))) 
        eff.ps.df<- cbind(prop.c,t(data.frame(eff.ps))) 
        
        return(list(eff.itt.df, eff.pp.df,eff.at.df,eff.iv.df,eff.mpp.df,eff.ps.df))
    }
    
    sim.df<-sim()
    
    itt.df<-as.data.frame(sim.df[1])
    colnames(itt.df)<- c('prop.c', 1:nIterations)
    row.names(itt.df)<- c(1:interval)
    itt.melt<- melt(itt.df, id.vars='prop.c')
    itt<-ggplot(itt.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=itt.melt$value, colour="Intention to treat", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=itt.melt$value, colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.itt)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    pp.df<-as.data.frame(sim.df[2])
    colnames(pp.df)<- c('prop.c', 1:nIterations)
    row.names(pp.df)<- c(1:interval)
    pp.melt<- melt(pp.df, id.vars='prop.c')
    pp<-ggplot(pp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=pp.melt$value, colour="Per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=pp.melt$value, colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.pp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    at.df<-as.data.frame(sim.df[3])
    colnames(at.df)<- c('prop.c', 1:nIterations)
    row.names(at.df)<- c(1:interval)
    at.melt<- melt(at.df, id.vars='prop.c')
    at<-ggplot(at.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=at.melt$value, colour="As treated", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=at.melt$value, colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.at)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    iv.df<-as.data.frame(sim.df[4])
    colnames(iv.df)<- c('prop.c', 1:nIterations)
    row.names(iv.df)<- c(1:interval)
    iv.melt<- melt(iv.df, id.vars='prop.c')
    iv<-ggplot(iv.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=iv.melt$value, colour="Instrumental variable", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=iv.melt$value, colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.iv)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    mpp.df<-as.data.frame(sim.df[5])
    colnames(mpp.df)<- c('prop.c', 1:nIterations)
    row.names(mpp.df)<- c(1:interval)
    mpp.melt<- melt(mpp.df, id.vars='prop.c')
    mpp<-ggplot(mpp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=mpp.melt$value, colour="Modified per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=mpp.melt$value, colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.mpp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    ps.df<-as.data.frame(sim.df[6])
    colnames(ps.df)<- c('prop.c', 1:nIterations)
    row.names(ps.df)<- c(1:interval)
    ps.melt<- melt(ps.df, id.vars='prop.c')
    ps<-ggplot(ps.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=ps.melt$value, colour="Propensity score matching", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=ps.melt$value, colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.ps)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    bias.plot<-ggarrange(itt, pp, at,
                         iv, mpp, ps,
                         ncol = 3, nrow = 2, common.legend = FALSE )
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.1.2<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #alpha error and critical value 
    alpha=0.025
    z = qnorm(1-0.05/2)
    
    #number of data points in simulation 
    interval<-50 
    
    #make up vectors for simulations 
    .type1.error.at<-.type1.error.itt<-.type1.error.iv<-.type1.error.mpp<-.type1.error.pp<-.type1.error.ps<-c() 
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) {  
                    
                    tryCatch({ #allow the function to run in case of errors
                        
                        #simulate data frame 
                        
                        #RANDOMISATION ratio 1:1
                        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                        simdata$id<- seq(1,(2*n), by=1) #create participant id 
                        
                        #CONFOUNDER normal distribution ranging 0-1 
                        simdata$confounder<- rep(NA,2*n) 
                        simdata$confounder<- rnorm(2*n) 
                        simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #CONFOUNDER nknown normal distribution ranging 0-1 
                        simdata$confounder.u<- rep(NA,2*n) 
                        simdata$confounder.u<- rnorm(2*n) 
                        simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                        simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                        simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                        
                        prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                        simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                        
                        prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                        simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                        
                        #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                        simdata$intervention <- rep(NA,2*n)
                        
                        for (d in 1: (2*n)) {        
                            if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                                simdata$intervention[d] = simdata$randomisation[d]
                            } 
                        }
                        
                        compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                        compliance<- (sample(1:n, compliance[i]))
                        
                        for (d in 1:n) { 
                            if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                                simdata$intervention[d] = simdata$randomisation[d]
                            } else {
                                simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                            }
                        }
                        
                        #ACTUAL OUTCOMES depend on intervention
                        simdata$outcome<- rep(NA,2*n)
                        for (y in 1:(2*n)) {
                            if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                                simdata$outcome[y]<-simdata$outcome0[y]}
                        }
                        
                        #estimating type 1 error
                        ## intention to treat 
                        pz1.vector<- (filter(simdata, randomisation==1))$outcome
                        pz1.value<- mean(pz1.vector)  
                        pz0.vector<- (filter(simdata, randomisation==0))$outcome
                        pz0.value<- mean(pz0.vector)
                        eff.itt<- pz1.value-pz0.value    
                        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
                        CI.itt<- eff.itt + c(-1,1)*z*sqrt(var.eff.itt)
                        .type1.error.itt[l]<-CI.itt[2]<NImargin
                        
                        ## per protocol 
                        pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                        .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome
                        p11.value<- mean(p11.vector)
                        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome
                        p00.value<- mean(p00.vector) 
                        eff.pp <- p11.value-p00.value   
                        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
                        CI.pp<- eff.pp + c(-1,1)*z*sqrt(var.eff.pp)
                        .type1.error.pp[l]<- CI.pp[2]<NImargin
                        
                        ## as treated
                        pd1.vector<- (filter(simdata, intervention==1))$outcome
                        pd1.value<- mean(pd1.vector)
                        pd0.vector<- (filter(simdata, intervention==0))$outcome
                        pd0.value<- mean(pd0.vector)
                        eff.at <- pd1.value-pd0.value        
                        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
                        CI.at<- eff.at + c(-1,1)*z*sqrt(var.eff.at)
                        .type1.error.at[l]<- CI.at[2]<NImargin
                        
                        ## linear binomial generalized linear models with inverse probability weights
                        E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                        ps<-predict(E.out, type="response")
                        sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                        mpp<- geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)
                        eff.mpp<-mpp$coefficients['intervention']# fit linear binomial model to the weighted data
                        se<-summary(mpp)$coefficients[2,2] 
                        CI.mpp<-eff.mpp + c(-1,1)*z*se
                        .type1.error.mpp[l]<- CI.mpp[2]<NImargin
                        
                        ## propensity score matching 
                        m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                        pscore<- m$fitted.values
                        mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                        match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                        matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                        y_trt<-matched$outcome[matched$intervention==1]
                        y_con<-matched$outcome[matched$intervention==0]
                        diffy<-y_trt-y_con
                        eff.ps<-(t.test(diffy))$estimate
                        CI.psm<-(t.test(diffy))$conf.int
                        .type1.error.ps[l]<- CI.psm[2]<NImargin
                        
                        # iv with 2 stage regression
                        y<- simdata$outcome
                        x<- simdata$intervention
                        z<-simdata$randomisation
                        data<-data.frame(y,x,z)
                        asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                        eff.iv<-(summary(asmm))$ coefficients [2,1]
                        se<-(summary(asmm))$ coefficients [2,2]
                        CI.iv<-eff.iv + c(-1,1)*z*se
                        .type1.error.iv[l]<- CI.iv[2]<NImargin
                        #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                        #eff.iv<-ivmodel$coef[2] 
                        #se<-robust.se(ivmodel)[2,2]
                        #CI.iv<-eff.iv + c(-1,1)*z*se
                        # .type1.error.iv[l]<- CI.iv[2]<NImargin
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(.type1.error.itt, na.rm=TRUE)
                type1.error.pp[i] <- mean(.type1.error.pp, na.rm=TRUE)
                type1.error.at[i] <- mean(.type1.error.at, na.rm=TRUE)
                type1.error.mpp[i]<- mean(.type1.error.mpp, na.rm=TRUE)
                type1.error.iv[i] <- mean(.type1.error.iv, na.rm=TRUE)
                type1.error.ps[i] <- mean(.type1.error.ps, na.rm=TRUE)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.1.2<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-50
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample(1:n, compliance[i]))
                    
                    for (d in 1:n) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    y<- simdata$outcome
                    x<- simdata$intervention
                    z<-simdata$randomisation
                    data<-data.frame(y,x,z)
                    asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                    eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                    
                    #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    #.eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.2
b1.2<-simdata.bias.1.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=200)
t1.2<-simdata.type1.1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.2<- simdata.power.1.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case1.2<-ggarrange(b1.2, 
                   ggarrange (t1.2, p1.2, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case1.2, top = text_grob("Case 1.2: Non-compliance caused by random processes affecting intervention group", face = "bold", size = 12))
ggsave(filename = paste("case1.2", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

#################CASE 1.3 Non-compliance caused by random process (affect control group) #######################
##################################################Not shown#####################################################
#BIAS
simdata.bias.1.3<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() #output from each simulation
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-50 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                for (d in 1: (2*n)) {        
                    if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                compliance<- (sample((n+1):(2*n), compliance[i]))
                
                for (d in (n+1):(2*n)) { 
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } else {
                        simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                    }
                }
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                y<- simdata$outcome
                x<- simdata$intervention
                z<-simdata$randomisation
                data<-data.frame(y,x,z)
                asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                .eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                
                #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                #.eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
            #save results of every iteration 
            eff.itt[[i]]<- .eff.itt
            eff.pp[[i]]<- .eff.pp
            eff.at[[i]]<- .eff.at
            eff.iv[[i]]<- .eff.iv
            eff.mpp[[i]]<- .eff.mpp
            eff.ps[[i]]<- .eff.ps
            
        }
        
        # make a dataframe of the outcomes of the iterations 
        eff.itt.df<- cbind(prop.c,t(data.frame(eff.itt))) 
        eff.pp.df<- cbind(prop.c,t(data.frame(eff.pp))) 
        eff.at.df<- cbind(prop.c,t(data.frame(eff.at))) 
        eff.iv.df<- cbind(prop.c,t(data.frame(eff.iv))) 
        eff.mpp.df<- cbind(prop.c,t(data.frame(eff.mpp))) 
        eff.ps.df<- cbind(prop.c,t(data.frame(eff.ps))) 
        
        return(list(eff.itt.df, eff.pp.df,eff.at.df,eff.iv.df,eff.mpp.df,eff.ps.df))
    }
    
    sim.df<-sim()
    
    itt.df<-as.data.frame(sim.df[1])
    colnames(itt.df)<- c('prop.c', 1:nIterations)
    row.names(itt.df)<- c(1:interval)
    itt.melt<- melt(itt.df, id.vars='prop.c')
    itt<-ggplot(itt.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=itt.melt$value, colour="Intention to treat", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=itt.melt$value, colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.itt)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    pp.df<-as.data.frame(sim.df[2])
    colnames(pp.df)<- c('prop.c', 1:nIterations)
    row.names(pp.df)<- c(1:interval)
    pp.melt<- melt(pp.df, id.vars='prop.c')
    pp<-ggplot(pp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=pp.melt$value, colour="Per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=pp.melt$value, colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.pp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    at.df<-as.data.frame(sim.df[3])
    colnames(at.df)<- c('prop.c', 1:nIterations)
    row.names(at.df)<- c(1:interval)
    at.melt<- melt(at.df, id.vars='prop.c')
    at<-ggplot(at.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=at.melt$value, colour="As treated", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=at.melt$value, colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.at)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    iv.df<-as.data.frame(sim.df[4])
    colnames(iv.df)<- c('prop.c', 1:nIterations)
    row.names(iv.df)<- c(1:interval)
    iv.melt<- melt(iv.df, id.vars='prop.c')
    iv<-ggplot(iv.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=iv.melt$value, colour="Instrumental variable", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=iv.melt$value, colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.iv)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    mpp.df<-as.data.frame(sim.df[5])
    colnames(mpp.df)<- c('prop.c', 1:nIterations)
    row.names(mpp.df)<- c(1:interval)
    mpp.melt<- melt(mpp.df, id.vars='prop.c')
    mpp<-ggplot(mpp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=mpp.melt$value, colour="Modified per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=mpp.melt$value, colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.mpp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    ps.df<-as.data.frame(sim.df[6])
    colnames(ps.df)<- c('prop.c', 1:nIterations)
    row.names(ps.df)<- c(1:interval)
    ps.melt<- melt(ps.df, id.vars='prop.c')
    ps<-ggplot(ps.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=ps.melt$value, colour="Propensity score matching", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=ps.melt$value, colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.ps)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    bias.plot<-ggarrange(itt, pp, at,
                         iv, mpp, ps,
                         ncol = 3, nrow = 2,
                         common.legend = FALSE )
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.1.3<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #alpha error and critical value 
    alpha=0.025
    z = qnorm(1-0.05/2)
    
    #number of data points in simulation 
    interval<-50 
    
    #make up vectors for simulations 
    .type1.error.at<-.type1.error.itt<-.type1.error.iv<-.type1.error.mpp<-.type1.error.pp<-.type1.error.ps<-c()                  
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    
                    tryCatch({ #allow the function to run in case of errors
                        
                        #simulate data frame 
                        
                        #RANDOMISATION ratio 1:1
                        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                        simdata$id<- seq(1,(2*n), by=1) #create participant id 
                        
                        #CONFOUNDER normal distribution ranging 0-1 
                        simdata$confounder<- rep(NA,2*n) 
                        simdata$confounder<- rnorm(2*n) 
                        simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #CONFOUNDER nknown normal distribution ranging 0-1 
                        simdata$confounder.u<- rep(NA,2*n) 
                        simdata$confounder.u<- rnorm(2*n) 
                        simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                        simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                        simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                        
                        prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                        simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                        
                        prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                        simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                        
                        #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                        simdata$intervention <- rep(NA,2*n)
                        
                        for (d in 1: (2*n)) {        
                            if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                                simdata$intervention[d] = simdata$randomisation[d]
                            } 
                        }
                        
                        compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                        compliance<- (sample((n+1):(2*n), compliance[i]))
                        
                        for (d in (n+1):(2*n)) { 
                            if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                                simdata$intervention[d] = simdata$randomisation[d]
                            } else {
                                simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                            }
                        }
                        
                        #ACTUAL OUTCOMES depend on intervention
                        simdata$outcome<- rep(NA,2*n)
                        for (y in 1:(2*n)) {
                            if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                                simdata$outcome[y]<-simdata$outcome0[y]}
                        }
                        
                        #estimating type 1 error
                        ## intention to treat 
                        pz1.vector<- (filter(simdata, randomisation==1))$outcome
                        pz1.value<- mean(pz1.vector)  
                        pz0.vector<- (filter(simdata, randomisation==0))$outcome
                        pz0.value<- mean(pz0.vector)
                        eff.itt<- pz1.value-pz0.value    
                        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
                        CI.itt<- eff.itt + c(-1,1)*z*sqrt(var.eff.itt)
                        .type1.error.itt[l]<-CI.itt[2]<NImargin
                        
                        ## per protocol 
                        pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                        .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome
                        p11.value<- mean(p11.vector)
                        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome
                        p00.value<- mean(p00.vector) 
                        eff.pp <- p11.value-p00.value   
                        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
                        CI.pp<- eff.pp + c(-1,1)*z*sqrt(var.eff.pp)
                        .type1.error.pp[l]<- CI.pp[2]<NImargin
                        
                        ## as treated
                        pd1.vector<- (filter(simdata, intervention==1))$outcome
                        pd1.value<- mean(pd1.vector)
                        pd0.vector<- (filter(simdata, intervention==0))$outcome
                        pd0.value<- mean(pd0.vector)
                        eff.at <- pd1.value-pd0.value        
                        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
                        CI.at<- eff.at + c(-1,1)*z*sqrt(var.eff.at)
                        .type1.error.at[l]<- CI.at[2]<NImargin
                        
                        ## linear binomial generalized linear models with inverse probability weights
                        E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                        ps<-predict(E.out, type="response")
                        sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                        mpp<- geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)
                        eff.mpp<-mpp$coefficients['intervention']# fit linear binomial model to the weighted data
                        se<-summary(mpp)$coefficients[2,2] 
                        CI.mpp<-eff.mpp + c(-1,1)*z*se
                        .type1.error.mpp[l]<- CI.mpp[2]<NImargin
                        
                        ## propensity score matching 
                        m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                        pscore<- m$fitted.values
                        mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                        match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                        matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                        y_trt<-matched$outcome[matched$intervention==1]
                        y_con<-matched$outcome[matched$intervention==0]
                        diffy<-y_trt-y_con
                        eff.ps<-(t.test(diffy))$estimate
                        CI.psm<-(t.test(diffy))$conf.int
                        .type1.error.ps[l]<- CI.psm[2]<NImargin
                        
                        # iv with 2 stage regression
                        y<- simdata$outcome
                        x<- simdata$intervention
                        z<-simdata$randomisation
                        data<-data.frame(y,x,z)
                        asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                        eff.iv<-(summary(asmm))$ coefficients [2,1]
                        se<-(summary(asmm))$ coefficients [2,2]
                        CI.iv<-eff.iv + c(-1,1)*z*se
                        .type1.error.iv[l]<- CI.iv[2]<NImargin
                        #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                        #eff.iv<-ivmodel$coef[2] 
                        #se<-robust.se(ivmodel)[2,2]
                        #CI.iv<-eff.iv + c(-1,1)*z*se
                        # .type1.error.iv[l]<- CI.iv[2]<NImargin
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(.type1.error.itt, na.rm=TRUE)
                type1.error.pp[i] <- mean(.type1.error.pp, na.rm=TRUE)
                type1.error.at[i] <- mean(.type1.error.at, na.rm=TRUE)
                type1.error.mpp[i]<- mean(.type1.error.mpp, na.rm=TRUE)
                type1.error.iv[i] <- mean(.type1.error.iv, na.rm=TRUE)
                type1.error.ps[i] <- mean(.type1.error.ps, na.rm=TRUE)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.1.3<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-50
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample((n+1):(2*n), compliance[i]))
                    
                    for (d in (n+1):(2*n)) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    y<- simdata$outcome
                    x<- simdata$intervention
                    z<-simdata$randomisation
                    data<-data.frame(y,x,z)
                    asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                    eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                    
                    #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    #.eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.3
b1.3<-simdata.bias.1.3   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=200)
t1.3<-simdata.type1.1.3  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.3<- simdata.power.1.3 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case1.3<-ggarrange(b1.3, 
                   ggarrange (t1.3, p1.3, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case1.3, top = text_grob("Case 1.3: Non-compliance caused by random processes affecting control group", face = "bold", size = 12))
ggsave(filename = paste("case1.3", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

############################CASE 2 Non-compliance caused by non random process  ################################
################################################################################################################
#BIAS
simdata.bias.2<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() 
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-50 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                .b0<-logit (i0)
                b0<- seq (.b0,logit(0.01), length.out = interval)
                .b1<-log (confounder.eff.i)
                b1<- seq (.b1,0, length.out = interval)
                .b2<-log (confounder.u.eff.i)
                b2<- seq (.b2,0, length.out = interval)
                .b3<-logit (i1)- logit (i0)
                b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                
                if (b1[i]==0) { pi<- simdata$randomisation} else {
                    logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                    pi<-exp(logit.pi)/(1+exp(logit.pi))
                }
                simdata$intervention <- rbinom(2*n,1,prob=pi)
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                y<- simdata$outcome
                x<- simdata$intervention
                z<-simdata$randomisation
                data<-data.frame(y,x,z)
                asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                .eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                
                #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                #.eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
            #save results of every iteration 
            eff.itt[[i]]<- .eff.itt
            eff.pp[[i]]<- .eff.pp
            eff.at[[i]]<- .eff.at
            eff.iv[[i]]<- .eff.iv
            eff.mpp[[i]]<- .eff.mpp
            eff.ps[[i]]<- .eff.ps
            
        }
        
        # make a dataframe of the outcomes of the iterations 
        eff.itt.df<- cbind(prop.c,t(data.frame(eff.itt))) 
        eff.pp.df<- cbind(prop.c,t(data.frame(eff.pp))) 
        eff.at.df<- cbind(prop.c,t(data.frame(eff.at))) 
        eff.iv.df<- cbind(prop.c,t(data.frame(eff.iv))) 
        eff.mpp.df<- cbind(prop.c,t(data.frame(eff.mpp))) 
        eff.ps.df<- cbind(prop.c,t(data.frame(eff.ps))) 
        
        return(list(eff.itt.df, eff.pp.df,eff.at.df,eff.iv.df,eff.mpp.df,eff.ps.df))
    }
    
    sim.df<-sim()
    
    itt.df<-as.data.frame(sim.df[1])
    colnames(itt.df)<- c('prop.c', 1:nIterations)
    row.names(itt.df)<- c(1:interval)
    itt.melt<- melt(itt.df, id.vars='prop.c')
    itt<-ggplot(itt.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=itt.melt$value, colour="Intention to treat", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=itt.melt$value, colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.itt)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    pp.df<-as.data.frame(sim.df[2])
    colnames(pp.df)<- c('prop.c', 1:nIterations)
    row.names(pp.df)<- c(1:interval)
    pp.melt<- melt(pp.df, id.vars='prop.c')
    pp<-ggplot(pp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=pp.melt$value, colour="Per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=pp.melt$value, colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.pp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    at.df<-as.data.frame(sim.df[3])
    colnames(at.df)<- c('prop.c', 1:nIterations)
    row.names(at.df)<- c(1:interval)
    at.melt<- melt(at.df, id.vars='prop.c')
    at<-ggplot(at.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=at.melt$value, colour="As treated", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=at.melt$value, colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.at)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    iv.df<-as.data.frame(sim.df[4])
    colnames(iv.df)<- c('prop.c', 1:nIterations)
    row.names(iv.df)<- c(1:interval)
    iv.melt<- melt(iv.df, id.vars='prop.c')
    iv<-ggplot(iv.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=iv.melt$value, colour="Instrumental variable", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=iv.melt$value, colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.iv)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    mpp.df<-as.data.frame(sim.df[5])
    colnames(mpp.df)<- c('prop.c', 1:nIterations)
    row.names(mpp.df)<- c(1:interval)
    mpp.melt<- melt(mpp.df, id.vars='prop.c')
    mpp<-ggplot(mpp.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=mpp.melt$value, colour="Modified per protocol", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=mpp.melt$value, colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.mpp)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    ps.df<-as.data.frame(sim.df[6])
    colnames(ps.df)<- c('prop.c', 1:nIterations)
    row.names(ps.df)<- c(1:interval)
    ps.melt<- melt(ps.df, id.vars='prop.c')
    ps<-ggplot(ps.melt, aes(x=prop.c, y=value))+ 
        geom_point(aes(y=ps.melt$value, colour="Propensity score matching", alpha=0.01), size=0.2) + 
        geom_smooth(aes(y=ps.melt$value, colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=col.ps)+
        theme(legend.title=element_blank(), legend.position="none")+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.2, 0.4))
    
    bias.plot<-ggarrange(itt, pp, at,
                         iv, mpp, ps,
                         ncol = 3, nrow = 2,
                         common.legend = FALSE )
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.2<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #alpha error and critical value 
    alpha=0.025
    z = qnorm(1-0.05/2)
    
    #number of data points in simulation 
    interval<-50 
    
    #make up vectors for simulations 
    .type1.error.at<-.type1.error.itt<-.type1.error.iv<-.type1.error.mpp<-.type1.error.pp<-.type1.error.ps<-c()              
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) {  
                    
                    tryCatch({ #allow the function to run in case of errors
                        
                        #simulate data frame 
                        
                        #RANDOMISATION ratio 1:1
                        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                        simdata$id<- seq(1,(2*n), by=1) #create participant id 
                        
                        #CONFOUNDER normal distribution ranging 0-1 
                        simdata$confounder<- rep(NA,2*n) 
                        simdata$confounder<- rnorm(2*n) 
                        simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #CONFOUNDER nknown normal distribution ranging 0-1 
                        simdata$confounder.u<- rep(NA,2*n) 
                        simdata$confounder.u<- rnorm(2*n) 
                        simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                        
                        #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                        simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                        simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                        
                        prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                        simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                        
                        prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                        simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                        
                        #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                        simdata$intervention <- rep(NA,2*n)
                        
                        .b0<-logit (i0)
                        b0<- seq (.b0,logit(0.01), length.out = interval)
                        .b1<-log (confounder.eff.i)
                        b1<- seq (.b1,0, length.out = interval)
                        .b2<-log (confounder.u.eff.i)
                        b2<- seq (.b2,0, length.out = interval)
                        .b3<-logit (i1)- logit (i0)
                        b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                        
                        if (b1[i]==0) { pi<- simdata$randomisation} else {
                            logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                            pi<-exp(logit.pi)/(1+exp(logit.pi))
                        }
                        simdata$intervention <- rbinom(2*n,1,prob=pi)
                        
                        #ACTUAL OUTCOMES depend on intervention
                        simdata$outcome<- rep(NA,2*n)
                        for (y in 1:(2*n)) {
                            if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                                simdata$outcome[y]<-simdata$outcome0[y]}
                        }
                        
                        #estimating type 1 error
                        ## intention to treat 
                        pz1.vector<- (filter(simdata, randomisation==1))$outcome
                        pz1.value<- mean(pz1.vector)  
                        pz0.vector<- (filter(simdata, randomisation==0))$outcome
                        pz0.value<- mean(pz0.vector)
                        eff.itt<- pz1.value-pz0.value    
                        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
                        CI.itt<- eff.itt + c(-1,1)*z*sqrt(var.eff.itt)
                        .type1.error.itt[l]<-CI.itt[2]<NImargin
                        
                        ## per protocol 
                        pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                        .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome
                        p11.value<- mean(p11.vector)
                        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome
                        p00.value<- mean(p00.vector) 
                        eff.pp <- p11.value-p00.value   
                        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
                        CI.pp<- eff.pp + c(-1,1)*z*sqrt(var.eff.pp)
                        .type1.error.pp[l]<- CI.pp[2]<NImargin
                        
                        ## as treated
                        pd1.vector<- (filter(simdata, intervention==1))$outcome
                        pd1.value<- mean(pd1.vector)
                        pd0.vector<- (filter(simdata, intervention==0))$outcome
                        pd0.value<- mean(pd0.vector)
                        eff.at <- pd1.value-pd0.value        
                        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
                        CI.at<- eff.at + c(-1,1)*z*sqrt(var.eff.at)
                        .type1.error.at[l]<- CI.at[2]<NImargin
                        
                        ## linear binomial generalized linear models with inverse probability weights
                        E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                        ps<-predict(E.out, type="response")
                        sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                        mpp<- geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)
                        eff.mpp<-mpp$coefficients['intervention']# fit linear binomial model to the weighted data
                        se<-summary(mpp)$coefficients[2,2] 
                        CI.mpp<-eff.mpp + c(-1,1)*z*se
                        .type1.error.mpp[l]<- CI.mpp[2]<NImargin
                        
                        ## propensity score matching 
                        m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                        pscore<- m$fitted.values
                        mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                        match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                        matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                        y_trt<-matched$outcome[matched$intervention==1]
                        y_con<-matched$outcome[matched$intervention==0]
                        diffy<-y_trt-y_con
                        eff.ps<-(t.test(diffy))$estimate
                        CI.psm<-(t.test(diffy))$conf.int
                        .type1.error.ps[l]<- CI.psm[2]<NImargin
                        
                        # iv with 2 stage regression
                        y<- simdata$outcome
                        x<- simdata$intervention
                        z<-simdata$randomisation
                        data<-data.frame(y,x,z)
                        asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                        eff.iv<-(summary(asmm))$ coefficients [2,1]
                        se<-(summary(asmm))$ coefficients [2,2]
                        CI.iv<-eff.iv + c(-1,1)*z*se
                        .type1.error.iv[l]<- CI.iv[2]<NImargin
                        #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                        #eff.iv<-ivmodel$coef[2] 
                        #se<-robust.se(ivmodel)[2,2]
                        #CI.iv<-eff.iv + c(-1,1)*z*se
                        # .type1.error.iv[l]<- CI.iv[2]<NImargin
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(.type1.error.itt, na.rm=TRUE)
                type1.error.pp[i] <- mean(.type1.error.pp, na.rm=TRUE)
                type1.error.at[i] <- mean(.type1.error.at, na.rm=TRUE)
                type1.error.mpp[i]<- mean(.type1.error.mpp, na.rm=TRUE)
                type1.error.iv[i] <- mean(.type1.error.iv, na.rm=TRUE)
                type1.error.ps[i] <- mean(.type1.error.ps, na.rm=TRUE)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.2<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-50
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    .b0<-logit (i0)
                    b0<- seq (.b0,logit(0.01), length.out = interval)
                    .b1<-log (confounder.eff.i)
                    b1<- seq (.b1,0, length.out = interval)
                    .b2<-log (confounder.u.eff.i)
                    b2<- seq (.b2,0, length.out = interval)
                    .b3<-logit (i1)- logit (i0)
                    b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                    
                    if (b1[i]==0) { pi<- simdata$randomisation} else {
                        logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                        pi<-exp(logit.pi)/(1+exp(logit.pi))
                    }
                    simdata$intervention <- rbinom(2*n,1,prob=pi)
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    y<- simdata$outcome
                    x<- simdata$intervention
                    z<-simdata$randomisation
                    data<-data.frame(y,x,z)
                    asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"z"], vcov="iid")
                    eff.iv[l]<-(summary(asmm))$ coefficients [2,1]
                    
                    #ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    #.eff.iv[l]<-ivmodel$coef[2]  
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.01)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.01)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.01)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.01)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.01)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.01)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 2.1: Non-compliance caused by non random process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2.11<-simdata.bias.2   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 0.5,     confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.11<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.11<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.11<-ggarrange(b2.11, 
                   ggarrange (t2.11, p2.11, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.11, top = text_grob("Case 2.11: Non-compliance caused by non random process (affect both groups)-higher value of confounder makes intervention less likely, outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.11", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

## higher value of confounder makes intervention more likely, outcome more likely
b2.12<-simdata.bias.2   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 2,     confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.12<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.12<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.12<-ggarrange(b2.12, 
                   ggarrange (t2.12, p2.12, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.12, top = text_grob("Case 2.12: Non-compliance caused by non random process (affect both groups)-higher value of confounder makes intervention more likely, outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.12", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

## higher value of confounder makes intervention less likely, outcome less likely
b2.13<-simdata.bias.2   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 0.5,     confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.13<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.13<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.13<-ggarrange(b2.13, 
                    ggarrange (t2.13, p2.13, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.13, top = text_grob("Case 2.13: Non-compliance caused by non random process (affect both groups)-higher value of confounder makes intervention less likely, outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.13", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

## higher value of confounder makes intervention more likely, outcome less likely
b2.14<-simdata.bias.2   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 2,     confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.14<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.14<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.14<-ggarrange(b2.14, 
                    ggarrange (t2.14, p2.14, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.14, top = text_grob("Case 2.14: Non-compliance caused by non random process (affect both groups)-higher value of confounder makes intervention more likely, outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.14", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

#Run simulations for case 2.2: Non-compliance caused by non random process (affect intervention groups)
b2.21<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.21<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.21<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.21<-ggarrange(b2.21, 
                    ggarrange (t2.21, p2.21, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.21, top = text_grob("Case 2.21: Non-compliance caused by non random process (affect intervention group)-higher value of confounder makes outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.21", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

b2.22<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.22<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.22<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.22<-ggarrange(b2.22, 
                   ggarrange (t2.22, p2.22, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.22, top = text_grob("Case 2.22: Non-compliance caused by non random process (affect intervention group)-higher value of confounder makes outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.22", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

#Run simulations for case 2.3: Non-compliance caused by non random process (affect control groups)
b2.31<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.31<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.31<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.31<-ggarrange(b2.31, 
                    ggarrange (t2.31, p2.31, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.31, top = text_grob("Case 2.31: Non-compliance caused by non random process (affect control group)-higher value of confounder makes outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.31", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

b2.32<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=200)
t2.32<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.32<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2.32<-ggarrange(b2.32, 
                    ggarrange (t2.32, p2.32, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.32, top = text_grob("Case 2.32: Non-compliance caused by non random process (affect control group)-higher value of confounder makes outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.32", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

#Run simulations for case 2.4: Non-compliance caused by non random process (affect control groups) with undocumented confounders
b2.41<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 10, confounder.eff.o = 0.01,  confounder.u.eff.i=100000,confounder.u.eff.o=0.2,nIterations=200)
t2.41<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 10, confounder.eff.o = 0.01, confounder.u.eff.i=100000, confounder.u.eff.o=0.2, nIterations=1000)
p2.41<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 10, confounder.eff.o = 0.01,  confounder.u.eff.i=100000,confounder.u.eff.o=0.2,nIterations=1000)
case2.41<-ggarrange(b2.41, 
                    ggarrange (t2.41, p2.41, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.41, top = text_grob("Case 2.41: Non-compliance caused by non random process (affect control group) with undocumented confounder", face = "bold", size = 12))
ggsave(filename = paste("case2.41", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

b2.42<-simdata.bias.2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2,  confounder.u.eff.i=0.01,confounder.u.eff.o=-0.1,nIterations=200)
t2.42<-simdata.type1.2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2, confounder.u.eff.i=0.01, confounder.u.eff.o=-0.1, nIterations=1000)
p2.42<- simdata.power.2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2,  confounder.u.eff.i=0.01,confounder.u.eff.o=-0.1,nIterations=1000)
case2.42<-ggarrange(b2.42, 
                    ggarrange (t2.42, p2.42, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.42, top = text_grob("Case 2.42: Non-compliance caused by non random process (affect control group) with undocumented confounder", face = "bold", size = 12))
ggsave(filename = paste("case2.42", Sys.Date(), ".pdf"), width = 13.9, height = 10.2, units = c("in"))

