###Calculation of sample size with instrumental variable in presence of non-compliance###
#########################################################################################


rm(list=ls()) # Clean working environment 
require(ivpack) #load ivpack package
require(gmm) #load gmm package

#Define parameters 
n<- 200                 # n:                  number of participants per group
NImargin<- 0.2          # NImargin:           non inferiority margin 
p.experiment<-0.4       # p.experiment:       proportion of participants with outcome in experimental arm 
p.stdcare<- 0.4         # p.stdcare:          proportion of participants with outcome in standard-of-care arm 
comply.experiment<-0.8  # comply.experiment:  proportion of participants who complied to allocated intervention in experimental arm 
comply.stdcare<- 0.9    # comply.stdcare:     proportion of participants who complied to allocated intervention in standard-of-care arm
significance<-1         # significance:       level of significance / alpha 

#Simulate a two-arm randomised controlled non-inferiority trial 

#alpha error and critical value 
iv.power<- function(n, p.experiment, p.stdcare, comply.experiment, comply.stdcare, significance, NImargin){  
    
    if (significance=="1 sided 97.5%") {z <- qnorm(1-0.05/2)} else {z <- qnorm(1-0.1/2)}
    true.eff<-p.experiment-p.stdcare
    .power.iv <- c()
    nIterations<- 50 #number of iterations per simulation run
    
    if ((NImargin>0) & (true.eff < NImargin)) {
        
        for(l in 1:nIterations) { 
            tryCatch({ #allow the function to run in case of errors
                
                #simulate trial data frame 
                #RANDOMISATION ratio 1:1
                simdata <- data.frame(seq(1,(2*n), by=1)) #create participant ID 
                simdata$randomization<-  c(rep(1,n), rep(0,n)) # randomization in 1:1 ratio 
                
                #INTERVENTION - according to predefined compliance in each group 
                n.comply.experiment<- round(comply.experiment*n)
                n.nocomply.experiment<- n-n.comply.experiment
                n.comply.stdcare<- round(comply.stdcare*n)
                n.nocomply.stdcare<- n-n.comply.stdcare
                simdata$intervention <- c(rep(1,n.comply.experiment), rep(0, n.nocomply.experiment), rep(0,n.comply.stdcare), rep(1,n.nocomply.stdcare))
                
                #OUTCOMES counterfactual, depend on actual intervention
                simdata$outcome1 <- rbinom(2*n,1,prob=p.experiment) #outcome for all participants if they take up intervention 
                simdata$outcome0 <- rbinom(2*n,1,prob=p.stdcare) #outcome for all participants if they take up standard of care
                
                for (k in 1:(2*n)) {
                    if (simdata$intervention[k]==1) {simdata$outcome[k]<-simdata$outcome1[k] } else 
                        simdata$outcome[k]<-simdata$outcome0[k] 
                } 
                
                # IV with 2 stage regression
                y<- simdata$outcome
                x<- simdata$intervention
                iv<-simdata$randomization
                data<-data.frame(y,x,iv)
                asmm <- gmm(data[,"y"] ~ data[,"x"], x=data[,"iv"], vcov="iid")
                eff.iv<-(summary(asmm))$ coefficients [2,1]
                se<-(summary(asmm))$ coefficients [2,2]
                CI.iv<-eff.iv + c(-1,1)*z*se
                .power.iv[l]<- CI.iv[2]<NImargin
                
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
        }
        # mean of power from iterated data 
        power.iv<- mean(.power.iv)
        return (power.iv)
    }
    else (print ("NI margin must be positive, and true effect (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

iv.power(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, comply.experiment=comply.experiment, comply.stdcare=comply.stdcare,significance=significance, NImargin=NImargin) 

