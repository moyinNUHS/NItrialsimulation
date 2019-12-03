###Calculation of sample size in presence of non-compliance##############################
#########################################################################################
rm(list = ls())
library(shiny)
library(Rcpp)
library(gmm)
library(shinythemes)
sourceCpp(file='rcpp.cpp')

# Define UI 
ui <-fluidPage(
  theme = shinytheme("sandstone"),
    navbarPage("Power calculator accounting for non-adherence in a non-inferiority trial",
               
               tabPanel("Introduction",
                        mainPanel(
                            h3("This is a power calculator based on simulations of a two-arm non-inferiority trial with a binary outcome and time-fixed treatment."),
                            h2(),
                            h3("How to use"),
                            h4("1. Consider the potential factors which may cause non-adherence. Confounders are factors that affect both adherence to allocated treatment and the outcome.
                               Factors that affect adherence but do not affect the outcome are non-confounding. Choose the appropriate tabs above by considering the major drivers of non-adherence in the study."),
                            h4("2. Enter the number of participants per group that you would like to use to calculate power."),
                            h4("   3. Choose the non-inferiority margin."),
                            h4("   4. Choose the estimated proportion of adherent participants in each group."),
                            h4("   5. Choose the estimated proportion of participants reaching the pre-defined outcome in each group."),
                            h4("   6. Choose the level of significance."),
                            h4("   7. If the driver of non-adherence is mainly due to confounding factors, choose the confounder's direction of effect on adherence and outcome"),
                            h2(),
                            h3("Simulation mechanism"),
                            h4("Participants are randomised in a 1:1 ratio. Adherence to assigned treatment may be dependent or independent of the participants' characteristics, depending on if non-adherence is caused by confounding or non-confounding factors respectively.\n
                               Outcome is measured by the absolute risk difference between treatment failures in the experimental and control groups. 
                               Power is estimated through simulating trial data based on the alternative hypothesis with treatment effect is less than the non-inferiority margin. 
                               Iterations with the upper 95% confidence interval boundary less than the non-inferiority margin are considered to have made the correct conclusion, hence contributing to power."),
                            h2(),
                            h3("Feedback"),
                            h4("Please contact Mo Yin (moyin@tropmedres.ac) to report any issues."),
                            h2(),
                            h3("Lastest update"),
                            h4(Sys.time()),
                            width = 15
                        )
               ),
               
               tabPanel("Non-confounding factors",
                        
                        # Sidebar layout with input and output definitions ----
                        sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                                
                                # Input
                                numericInput(inputId = "nNC",
                                             label = "Number of participants per group",
                                             value = 100, 
                                             min = 1, max = NA, step = NA),
                                
                                # Input
                                sliderInput(inputId = "NImarginNC",
                                            label = "Non-inferiority margin ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0.1,
                                            step = 0.01),
                                
                                # Input
                                sliderInput(inputId = "comply.experimentNC",
                                            label = "Expected proporion of adherent participants in experimental arm ",
                                            min = 0.5,
                                            max = 1,
                                            value = 1,
                                            step = 0.05),
                                
                                # Input
                                sliderInput(inputId = "comply.stdcareNC",
                                            label = "Expected proportion of adherent participants in standard-of-care arm ",
                                            min = 0.5,
                                            max = 1,
                                            value = 1,
                                            step = 0.05),
                                
                                # Input
                                sliderInput(inputId = "p.experimentNC",
                                            label = "Expected proportion of participants to have treatment failure in experimental arm ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0,
                                            step = 0.05),
                                # Input
                                sliderInput(inputId = "p.stdcareNC",
                                            label = "Expected proportion of participants to have treatment failure in standard-of-care arm ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0,
                                            step = 0.05),
                                
                                #input 
                                selectInput(inputId="significanceNC", 
                                            label="Level of significance", 
                                            choices= list("1 sided 97.5%" = "1 sided 97.5%", "1 sided 95%" = "1 sided 95%"),  
                                            selected = "1 sided 97.5%", 
                                            multiple = FALSE, 
                                            selectize = TRUE),
                                
                                #calculate! 
                                actionButton(inputId="runNC", 
                                             label="Calculate!")
                            ),
                            
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                                tableOutput("DisplayNC"),
                                p("Note: The expected proportion of participants to have treatment failure in experimental arm should not exceed that in the standard-of-care arm by non-inferiority margin or more")
                            )
                        )), 
               
               tabPanel("Confounding factors",
                        
                        # Sidebar layout with input and output definitions ----
                        sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                                
                                # Input
                                numericInput(inputId = "nC",
                                             label = "Number of participants per group",
                                             value = 100, 
                                             min = 1, max = NA, step = NA),
                                
                                # Input
                                sliderInput(inputId = "NImarginC",
                                            label = "Non-inferiority margin ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0.1,
                                            step = 0.01),
                                
                                # Input
                                sliderInput(inputId = "comply.experimentC",
                                            label = "Expected proportion of adherent participants in experimental arm ",
                                            min = 0.5,
                                            max = 1,
                                            value = 1,
                                            step = 0.05),
                                
                                # Input
                                sliderInput(inputId = "comply.stdcareC",
                                            label = "Expected proportion of adherent participants in standard-of-care arm ",
                                            min = 0.5,
                                            max = 1,
                                            value = 1,
                                            step = 0.05),
                                
                                # Input
                                sliderInput(inputId = "p.experimentC",
                                            label = "Expected proportion of participants to have treatment failure in experimental arm ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0,
                                            step = 0.05),
                                # Input
                                sliderInput(inputId = "p.stdcareC",
                                            label = "Expected proportion of participants to have treatment failure in standard-of-care arm ",
                                            min = 0,
                                            max = 0.5,
                                            value = 0,
                                            step = 0.05),
                                
                                #input 
                                selectInput(inputId="significanceC", 
                                            label="Level of significance", 
                                            choices= list("1 sided 97.5%" = "1 sided 97.5%", "1 sided 95%" = "1 sided 95%"),  
                                            selected = "1 sided 97.5%", 
                                            multiple = FALSE, 
                                            selectize = TRUE),
                                
                                #input 
                                selectInput(inputId="confounder.intervention", 
                                            label="Effect of confounder on taking up the experimental treatment", 
                                            choices= list("Increase probability" = "Increase probability", "Decrease probability" = "Decrease probability"),  
                                            selected = "Increase probability", 
                                            multiple = FALSE, 
                                            selectize = TRUE),
                                
                                #input 
                                selectInput(inputId="confounder.outcome", 
                                            label="Effect of confounder on treatment failure", 
                                            choices= list("Increase probability" = "Increase probability", "Decrease probability" = "Decrease probability"),  
                                            selected = "Increase probability", 
                                            multiple = FALSE, 
                                            selectize = TRUE),
                                
                                #calculate! 
                                actionButton(inputId="runC", 
                                             label="Calculate!")
                            ),
                            
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                                tableOutput("DisplayC"),
                                p("Note: The expected proportion of participants to have treatment failure in experimental arm should not exceed that in the standard-of-care arm by non-inferiority margin or more")
                            )
                        ))
    ))

# Define server logic required ----
server <- function(input, output,session) {
    
    tablevaluesNC<- eventReactive (input$runNC, {
        
        n=input$nNC;
        p.experiment=input$p.experimentNC;
        p.stdcare=input$p.stdcareNC;
        comply.experiment=input$comply.experimentNC;
        comply.stdcare=input$comply.stdcareNC;
        significance=input$significanceNC;
        NImargin=input$NImarginNC
        
        withProgress(message='Calculating Power', value=0,{
            
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
                    confounder = rep(rbeta(n=n,shape1=2,shape2=2),2) #confounder beta distribution ranging 0-1 
                    
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
                    weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
                    outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
                    eff.mpp=coef(outcomemodel)[2]
                    se<-sqrt(diag(vcovHC(outcomemodel)))[2]
                    CI.mpp<-eff.mpp+z*se
                    mpp<- CI.mpp<NImargin
                    
                    # iv with 2 stage regression
                    asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
                    eff.iv=summary(asmm)$ coefficients [2,1]
                    se<-summary(asmm)$coefficients [2,2]
                    CI.iv<-eff.iv + z*se
                    iv<- CI.iv<NImargin
                    
                    power.iter[[l]]<-c(itt, ppp, mpp, iv)
                    
                    # Increment the progress bar, and update the detail text.
                    incProgress(1/nIterations, detail = paste("Iteration number (out of 1000 iterations)", l))
                    
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
            }
            
            # mean of power from iterated data 
            power.matrix=matrix(as.numeric(unlist(power.iter)), nrow=nIterations, ncol=4,byrow = TRUE)
            power.cal= colMeans(power.matrix, na.rm = TRUE)
        })
        
        Name = c("Number of participants per group", 
                 "Non-inferiority margin", 
                 "Proportion of participants with outcome in experimental arm", 
                 "Proportion of participants with outcome in standard-of-care arm",
                 "Proportion of participants who complied to allocated treatment in experimental arm",
                 "Proportion of participants who complied to allocated treatment in standard-of-care arm",
                 "Level of significance",
                 "Power using intention to treat analysis",
                 "Power using per-protocol analysis",
                 "Power using inverse probability weighting analysis",
                 "Power using instrumental varible analysis")
        
        Value = as.character(c(input$nNC, 
                               input$NImarginNC, 
                               input$p.experimentNC, 
                               input$p.stdcareNC, 
                               input$comply.experimentNC, 
                               input$comply.stdcareNC, 
                               input$significanceNC,
                               power.cal[1],
                               power.cal[2],
                               power.cal[3],
                               power.cal[4]))
        cbind(Name,Value)
        
    })
    
    output$DisplayNC <- renderTable({
        tablevaluesNC()
        
    })
    
    tablevaluesC<- eventReactive (input$runC, { 
        
        n=input$nC;
        p.experiment=input$p.experimentC;
        p.stdcare=input$p.stdcareC;
        comply.experiment=input$comply.experimentC;
        comply.stdcare=input$comply.stdcareC;
        significance=input$significanceC;
        confounder.intervention=input$confounder.intervention;
        confounder.outcome=input$confounder.outcome
        NImargin=input$NImarginC
        
        withProgress(message='Calculating Power', value=0,{
            
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
                    if (confounder.outcome=="Increase probability") {
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
                    if (confounder.intervention=="Increase probability") {
                        shape2<-runif(1, min=1.5, max=5)
                        shape1<- shape2*comply.experiment/(1-comply.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
                        comply.experiment.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
                        int.intervention<- rbinom(n, 1, prob=comply.experiment.ind) #increasing confounder value will have decreasing probability for outcome 
                        
                        shape1<- shape2*(1-comply.stdcare)/(1-(1-comply.stdcare)) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
                        comply.stdcare.ind<-sort(rbeta(n=n, shape1 = shape1, shape2 = shape2)) #individual probability with mean of p.experiment, in increasing order
                        cont.intervention<- rbinom(n, 1, prob=comply.stdcare.ind)
                        
                    } else {
                        shape2<-runif(1, min=1.5, max=5)
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
                    weight= pp$intervention*mean(pp$intervention)/score+(1-pp$intervention)*(1-mean(pp$intervention))/(1-score)#create stabilized weights, using a null model with intervention as the dependent variable
                    outcomemodel=glm(outcome~intervention, family=binomial(link="identity"), weights=weight, data=pp) #identity link for risk difference
                    eff.mpp=coef(outcomemodel)[2]
                    se<-sqrt(diag(vcovHC(outcomemodel,type="HC0")))[2]
                    CI.mpp<-eff.mpp+z*se
                    mpp<- CI.mpp<NImargin
                    
                    # iv with 2 stage regression
                    asmm <- gmm(simdata[,5] ~ simdata[,4], x=simdata[,2], vcov="iid")
                    eff.iv=summary(asmm)$ coefficients [2,1]
                    se<-summary(asmm)$coefficients [2,2]
                    CI.iv<-eff.iv + z*se
                    iv<- CI.iv<NImargin
                    
                    power.iter[[l]]<-c(itt,ppp, mpp,iv)
                    
                    # Increment the progress bar, and update the detail text.
                    incProgress(1/nIterations, detail = paste("Iteration number (out of 1000 iterations)", l))
                    
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
            }
            
            # mean of power from iterated data 
            power.matrix=matrix(as.numeric(unlist(power.iter)), nrow=nIterations, ncol=4,byrow = TRUE)
            power.cal= colMeans(power.matrix, na.rm = TRUE)
        })
        
        Name = c("Number of participants per group", 
                 "Non-inferiority margin", 
                 "Proportion of participants with outcome in experimental arm", 
                 "Proportion of participants with outcome in standard-of-care arm",
                 "Proportion of participants who complied to allocated treatment in experimental arm",
                 "Proportion of participants who complied to allocated treatment in standard-of-care arm",
                 "Level of significance",
                 "Effect of confounder on taking up the experimental treatment",
                 "Effect of confounder on treatment failure",
                 "Power using intention to treat analysis",
                 "Power using per-protocol analysis",
                 "Power using inverse probability weighting analysis",
                 "Power using instrumental varible analysis")
        
        Value = as.character(c(input$nC, 
                               input$NImarginC, 
                               input$p.experimentC, 
                               input$p.stdcareC, 
                               input$comply.experimentC, 
                               input$comply.stdcareC, 
                               input$significanceC,
                               input$confounder.intervention,
                               input$confounder.outcome,
                               power.cal[1],
                               power.cal[2],
                               power.cal[3],
                               power.cal[4]))
        cbind(Name,Value)
        
    })
    
    output$DisplayC <- renderTable({
        tablevaluesC()
    })
    
}

shinyApp(ui = ui, server = server)

