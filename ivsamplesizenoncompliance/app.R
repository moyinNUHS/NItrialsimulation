###Calculation of sample size with instrumental variable in presence of non-compliance###
#########################################################################################
rm(list = ls())
require(shiny)
require(ivpack) #load ivpack package
require(gmm) #load gmm package

iv.power<- function(n, p.experiment, p.stdcare, comply.experiment, comply.stdcare, significance, NImargin){  
    
    if (significance=="1 sided 97.5%") {z <- qnorm(1-0.05/2)} else {z <- qnorm(1-0.1/2)}
    true.eff<-p.experiment-p.stdcare
    .power.iv <- c()
    nIterations<- 500 #number of iterations per simulation run
    
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
    else (print ("Non-inferiority margin must be positive, and true effect (in terms of negative outcomes) must be less than non-inferiority margin in this simulation"))
}

# Define UI 
ui <- fluidPage(

    # App title ----
    titlePanel("Power calculator accounting for non-compliance in a non-inferiority trial using instrumental variable estimation as the analysis method"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input
            numericInput(inputId = "n",
                         label = "Number of participants per group",
                         value = 100, 
                         min = 1, max = NA, step = NA),
            
            # Input
            sliderInput(inputId = "NImargin",
                        label = "Non-inferiority margin ",
                        min = 0.01,
                        max = 0.5,
                        value = 0.1,
                        step = 0.01),
            
            # Input
            sliderInput(inputId = "comply.experiment",
                        label = "Expected proportion of participants to comply in experimental arm ",
                        min = 0.5,
                        max = 1,
                        value = 1,
                        step = 0.05),
            
            # Input
            sliderInput(inputId = "comply.stdcare",
                        label = "Expected proportion of participants to comply in standard-of-care arm ",
                        min = 0.5,
                        max = 1,
                        value = 1,
                        step = 0.05),
            
            # Input
            sliderInput(inputId = "p.experiment",
                        label = "Expected proportion of participants to have outcome in experimental arm ",
                        min = 0,
                        max = 0.5,
                        value = 0,
                        step = 0.05),
            # Input
            sliderInput(inputId = "p.stdcare",
                        label = "Expected proportion of participants to have outcome in standard-of-care arm ",
                        min = 0,
                        max = 0.5,
                        value = 0,
                        step = 0.05),
            
            #input 
            selectInput(inputId="significance", 
                        label="Level of significance", 
                        choices= list("1 sided 97.5%" = "1 sided 97.5%", "1 sided 95%" = "1 sided 95%"),  
                        selected = "1 sided 97.5%", 
                        multiple = FALSE, 
                        selectize = TRUE),
            
            #calculate! 
            actionButton(inputId="run", 
                         label="Calculate!")
        ),
        
        
        # Main panel for displaying outputs ----
        mainPanel(
            tableOutput("Display"),
            
            p("Note: Power is calculated with 1000 iterations of simulated data based on entered values.")

        )
    )
)

# Define server logic required to draw a histogram ----
server <- function(input, output,session) {

    tablevalues<- eventReactive (input$run, {
        
        ivpower<-iv.power(n=input$n,
                          p.experiment=input$p.experiment,
                          p.stdcare=input$p.stdcare,
                          comply.experiment=input$comply.experiment,
                          comply.stdcare=input$comply.stdcare,
                          significance=input$significance,
                          NImargin=input$NImargin)
        
        Name = c("Number of participants per group", 
                 "Non-inferiority margin", 
                 "Proportion of participants with outcome in experimental arm", 
                 "Proportion of participants with outcome in standard-of-care arm",
                 "Proportion of participants who complied to allocated intervention in experimental arm",
                 "Proportion of participants who complied to allocated intervention in standard-of-care arm",
                 "Level of significance",
                 "Power")
        Value = as.character(c(input$n, 
                               input$NImargin, 
                               input$p.experiment, 
                               input$p.stdcare, 
                               input$comply.experiment, 
                               input$comply.stdcare, 
                               input$significance,
                               ivpower))
        cbind(Name,Value)
        
        })

    output$Display <- renderTable({
        tablevalues()
        
    })

}

shinyApp(ui = ui, server = server)

