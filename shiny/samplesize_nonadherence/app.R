###Calculation of sample size in presence of non-compliance###
##############################################################
rm(list = ls())

#load required libraries
library(shiny); library(shinythemes); library(shinyjs)
library(Rcpp)
library(gmm); library(speedglm)

#source supporting codes
sourceCpp(file = 'rcpp.cpp')
source(file = 'analysisestimate.R')
source(file = 'analysispower.R')
eff.conf.df = read.csv(file = 'eff.conf.outcome.tab.csv')[,-1]

# Define UI
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  theme = shinytheme("sandstone"),
  
  navbarPage(
    title = "Power calculator accounting for non-adherence in a non-inferiority trial",
    
    tabPanel(
      "Introduction",
      mainPanel(
        h4(
          "This is a power calculator based on simulations of a two-arm non-inferiority trial with a binary outcome and time-fixed treatment."
        ),
        h2(),
        h4("How to use"),
        h5(
          "1. Consider the potential factors which may cause non-adherence. Confounders are factors that affect both adherence to allocated treatment and the outcome.
          Factors that affect adherence but do not affect the outcome are non-confounding. Choose the appropriate tabs above by considering the major drivers of non-adherence in the study."
        ),
        h5(
          "2. Enter the number of participants per arm that you would like to use to calculate power."
        ),
        h5(
          "   3. Choose the non-inferiority margin. This is the absolute difference between probability of treatment failure in experimental and standard-of-care arms below which the experimental treatment is considered non-inferior."
        ),
        h5(
          "   4. Choose the actual treatment that non-adherent participants are likely to take up. Non-adherent participants may either cross-over to the opposite arm, or take up an alternative treatment with a different treatment efficacy."
        ),
        h5(
          "   5. Choose the estimated proportion of adherent participants in each arm."
        ),
        h5(
          "   6. Choose the estimated proportion of participants with treatment failure in each arm."
        ),
        h5("  7. Choose the level of significance."
        ),
        h5(
          "   8. If the driver of non-adherence is mainly due to confounding factors," 
        ),
        h5( "i) choose the confounder's direction of effect on adherence and outcome; "
        ),
        h5(" ii) choose the confounder's magnitude of effect on outcome (estimated with treatment failure as the dependent variable and confounder as the independent variable)."
        ),
        h2(),
        h4("Simulation mechanism"),
        h5(
          "Participants are randomised in a 1:1 ratio. Adherence to assigned treatment may be dependent or independent of the participants' characteristics, depending on if non-adherence is caused by confounding or non-confounding factors respectively."
        ),
        h5("The treatment effect is measured by the absolute risk difference between treatment failures in the experimental and control arms (so lower proportion of treatment failure in the experimental arm correspond to lower treatment effect)."
        ),
        h5("Power is estimated through simulating trial data based on the alternative hypothesis that the treatment effect is less than the non-inferiority margin.
          Iterations with the upper 95% confidence interval boundary less than the non-inferiority margin are considered to have made the correct conclusion, hence contributing to power."
        ),
        h4(),
        h5(
          "Further details of the simulation mechanism are shared in our publication, Mo Y, Lim C, Mukaka M and Cooper BS. Statistical considerations in the design and analysis of non-inferiority trials with binary endpoints in the presence of non-adherence: a simulation study. Wellcome Open Res 2019, 4:207 (https://doi.org/10.12688/wellcomeopenres.15636.1)."
        ),
        h2(),
        h4("Feedback"),
        h5(
          "Please contact Mo Yin (moyin@tropmedres.ac) to report any issues."
        ),
        h2(),
        h4("Lastest update"),
        h5("25th February, 2020"),
        width = 15
      )
    ),
    
    tabPanel(
      "Non-confounding factors",
      
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
          # Input
          numericInput(
            inputId = "nNC",
            label = "Number of participants per arm",
            value = 100,
            min = 1,
            max = NA,
            step = NA
          ),
          
          tags$head(tags$style(HTML('.js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          
          # Input
          sliderInput(
            inputId = "NImarginNC",
            label = "Non-inferiority margin (absolute scale)",
            min = 0.01,
            max = 0.3,
            value = 0.1,
            step = 0.01
          ),
          
          # Input
          radioButtons(
            inputId = "nonadhere.txNC",
            label = "Actual treatment likely to be taken up by the non-adherent participants",
            choices = list("Cross-over to the opposite arm" = 1, 
                           "Alternative treatment (e.g. default all treatments)" = 2),
            selected = 1
          ),
          
          # Input
          sliderInput(
            inputId = "adhere.experimentNC",
            label = "Expected proportion of adherent participants in experimental arm ",
            min = 0.5,
            max = 1,
            value = 1,
            step = 0.05
          ),
          
          # Input
          sliderInput(
            inputId = "adhere.stdcareNC",
            label = "Expected proportion of adherent participants in standard-of-care arm ",
            min = 0.5,
            max = 1,
            value = 1,
            step = 0.05
          ),
          
          # Input
          sliderInput(
            inputId = "p.experimentNC",
            label = "Expected proportion of participants to have treatment failure in the experimental arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          #Input
          sliderInput(
            inputId = "p.stdcareNC",
            label = "Expected proportion of participants to have treatment failure in the standard-of-care arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          #Input
          sliderInput(
            inputId = "p.altNC",
            label = "Expected proportion of participants to have treatment failure in the alternative treatment arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          #input
          selectInput(
            inputId = "significanceNC",
            label = "Level of significance",
            choices = list("1 sided 97.5%" = "1 sided 97.5%", "1 sided 95%" = "1 sided 95%"),
            selected = "1 sided 97.5%",
            multiple = FALSE,
            selectize = TRUE
          ),
          
          #calculate!
          actionButton(inputId = "runNC",
                       label = "Calculate!")
        ),
        
        
        # Main panel for displaying outputs ----
        mainPanel(
          tableOutput("DisplayNC"),
          p(
            "Note: proportion of participants with treatment failure in the experimental arm has been constrained so as not to exceed that in the standard-of-care arm by the non-inferiority margin or more"
          )
        )
      )
    ),
    
    tabPanel(
      "Confounding factors",
      
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
          # Input
          numericInput(
            inputId = "nC",
            label = "Number of participants per arm",
            value = 100,
            min = 1,
            max = NA,
            step = NA
          ),
          
          tags$head(tags$style(HTML('.js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-7 .irs-single, .js-irs-7 .irs-bar-edge, .js-irs-7 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-8 .irs-single, .js-irs-8 .irs-bar-edge, .js-irs-8 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-9 .irs-single, .js-irs-9 .irs-bar-edge, .js-irs-9 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-10 .irs-single, .js-irs-10 .irs-bar-edge, .js-irs-10 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-11 .irs-single, .js-irs-11 .irs-bar-edge, .js-irs-11 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          tags$head(tags$style(HTML('.js-irs-12 .irs-single, .js-irs-12 .irs-bar-edge, .js-irs-12 .irs-bar {
                                                  background: #40798C;
                                                  border-top: 1px solid #40798C;
                                                  border-bottom: 1px solid #40798C;}'))),
          
          # Input
          sliderInput(
            inputId = "NImarginC",
            label = "Non-inferiority margin (absolute scale)",
            min = 0.01,
            max = 0.3,
            value = 0.1,
            step = 0.01
          ),
          
          # Input
          radioButtons(
            inputId = "nonadhere.txC",
            label = "Actual treatment likely to be taken up by the non-adherent participants",
            choices = list("Cross-over to the opposite arm" = 1, 
                           "Alternative treatment (e.g. default all treatments)" = 2),
            selected = 1
          ),
          
          # Input
          sliderInput(
            inputId = "adhere.experimentC",
            label = "Expected proportion of adherent participants in experimental arm ",
            min = 0.5,
            max = 1,
            value = 1,
            step = 0.05
          ),
          
          # Input
          sliderInput(
            inputId = "adhere.stdcareC",
            label = "Expected proportion of adherent participants in standard-of-care arm ",
            min = 0.5,
            max = 1,
            value = 1,
            step = 0.05
          ),
          
          # Input
          sliderInput(
            inputId = "p.experimentC",
            label = "Expected proportion of participants to have treatment failure in the experimental arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          # Input
          sliderInput(
            inputId = "p.stdcareC",
            label = "Expected proportion of participants to have treatment failure in the standard-of-care arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          # Input
          sliderInput(
            inputId = "p.altC",
            label = "Expected proportion of participants to have treatment failure in the alternative treatment arm ",
            min = 0,
            max = 0.5,
            value = 0.1,
            step = 0.01
          ),
          
          #input
          selectInput(
            inputId = "significanceC",
            label = "Level of significance",
            choices = list("1 sided 97.5%" = "1 sided 97.5%", "1 sided 95%" = "1 sided 95%"),
            selected = "1 sided 97.5%",
            multiple = FALSE,
            selectize = TRUE
          ),
          
          #input
          selectInput(
            inputId = "confounder.intervention",
            label = "Effect of confounder on taking up the experimental treatment",
            choices = list(
              "Participants with high confounder values from the standard-of-care arm tend to be non-adherent; 
              participants with low confounder values from from the experimental arm tend to be non-adherent" = "Increase probability",
              "Participants with low confounder values from the standard-of-care arm tend to be non-adherent; 
              participants with high confounder values from from the experimental arm tend to be non-adherent" = "Decrease probability"
            ),
            selected = "Participants with high confounder values from the standard-of-care arm tend to be non-adherent; 
              participants with low confounder values from from the experimental arm tend to be non-adherent",
            multiple = FALSE,
            selectize = TRUE
          ),
          
          #input
          selectInput(
            inputId = "confounder.outcome",
            label = "Effect of confounder on treatment failure",
            choices = list(
              "Increase probability" = "Increase probability",
              "Decrease probability" = "Decrease probability"
            ),
            selected = "Increase probability",
            multiple = FALSE,
            selectize = TRUE
          ),
          
          #input
          sliderInput(
            inputId = "confounder.eff",
            label = "Magnitude of direct confounding effect on treatment failure 
            (when estimated using treatment failure as the dependent variable, and confounder as the independent variable in a linear regression)",
            min = 0.5,
            max = 10,
            value = 1,
            step = 0.5
          ),
          
          #calculate!
          actionButton(inputId = "runC",
                       label = "Calculate!")
        ),
        
        
        # Main panel for displaying outputs ----
        mainPanel(
          tableOutput("DisplayC"),
          p(
            "Note: proportion of participants with treatment failure in the experimental arm has been constrained so as not to exceed that in the standard-of-care arm by the non-inferiority margin or more"
          )
        )
      )
    )
  )
)

# Define server logic required ----
server <- function(input, output, session) {
  
  #make dynamic slider
  observe({
    # Control the value, min, max, and step
    updateSliderInput(
      session,
      "p.stdcareNC",
      value = input$p.experimentNC,
      min = ifelse (
        input$p.experimentNC - input$NImarginNC + 0.01 <= 0, 0,
        input$p.experimentNC - input$NImarginNC + 0.01),
      max = input$p.experimentNC + 0.25,
      step = 0.01
    )
    
    updateSliderInput(
      session,
      "p.stdcareC",
      value = input$p.experimentC,
      min = ifelse (
        input$p.experimentC - input$NImarginC + 0.01 <= 0, 0,
        input$p.experimentC - input$NImarginC + 0.01),
      max = input$p.experimentC + 0.25,
      step = 0.01
    )
    
    updateSliderInput(
      session,
      "confounder.eff",
      value = input$confounder.outcome,
      min = ifelse (input$confounder.outcome == "Increase probability", 
                    0.5, -10),
      max = ifelse (input$confounder.outcome == "Increase probability", 
                    10, -0.5),
      step = 0.5
    )
    
  })
  
  observeEvent(input$nonadhere.txNC, {
    if(input$nonadhere.txNC == 1){
      shinyjs::disable("p.altNC")
    }else{
      shinyjs::enable("p.altNC")
    }})
  
  observeEvent(input$nonadhere.txC, {
    if(input$nonadhere.txC == 1){
      shinyjs::disable("p.altC")
    }else{
      shinyjs::enable("p.altC")
    }})
  
  tablevaluesNC <- eventReactive(input$runNC, {
    n = input$nNC
    NImargin = input$NImarginNC
    nonadhere.tx = input$nonadhere.txNC
    p.experiment = input$p.experimentNC
    p.stdcare = input$p.stdcareNC
    p.alt = input$p.altNC
    adhere.experiment = input$adhere.experimentNC
    adhere.stdcare = input$adhere.stdcareNC
    significance = input$significanceNC
    
    withProgress(message = 'Calculating Power', value = 0, {
      
      #make up vectors for simulations
      estimate.iter = c()
      
      #make dummy outcome2 if crossover
      if (nonadhere.tx == 1) { p.alt = 0 }
      
      #number of iterations
      nIterations = 1000
      
      #significance
      ifelse (significance == "1 sided 97.5%",
              z <- qnorm(0.975),
              z <- qnorm(0.95))
      
      id = 1: (2 * n) #create participant id
      randomisation = sample(rep(0:1, n)) #randomisation
      confounder = rbeta(n = 2 * n, shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1
      
      #simulate and derive treatment effect
      for (l in 1:nIterations) {
        tryCatch({
          #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
          outcome0 = rbinom(2 * n, prob = p.stdcare, size = 1)
          outcome1 = rbinom(2 * n, prob = p.experiment, size = 1)
          outcome2 = rbinom(2 * n, prob = p.alt, size = 1)
          
          #INTERVENTION dependent on adherence
          intervention = rep(NA, 2 * n)
          intervention[sample(which(randomisation == 1), size = adhere.experiment * n)] = 1
          intervention[sample(which(randomisation == 0), size = adhere.stdcare * n)] = 0
          
          if (nonadhere.tx == 1){
            intervention[which(randomisation == 1) %in% which(is.na(intervention))] = 0
            intervention[which(randomisation == 0) %in% which(is.na(intervention))] = 1
          } else {
            intervention[which(is.na(intervention))] = 2
          }
          
          #ACTUAL OUTCOMES depend on intervention
          outcome = getoutcome(outcome0, outcome1, outcome2, intervention)
          
          simdata = matrix(data = c(id, randomisation, confounder, intervention, outcome),
                           nrow = (2 * n))
          
          estimate.iter[[l]] = analysis.estimate(simdata = simdata)
          
          # Increment the progress bar, and update the detail text.
          incProgress(1 / nIterations,
                      detail = paste("Iteration number (out of 1000 iterations)", l))
          
        }, error = function(e) {
          cat("ERROR :", conditionMessage(e), "\n")
        }) #receive error message if there is an error
      }
      
      # mean of power from iterated data
      power.calNC = analysis.power(estimate.iter = estimate.iter, z = z, NImargin = NImargin, nIterations = nIterations)
    })
    
    Name = c(
      "Number of participants per arm",
      "Non-inferiority margin (absolute scale)",
      "Actual treatment likely to be taken up by the non-adherent participants",
      "Proportion of participants with treatment failure in the experimental arm",
      "Proportion of participants with treatment failure in the standard-of-care arm",
      "Proportion of participants with treatment failure in the alternative treatment arm",
      "Proportion of adherent participants in the experimental arm",
      "Proportion of adherent participants in the standard-of-care arm",
      "Level of significance",
      "Power using intention to treat analysis",
      "Power using per-protocol analysis",
      "Power using inverse probability weighting estimation",
      "Power using instrumental varible analysis"
    )
    
    Value = as.character(
      c(
        input$nNC,
        input$NImarginNC,
        input$nonadhere.txNC,
        input$p.experimentNC,
        input$p.stdcareNC,
        input$p.altNC,
        input$adhere.experimentNC,
        input$adhere.stdcareNC,
        input$significanceNC,
        round(power.calNC[1], 3),
        round(power.calNC[2], 3),
        round(power.calNC[3], 3),
        round(power.calNC[4], 3)
      )
    )
    
    t = cbind(Name, Value)
    if (nonadhere.tx != 1) {t = t[-nrow(t),]}
    t
  })
  
  output$DisplayNC <- renderTable({
    tablevaluesNC()},
    hover = TRUE, spacing = 'm')
  
  tablevaluesC <- eventReactive (input$runC, {
    n = input$nC
    NImargin = input$NImarginC
    nonadhere.tx = input$nonadhere.txC
    p.experiment = input$p.experimentC
    p.stdcare = input$p.stdcareC
    p.alt = input$p.altC
    adhere.experiment = input$adhere.experimentC
    adhere.stdcare = input$adhere.stdcareC
    significance = input$significanceC
    confounder.eff = input$confounder.eff
    confounder.intervention = input$confounder.intervention
    confounder.outcome = input$confounder.outcome
    
    withProgress(message = 'Calculating Power', value = 0, {
      
      #make up vectors for simulations
      estimate.iter = eff.conf.outcome = c()
      
      #make dummy outcome2 if crossover
      if (nonadhere.tx == 1) { p.alt = 0 }
      
      #number of iterations
      nIterations = 1000
      
      #significance
      ifelse (significance == "1 sided 97.5%",
              z <- qnorm(0.975),
              z <- qnorm(0.95))
      
      id = 1:(2 * n) #create participant id
      randomisation = sample(rep(0:1, n)) #randomisation
      
      #degree of effect of confounder on outcome
      p.experiment.rnd = round(p.experiment/0.05) * 0.05
      conf.df.inp = eff.conf.df[which(eff.conf.df$coef.lower == (abs(confounder.eff) - 0.25)),]
      conf.df.inpt = conf.df.inp[which(as.character(conf.df.inp$p) == as.character(p.experiment.rnd)),]
      shape2min = conf.df.inpt$shape.min
      shape2max = conf.df.inpt$shape.max
      
      #simulate and derive treatment effect
      for (l in 1:nIterations) {
        tryCatch({
          confounder = rbeta(n = 2 * n, shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1
          
          #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
          shape2 = runif(1, min = shape2min, max = shape2max)
          if (confounder.outcome == "Increase probability") {
            #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome
            shape1 = shape2 * p.experiment / (1 - p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            p.experiment.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2 ) #individual probability with mean of p.experiment, in increasing order
            p.experiment.ind = p.experiment.ind.gen[order(p.experiment.ind.gen)]
            
            outcome1 = rbinom(2 * n, 1, prob = p.experiment.ind) #increasing confounder value will have increasing probability for outcome
            
            shape1 = shape2 * p.stdcare / (1 - p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            p.stdcare.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2) #individual probability with mean of p.stdcare, in increasing order
            p.stdcare.ind = p.stdcare.ind.gen[order(p.stdcare.ind.gen)]
            
            outcome0 = rbinom(2 * n, 1, prob = p.stdcare.ind) #increasing confounder value will have increasing probability for outcome
            
            if (nonadhere.tx != 1) {
              shape1 = shape2 * p.alt / (1 - p.alt) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
              p.alt.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2) #individual probability with mean of p.alt, in increasing order
              p.alt.ind = p.alt.ind.gen[order(p.alt.ind.gen)]
              
              outcome2 = rbinom(2 * n, 1, prob = p.alt.ind) #increasing confounder value will have increasing probability for outcome
            } else {
              outcome2 = rbinom(2 * n, prob = p.alt, size = 1)
            }
            
          } else {
            shape1 = shape2 * p.experiment / (1 - p.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            p.experiment.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2) #individual probability with mean of p.experiment, in decreasing order
            p.experiment.ind = p.experiment.ind.gen[order(p.experiment.ind.gen, decreasing = TRUE)]
            
            outcome1 = rbinom(2 * n, 1, prob = p.experiment.ind) #increasing confounder value will have decreasing probability for outcome
            
            shape1 = shape2 * p.stdcare / (1 - p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            p.stdcare.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2) #individual probability with mean of p.experiment, in increasing order
            p.stdcare.ind = p.stdcare.ind.gen[order(p.stdcare.ind.gen, decreasing = TRUE)]
            
            outcome0 = rbinom(2 * n, 1, prob = p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
            
            if (nonadhere.tx != 1) {
              shape1 = shape2 * p.alt / (1 - p.alt) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
              p.alt.ind.gen = rbeta(n = (2 * n), shape1 = shape1, shape2 = shape2) #individual probability with mean of p.alt, in increasing order
              p.alt.ind = p.alt.ind.gen[order(p.alt.ind.gen, decreasing = TRUE)]
              
              outcome2 = rbinom(2 * n, 1, prob = p.alt.ind) #increasing confounder value will have decreasing probability for outcome
            } else {
              outcome2 = rbinom(2 * n, prob = p.alt, size = 1)
            }
          }
          
          d = matrix(data = c(id, randomisation, confounder),
                     nrow = (2 * n))
          d.ordered = matrix(cbind(d[order(d[, 3]), ], outcome0, outcome1, outcome2), ncol = 6) #order confounder in ascending order
          d.grouped = rbind(d.ordered[which(d.ordered[, 2] == 1), ], d.ordered[which(d.ordered[, 2] == 0), ])
          
          #INTERVENTION dependent on randomisation and confounders
          shape2 = runif(1, min = 2, max = 10)
          if (confounder.intervention == "Increase probability") {
            shape1 = shape2 * adhere.experiment / (1 - adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            adhere.experiment.ind.gen = rbeta(n = n, shape1 = shape1, shape2 = shape2) #individual probability with mean of p.experiment, in increasing order
            adhere.experiment.ind = adhere.experiment.ind.gen[order(adhere.experiment.ind.gen)]
            int.intervention = rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
            
            if (nonadhere.tx != 1) {
              int.intervention[int.intervention == 0] = 2
            }
            
            shape1 = shape2 * (1 - adhere.stdcare) / (1 - (1 - adhere.stdcare)) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            adhere.stdcare.ind.gen = rbeta(n = n, shape1 = shape1, shape2 = shape2) #individual probability with mean of p.experiment, in increasing order
            adhere.stdcare.ind = adhere.stdcare.ind.gen[order(adhere.stdcare.ind.gen)]
            cont.intervention = rbinom(n, 1, prob = adhere.stdcare.ind)
            
            if (nonadhere.tx != 1) {
              cont.intervention[cont.intervention == 1] = 2
            }
            
          } else {
            shape1 = shape2 * adhere.experiment / (1 - adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            adhere.experiment.ind.gen = rbeta(n = n, shape1 = shape1, shape2 = shape2)
            adhere.experiment.ind = adhere.experiment.ind.gen[order(adhere.experiment.ind.gen, decreasing = TRUE)] #individual probability with mean of p.experiment, in decreasing order
            int.intervention = rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome
            
            if (nonadhere.tx != 1) {
              int.intervention[int.intervention == 0] = 2
            }
            
            shape1 = shape2 * adhere.stdcare / (1 - adhere.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
            adhere.stdcare.ind.gen = rbeta(n = n, shape1 = shape1, shape2 = shape2) #individual probability with mean of p.experiment, in decreasing order
            adhere.stdcare.ind = adhere.stdcare.ind.gen[order(adhere.stdcare.ind.gen, decreasing = TRUE)]
            cont.intervention = rbinom(n, 1, prob = 1 - adhere.stdcare.ind)
            
            if (nonadhere.tx != 1) {
              cont.intervention[cont.intervention == 1] = 2
            }
          }
          
          intervention <- c(int.intervention, cont.intervention)
          
          #ACTUAL OUTCOMES depend on intervention
          outcome = getoutcome(d.grouped[, 4], d.grouped[, 5], d.grouped[, 6], intervention)
          
          simdata = matrix(cbind(d.grouped[, c(-6:-4)], intervention, outcome), ncol = 5)
          
          estimate.iter[[l]] = analysis.estimate(simdata = simdata)
          
          # Report effect of confounder on probability of outcome
          model = speedglm(V4 ~ V3, family = binomial(link = "logit"), data = as.data.frame(d.ordered)) #V3 = confounder, V4 = counterfactual outcome
          eff.conf.outcome[[l]] = coefficients(model)[2]
          
          # Increment the progress bar, and update the detail text.
          incProgress(1 / nIterations,
                      detail = paste("Iteration number (out of 1000 iterations)", l))
          
        }, error = function(e) {
          cat("ERROR :", conditionMessage(e), "\n")
        }) #receive error message if there is an error
      }
      
      # mean of power from iterated data
      power.calC = analysis.power(estimate.iter = estimate.iter, z = z, NImargin = NImargin, nIterations = nIterations)
      
      #mean of effect of confounder on outcome 
      eff.conf.outcomeC = sum(unlist(eff.conf.outcome))/nIterations
    })
    
    Name = c(
      "Number of participants per arm",
      "Non-inferiority margin (absolute scale)",
      "Actual treatment likely to be taken up by the non-adherent participants",
      "Proportion of participants with treatment failure in the experimental arm",
      "Proportion of participants with treatment failure in the standard-of-care arm",
      "Proportion of participants with treatment failure in the alternative treatment arm",
      "Proportion of adherent participants in the experimental arm",
      "Proportion of adherent participants in the standard-of-care arm",
      "Level of significance",
      "Effect of confounder on taking up the experimental treatment",
      "Effect of confounder on treatment failure",
      paste("Magnitude of direct confounding effect on treatment failure", "(estimated using treatment failure as the dependent variable, and confounder as the independent variable in a linear regression)", sep = '<br>'),
      "Power using intention to treat analysis",
      "Power using per-protocol analysis",
      "Power using inverse probability weighting estimation",
      "Power using instrumental varible analysis"
    )
    
    Value = as.character(
      c(
        input$nC,
        input$NImarginC,
        input$nonadhere.txC,
        input$p.experimentC,
        input$p.stdcareC,
        input$p.altC,
        input$adhere.experimentC,
        input$adhere.stdcareC,
        input$significanceC,
        input$confounder.intervention,
        input$confounder.outcome,
        format(round(eff.conf.outcomeC, 1), nsmall = 1),
        round(power.calC[1], 3),
        round(power.calC[2], 3),
        round(power.calC[3], 3),
        round(power.calC[4], 3)
      )
    )
    t = cbind(Name, Value)
    if (nonadhere.tx != 1) {t = t[-nrow(t),]}
    t
  })
  
  output$DisplayC <- renderTable({
    tablevaluesC()}, 
    hover = TRUE, spacing = 'm', 
    sanitize.text.function=identity
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)