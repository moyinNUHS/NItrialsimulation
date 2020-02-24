################################################################################################################
###################Using causal inference to address non-adherence in non inferiority trials####################
#################################################Set up#########################################################
################################################################################################################

# Required libraries 
library(Rcpp); library(speedglm)
library(Hmisc); library(rms); library(gsDesign)
library(dplyr); library(tableone); library(scales);
library(survey); library(gmm) #for analysis 
library(ggpubr); library(ggplot2); library(gridExtra); library(plotly) #for plots 

# source necessary codes
sourceCpp('shiny/samplesize_nonadherence/rcpp.cpp')
source('simdata_code100220.R')
source('analysis_code100220.R')
source('scenario_code100220.R')
source('plot_code100220.R')

# colour blind friendly palette 
cbPalette <- c("#009E73", "#CC79A7", "#E69F00", "#0072B2")
analysis.method <- c("Instrumental variable","Intention to treat","Inverse probability weighting","Per protocol")

# set number of data points at which simulated data is analysed  
start.interval <- 0.6
interval <- seq(from = start.interval, to = 1, by = 0.025)

getoutcome.unknownconfounding.multi<-function(vector.outcome1, vector.outcome0, intervention){
  outcome = c()
  for (i in 1:length(vector.outcome0)){
    if (intervention[i] == 1) {outcome[[i]] = vector.outcome1[i]} else {outcome[[i]] = vector.outcome0[i]}
  }
  return(unlist(outcome))
}




