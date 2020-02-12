################################################################################################################
###################Using causal inference to address non-adherence in non inferiority trials####################
#############################################Analysis functions#################################################
################################################################################################################

#####################Non-adherence caused by non confounding process#############################
#################################################################################################

#BIAS
bias.nonconfounding <- function(n, p.experiment, p.stdcare, p.alt, nIterations, interval, nonadhere.pop, cross.over, true.effect, ymin, ymax){  
  
  #get estimates for each interval over nIterations
  estimate.df = sim.analysis(nonconfounding = 'nonconfounding', bias = T, nonadhere.pop = nonadhere.pop, interval = interval, cross.over = cross.over, 
                             confounder.intervention = confounder.intervention, confounder.outcome = confounder.outcome, NImargin = NImargin)
  
  #plot
  bias.plot = plot.eff.multi(estimate = estimate.df, cross.over = cross.over, ymin = ymin, ymax = ymax)
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.nonconfounding <-  function(n, p.experiment, p.alt, nIterations, interval, NImargin, nonadhere.pop, cross.over){  
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare = p.experiment - NImargin 
  
  #get estimates for each interval over nIterations
  df = sim.analysis(nonconfounding = 'nonconfounding', bias = F, nonadhere.pop = nonadhere.pop, interval = interval, cross.over = cross.over, 
                    confounder.intervention = confounder.intervention, confounder.outcome = confounder.outcome, NImargin = NImargin)
  
  plot = plot.t1(df = df, cross.over = cross.over)
  
  return(plot)
}

##################### Non-adherence caused by confounding process  ################################
###################################################################################################
#BIAS
bias.confounding <- function(n, p.experiment, p.stdcare, p.alt, confounder.intervention, confounder.outcome, interval, nIterations, nonadhere.pop, cross.over, true.effect, ymin, ymax){  
  
  #get estimates for each interval over nIterations
  estimate.df = sim.analysis(nonconfounding = 'confounding', bias = T, nonadhere.pop = nonadhere.pop, interval = interval, cross.over = cross.over, 
                             confounder.intervention = confounder.intervention, confounder.outcome = confounder.outcome, NImargin = NImargin)
  
  #plot
  bias.plot = plot.eff.multi(estimate = estimate.df, cross.over = cross.over, ymin = ymin, ymax = ymax)
  
  return(bias.plot)
} 

#TYPE 1 ERROR
type1.confounding <- function(n, p.experiment, p.alt, NImargin, confounder.intervention, confounder.outcome, interval, nIterations, cross.over, nonadhere.pop){  
  
  ##build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p.stdcare = p.experiment-NImargin 
  
  #get estimates for each interval over nIterations
  df = sim.analysis(nonconfounding = 'confounding', bias = F, nonadhere.pop = nonadhere.pop, interval = interval, cross.over = cross.over, 
                    confounder.intervention = confounder.intervention, confounder.outcome = confounder.outcome, NImargin = NImargin)
  
  plot = plot.t1(df = df, cross.over = cross.over)
  
  return(plot)
  
}

############################ Non-adherence caused by unknown confounding process  ##################
####################################################################################################
#########analysis of multiple unknown covariates with IPW 

bias.unknownconfounding.multi =  function(n, p.experiment, p.stdcare, confounder.intervention, confounder.outcome,interval,nIterations,nonadhere.pop, ymin, ymax){  
  
  #get estimates for each interval over nIterations
  estimate = sim.analysis(nonconfounding = 'unknownconfounding', bias = T, nonadhere.pop = nonadhere.pop, interval = interval, cross.over = cross.over, 
                          confounder.intervention = confounder.intervention, confounder.outcome = confounder.outcome, NImargin = NImargin)
  
  #plot
  iv = plot.eff(df = estimate, method = analysis.method[1], nIterations = nIterations, true.effect = true.effect,  ymin = ymin, ymax = ymax)
  itt = plot.eff(df = estimate, method = analysis.method[2], nIterations = nIterations, true.effect = true.effect, ymin = ymin, ymax = ymax)
  pp = plot.eff(df = estimate, method = analysis.method[4], nIterations = nIterations, true.effect = true.effect, ymin =ymin, ymax = ymax)
  
  x = c()
  
  for (h in 1:length(interval)) {
    x[[h]]=estimate[[h]][,c(3,5,6,7)]
  }
  
  mppdata=as.data.frame(matrix(c(unlist(x),
                                 rep(interval, each=length(x[[1]])), 
                                 rep(rep(1:4, each=nIterations),length(interval))), ncol = 3))
  colnames(mppdata)=c('y','x','type')
  mppdata$type=as.factor(mppdata$type)
  
  # The colour palette 
  cbPalette.multi  =  c("#e69f00", "#E6AE00", "#E6BD00", "#E6CD00")
  
  mpp= ggplot(mppdata)+ 
    geom_point(aes(x=x, y=y, color=type, alpha=0.0001), size=0.2) +
    geom_smooth(aes(x=x,y=y, colour=type), method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants in each arm")+
    ylab("Treatment effect")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette.multi)+
    theme(legend.title=element_blank(), legend.text=element_text(size=legendfontsize), legend.position="none")+
    scale_x_continuous(limits=c(start.interval, 1))+
    scale_y_continuous(limits=c(ymin, ymax))+
    geom_hline(yintercept=true.effect, linetype='dashed', color='red', size=0.5)
  
  bias.plot = ggarrange(itt, pp, mpp, iv,
                        ncol = 2, nrow = 2,
                        legend = "none")
  
  return(bias.plot)
} 


