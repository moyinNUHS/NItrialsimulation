##need to place this file in main NItrialsimulation folder to source for codes
rm(list=ls())
source(file = 'setup_code100220.R')

sim.conf.eff.outcome = function(n = 505, p.experiment= 0.4, p.stdcare = 0.3, p.alt = 0.5, nIterations = nIterations, 
                                adhere.experiment = interval, adhere.stdcare = interval, 
                                NImargin = 0.1,
                                confounder.outcome = 'Increase likelihood',
                                confounder.intervention = 'Decrease likelihood',
                                cross.over = T, confounder.eff) {
  
  #make up vectors for simulations 
  .estimate = c() #output from each simulation
  estimate = c()  #for saving output from each interval
  
  for(i in 1:length(interval)) { print(paste("interval value", i,"out of",length(interval)))
    for(l in 1:nIterations) { 
      tryCatch({
        simdata = simdata.confounding(n = n, p.experiment = p.experiment, p.stdcare= p.stdcare, p.alt = p.alt, cross.over = cross.over, nonadhere.pop = nonadhere.pop, adhere.experiment = adhere.experiment, adhere.stdcare = adhere.stdcare, 
                                      confounder.outcome = confounder.outcome, confounder.intervention = confounder.intervention, 
                                      confounder.eff = confounder.eff, i=i)
        .estimate[[l]] = analysis.estimate(simdata = simdata) #gives a vector of estimates in the order of iv, itt, ipw, pp for each iteration
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")}) #receive error message if there is an error 
    }
    #save results of every interval - each column refers to each analysis method 
    estimate[[i]] = matrix(unlist(.estimate), ncol = length(analysis.method), byrow = TRUE)
  }
  
  x = y = c()
  for (i in 1:length(interval)) {
    x[[i]] = estimate[[i]][,4]
    y[[i]] = rep(interval[i],dim(estimate[[i]])[1])
  }
  
  plotdata = as.data.frame(matrix(c(unlist(x),unlist(y)), ncol = 2))
  plotdata$coeff.eff = as.factor(confounder.eff)
  
  return(plotdata)
}
t1.conf.eff.outcome = function(estimates, NImargin, nIterations, coeff.eff){
  t1 = c()
  #get type 1 error for each interval 
  for (i in 1:length(interval)) {
    df = estimates[estimates$V2 == interval[i],]
    #get sds for each interval 
    sd = sd(df$V1)
    #get upper bounds of 95%CI for each iteration of every interval
    ub = df$V1 + qnorm(0.975) * sd
    
    t1[[i]] = mean(ub < NImargin)
  }
  
  t1.df = as.data.frame(matrix(c(interval, unlist(t1)), ncol = 2))
  t1.df$coeff.eff = as.factor(coeff.eff)
  return(t1.df)
}


nIterations = 500
coeff1 = sim.conf.eff.outcome(n = 505, p.experiment= 0.4, p.stdcare = 0.3, p.alt = 0.5, nIterations = nIterations, 
                     adhere.experiment = interval, adhere.stdcare = interval, 
                     NImargin = 0.1,
                     confounder.outcome = 'Increase likelihood',
                     confounder.intervention = 'Decrease likelihood',
                     cross.over = T, confounder.eff = 1)
coeff5 = sim.conf.eff.outcome(n = 505, p.experiment= 0.4, p.stdcare = 0.3, p.alt = 0.5, nIterations = nIterations, 
                              adhere.experiment = interval, adhere.stdcare = interval, 
                              NImargin = 0.1,
                              confounder.outcome = 'Increase likelihood',
                              confounder.intervention = 'Decrease likelihood',
                              cross.over = T, confounder.eff = 3)
coeff9 = sim.conf.eff.outcome(n = 505, p.experiment= 0.4, p.stdcare = 0.3, p.alt = 0.5, nIterations = nIterations, 
                              adhere.experiment = interval, adhere.stdcare = interval, 
                              NImargin = 0.1,
                              confounder.outcome = 'Increase likelihood',
                              confounder.intervention = 'Decrease likelihood',
                              cross.over = T, confounder.eff = 9)

plotdata.e = rbind(coeff1, coeff5, coeff9)

width = 10
height = 7.5
legendfontsize=12
true.effect = 0.1
colours = c(cbPalette[4], "#6096BA", "#A3CEF1")
est = ggplot(plotdata.e, aes(x=V2, y=V1, color = coeff.eff))+ 
  geom_point(aes( alpha=0.0001), size=0.2) +
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  guides(alpha=FALSE)+
  xlab("Proportion of adherent participants in each arm")+
  ylab("Effect estimate")+
  theme_minimal()+
  scale_colour_manual(values=colours)+
  theme(legend.title=element_blank(), legend.position="none", legend.text=element_text(size=legendfontsize))+
  scale_x_continuous(limits=c(0.58, 1), breaks = seq(0.6, 1, by=0.1))+
  scale_y_continuous(limits=c(-0.1, 0.2))+
  geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = 0.2),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
  geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = -0.1),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
  annotate(geom="text", x=0.58, y=0.11, angle = 90, label="Favour control", color="grey40", hjust=0) +
  annotate(geom="text", x=0.58, y=0.09, angle = 90, label="Favour experiment", color="grey40", hjust=1) +
  annotate(geom="text", x=0.65, y=0.09,  label='(a)', color="grey40", hjust=0) +
  annotate(geom="text", x=0.7, y=0.06,label= '(b)', color="grey40", hjust=1) +
  annotate(geom="text", x=0.75, y=0.03,label= '(c)', color="grey40", hjust=1) +
  geom_hline(yintercept=true.effect, linetype='dashed', color='red', size=0.5)

t1.coef1 = t1.conf.eff.outcome(estimates = coeff1, NImargin=0.1, nIterations=nIterations, coeff.eff = 1)
t1.coef5 = t1.conf.eff.outcome(estimates = coeff5, NImargin=0.1, nIterations=nIterations, coeff.eff = 5)
t1.coef9 = t1.conf.eff.outcome(estimates = coeff9, NImargin=0.1, nIterations=nIterations, coeff.eff = 9)

plotdata.t = rbind(t1.coef1, t1.coef5, t1.coef9)
t1 = ggplot(plotdata.t, aes(x=V1, y=V2, color = coeff.eff))+ 
  geom_point(aes(alpha=0.0001), size=0.2) +
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  guides(alpha=FALSE)+
  xlab("Proportion of adherent participants in each arm")+
  ylab("Type 1 error")+
  theme_minimal()+
  scale_colour_manual(values=colours)+
  theme(legend.title=element_blank(), legend.position="none", legend.text=element_text(size=legendfontsize))+
  scale_x_continuous(limits=c(0.6, 1), breaks = seq(0.6, 1, by=0.1))+
  scale_y_continuous(limits=c(0, 0.8))+
  annotate(geom="text", x=0.65, y=0.09,  label='(a)', color="grey40", hjust=0) +
  annotate(geom="text", x=0.7, y=0.25,label= '(b)', color="grey40", hjust=1) +
  annotate(geom="text", x=0.78, y=0.4,label= '(c)', color="grey40", hjust=1) +
  geom_hline(yintercept=0.025, linetype='dashed', color='red')

p = ggarrange(est, t1, labels = c('A', 'B'))
p
ggsave(filename = paste("conf.eff.outcome", Sys.Date(), ".jpeg"), width = width, height = height, units = "in")

