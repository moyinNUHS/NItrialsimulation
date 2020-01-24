rm(list=ls())
setwd("/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/") #set working directory 
source(file = 'simulationcode_200120.R')

######parameters######
#assuming that the experimental and control treatments have the same treatment efficacy of 60%
#                  non-inferiority margin of 10% 
#                  type 1 error of 0.025
#                  90% power
#This required with 505 participants per arm 

n=505
p.experiment=0.4
p.stdcare=0.3
p.alt = 0.45
true.effect=p.experiment-p.stdcare
nIterations=1000
NImargin=0.1

width=15
height= 10

legendfontsize=12

#Run simulations for case 1 (both groups non adherent) 
b1both<-bias.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = T, nIterations=nIterations, interval = interval, nonadhere.pop = "both", true.effect=true.effect, ymin=-0.2, ymax=0.3)
t1both<-type1.nonconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  p.alt = p.alt, cross.over = T, nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case1.cross<-ggarrange(b1both, t1both,  ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case1.cross", Sys.Date(), ".pdf"), width =width, height = height, units = c("in"))

b1both<-bias.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = F, nIterations=nIterations, interval = interval, nonadhere.pop = "both", true.effect=true.effect, ymin=-0.2, ymax=0.3)
t1both<-type1.nonconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  p.alt = p.alt, cross.over = F, nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case1.inf<-ggarrange(b1both, t1both,  ncol = 2, labels = c('A','B'), 
                 common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case1.inf", Sys.Date(), ".pdf"), width =width, height = height, units = c("in"))

#Run simulations for case 2: Non-compliance caused by confounding process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2both1<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, nonadhere.pop = "both", true.effect=true.effect, ymin=-0.15, ymax=0.3)
t2both1<-type1.confounding(n=n, p.experiment=p.experiment, p.alt = p.alt, cross.over = T, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case2.cross<-ggarrange(b2both1, t2both1, ncol = 2, labels = c('A','B'), 
                    common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case2.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

b2both1<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, nonadhere.pop = "both", true.effect=true.effect, ymin=-0.15, ymax=0.3)
t2both1<-type1.confounding(n=n, p.experiment=p.experiment, p.alt = p.alt, cross.over = F, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case2.inf<-ggarrange(b2both1, t2both1, ncol = 2, labels = c('A','B'), 
                       common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case2.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## Run simulations for case 3: higher value of confounder makes intervention more likely, outcome more likely
b2both2<-bias.confounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, nonadhere.pop = "both", true.effect=true.effect,ymin=-0.2, ymax=0.4)
t2both2<-type1.confounding(n=n, p.experiment=p.experiment, p.alt = p.alt, cross.over = T, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case3.cross<-ggarrange(b2both2, t2both2, ncol = 2, labels = c('A','B'), 
                 common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case3.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

b2both2<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, nonadhere.pop = "both", true.effect=true.effect,ymin=-0.2, ymax=0.4)
t2both2<-type1.confounding(n=n, p.experiment=p.experiment, p.alt = p.alt, cross.over = F, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, nonadhere.pop = "both")
case3.inf<-ggarrange(b2both2, t2both2, ncol = 2, labels = c('A','B'), 
                       common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case3.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#Run simulations for case 4: Non compliance with unknown confounders 
#bknown<-bias.unknownconfounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, nonadhere.pop = "both", true.effect=true.effect)
#tknown<-type1.unknownconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, nonadhere.pop = "both")
#case4<-ggarrange(bknown, tknown, ncol = 2, labels = c('A','B'), 
#                 common.legend = TRUE, legend = "bottom" )
#ggsave(filename = paste("case4", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

bmulti<-bias.unknownconfounding.multi(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, nonadhere.pop = "both", ymin=-0.1, ymax=0.3)
ggsave(filename = paste("unknownmulti", Sys.Date(), ".pdf"), width = width, height = height*2/3, units = "in")

