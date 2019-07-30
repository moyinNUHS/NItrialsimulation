rm(list=ls())
setwd("/Users/moyin/Desktop/NItrialsimulation-master/") #set working directory 
source(file = 'causalinference 11 March 2019 (no AT and matching).R')

######parameters######
#assuming that the experimental and control treatments have the same treatment efficacy of 60%
#                  non-inferiority margin of 10% 
#                  type 1 error of 0.025
#                  90% power
#This required with 505 participants per arm 

n=505
p.experiment=0.4
p.stdcare=0.3
true.effect=p.experiment-p.stdcare
nIterations=1000
NImargin=0.1

width=10
height= 15

#Run simulations for case 1 (both groups non compliant)
b1both<-bias.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare,  nIterations=nIterations, interval = interval, noncomply = "both",true.effect=true.effect)
t1both<-type1.nonconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  nIterations=nIterations, interval = interval,noncomply = "both")
case1<-ggarrange(b1both, t1both,  ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case1", Sys.Date(), ".pdf"), width =width, height = height, units = c("in"))

#Run simulations for case 2: Non-compliance caused by confounding process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2both1<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both1<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
case2<-ggarrange(b2both1, t2both1, ncol = 2, labels = c('A','B'), 
                    common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case2", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## higher value of confounder makes intervention more likely, outcome more likely
b2both2<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both2<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
case3<-ggarrange(b2both2, t2both2, ncol = 2, labels = c('A','B'), 
                 common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case3", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#Run simulations for case 3: Non compliance with unknwon confounders 
bknown<-bias.unknownconfounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
tknown<-type1.unknownconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
case4<-ggarrange(bknown, tknown, ncol = 2, labels = c('A','B'), 
                 common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case4", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

bmulti<-bias.unknownconfounding.multi(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, noncomply = "both")
