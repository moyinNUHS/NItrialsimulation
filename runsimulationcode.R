rm(list=ls())
setwd("/Users/moyin/Desktop/angelsfly/NItrialsimulation/codes/") #set working directory 
source(file = 'causalinference 11 March 2019 (no AT and matching).R')

######parameters######
n=220
p.experiment=0.35
p.stdcare=0.4
true.effect=p.experiment-p.stdcare
nIterations=10000
NImargin=0.1

width=10
height= 15

#Run simulations for case 1.1 (both groups non compliant)
b1both<-bias.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare,  nIterations=nIterations, interval = interval, noncomply = "both",true.effect=true.effect)
t1both<-type1.nonconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  nIterations=nIterations, interval = interval,noncomply = "both")
p1both<-power.nonconfounding(n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, interval = interval,noncomply = "both")
case1both<-ggarrange(b1both, 
                     ggarrange (t1both, p1both, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                     nrow=2, 
                     labels ='A',
                     common.legend = TRUE, legend = "bottom" )
annotate_figure(case1both, top = text_grob("Case 1.1: Non-compliance caused by non confounding processes affecting both groups", face = "bold", size = 12))
ggsave(filename = paste("case1.1", Sys.Date(), ".pdf"), width =width, height = height, units = c("in"))

#Run simulations for case 1.2 (one groups non compliant)
b1one<-bias.nonconfounding  (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare,  nIterations=nIterations, interval = interval, noncomply = "experimental", true.effect=true.effect)
t1one<-type1.nonconfounding (n=n, p.experiment=p.experiment, NImargin=NImargin,  nIterations=nIterations, interval = interval,noncomply = "experimental")
p1one<-power.nonconfounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, interval = interval,noncomply = "experimental")
case1.2<-ggarrange(b1one, 
                   ggarrange (t1one, p1one, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                   nrow=2, 
                   labels ='A',
                   common.legend = TRUE, legend = "bottom" )
annotate_figure(case1.2, top = text_grob("Case 1.2: Non-compliance caused by non confounding processes affecting intervention or control group", face = "bold", size = 12))
ggsave(filename = paste("case1.2", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#Run simulations for case 2.1: Non-compliance caused by confounding process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2both1<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both1<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
p2both1<-power.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", interval = interval, noncomply = "both")
case2.11<-ggarrange(b2both1, 
                    ggarrange (t2both1, p2both1, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.11, top = text_grob("Case 2.11: Non-compliance caused by confounding process (affect both groups)-higher value of confounder makes intervention less likely, outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.11", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

###more sample size 
t2both1more<-type1.confounding(n=n+100, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
p2both1more<-power.confounding (n=n+100, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Decrease likelihood", confounder.outcome="Increase likelihood", interval = interval, noncomply = "both")
case2.11more<-ggarrange(ggarrange(t2both1, p2both1,ncol = 2, labels = c(paste(LETTERS[1])), legend=NULL),
                        ggarrange(t2both1more, p2both1more,ncol = 2, labels = c(paste(LETTERS[2])), legend=NULL), 
                    nrow=2,
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.11more, top = text_grob("Case 2.11A: More sample size (+50) from case 2.11", face = "bold", size = 12))
ggsave(filename = paste("case2.11more", Sys.Date(), ".pdf"), width = 7.5, height = 7.5, units = "in")


## higher value of confounder makes intervention more likely, outcome more likely
b2both2<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both2<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
p2both2<-power.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Increase likelihood", confounder.outcome="Increase likelihood", interval = interval, noncomply = "both")
case2.12<-ggarrange(b2both2, 
                    ggarrange (t2both2, p2both2, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.12, top = text_grob("Case 2.12: Non-compliance caused by confounding process (affect both groups)-higher value of confounder makes intervention more likely, outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.12", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## higher value of confounder makes intervention less likely, outcome less likely
b2both3<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Decrease likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both3<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Decrease likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
p2both3<-power.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Decrease likelihood", confounder.outcome="Decrease likelihood", interval = interval, noncomply = "both")
case2.13<-ggarrange(b2both3, 
                    ggarrange (t2both3, p2both3, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.13, top = text_grob("Case 2.13: Non-compliance caused by confounding process (affect both groups)-higher value of confounder makes intervention less likely, outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.13", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## higher value of confounder makes intervention more likely, outcome less likely
b2both4<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
t2both4<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
p2both4<-power.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", interval = interval, noncomply = "both")
case2.14<-ggarrange(b2both4, 
                    ggarrange (t2both4, p2both4, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.14, top = text_grob("Case 2.14: Non-compliance caused by confounding process (affect both groups)-higher value of confounder makes intervention more likely, outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case2.14", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))

#Run simulations for case 2.2: Non-compliance caused by confounding process (affect intervention groups)
b2one<-bias.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, noncomply = "experimental", true.effect=true.effect)
t2one<-type1.confounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, noncomply = "experimental")
p2one<-power.confounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", interval = interval, noncomply = "experimental")
case2.21<-ggarrange(b2one, 
                    ggarrange (t2one, p2one, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case2.21, top = text_grob("Case 2.21: Non-compliance caused by confounding process (affect intervention group)-higher value of confounder makes outcome more likely", face = "bold", size = 12))
ggsave(filename = paste("case2.21", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#Run simulations for case 3: Non compliance with unknwon confounders 
bknown<-bias.unknownconfounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, interval = interval,confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood",nIterations=nIterations, noncomply = "both", true.effect=true.effect)
tknown<-type1.unknownconfounding(n=n, p.experiment=p.experiment, NImargin=NImargin,  confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", nIterations=nIterations, interval = interval, noncomply = "both")
pknown<-power.unknownconfounding (n=n, p.experiment=p.experiment, p.stdcare=p.stdcare, NImargin=NImargin, nIterations=nIterations, confounder.intervention="Increase likelihood", confounder.outcome="Decrease likelihood", interval = interval, noncomply = "both")
case3<-ggarrange(bknown, 
                    ggarrange (tknown, pknown, ncol = 2, labels = c(paste(LETTERS[2:3], sep="")), legend=NULL), 
                    nrow=2, 
                    labels ='A',
                    common.legend = TRUE, legend = "bottom" )
annotate_figure(case3, top = text_grob("Case 2.21: Non-compliance caused by unknown confounding process (affect both groups)-higher value of confounder makes intervention more likely, outcome less likely", face = "bold", size = 12))
ggsave(filename = paste("case3", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

