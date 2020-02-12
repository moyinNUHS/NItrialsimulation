rm(list=ls())
source(file = 'setup_code100220.R')

######parameters######
#assuming that the experimental and control treatments have the same treatment efficacy of 60%
#                  non-inferiority margin of 10% 
#                  type 1 error of 0.025
#                  90% power
#This required with 505 participants per arm 

n = 505
p.experiment= 0.4
p.stdcare = 0.3
p.alt = 0.5
true.effect = p.experiment-p.stdcare
nIterations = 100
NImargin = 0.1

width = 15
height = 10

legendfontsize=12

#Run simulations for case 1 (both groups non adherent, non confounding) 
b1both.cross <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, nIterations = nIterations, interval = interval, nonadhere.pop = "both", true.effect = true.effect, ymin = -0.2, ymax = 0.3)
t1both.cross <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = NImargin,  p.alt = p.alt, cross.over = T, nIterations = nIterations, interval = interval, nonadhere.pop = "both")
case1.cross <- ggarrange(b1both.cross, t1both.cross,  ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case1.cross", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))

#Run simulations for case 2: Non-compliance caused by confounding process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2both1.cross <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "both", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.cross <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, interval = interval, nonadhere.pop = "both")
case2.cross <- ggarrange(b2both1.cross, t2both1.cross, ncol = 2, labels = c('A','B'), 
                    common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case2.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## Run simulations for case 3: higher value of confounder makes intervention more likely, outcome more likely
b2both2.cross <- bias.confounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "both", true.effect=true.effect,ymin = -0.2, ymax = 0.4)
t2both2.cross <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, interval = interval, nonadhere.pop = "both")
case3.cross <- ggarrange(b2both2.cross, t2both2.cross, ncol = 2, labels = c('A','B'), 
                 common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case3.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#Run simulations for case 4: Non compliance with unknown confounders 
bmulti <- bias.unknownconfounding.multi(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "both", ymin = -0.1, ymax = 0.3)
ggsave(filename = paste("unknownmulti", Sys.Date(), ".pdf"), width = width, height = height*2/3, units = "in")


##########################################################################
#########################Supplementary figures############################
###non confounding factors 
#both groups non adherent 
#non adherent participants receive inferior treatment 
b1both.inf.S2b <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, nIterations = nIterations, interval = interval, nonadhere.pop = "both", true.effect=true.effect, ymin = -0.2, ymax = 0.3)
t1both.inf.S2b <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = nImargin,  p.alt = p.alt, cross.over = F, nIterations=2000, interval = interval, nonadhere.pop = "both")
case1.inf.S2b <- ggarrange(b1both.inf.S2b, t1both.inf.S2b,  ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case1.inf.S2b", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))

#experimental group non adherent 
#cross over 
b1both.S2c <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, nIterations = nIterations, interval = interval, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.2, ymax = 0.3)
t1both.S2c <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = nImargin,  p.alt = p.alt, cross.over = T, nIterations = 3000, interval = interval, nonadhere.pop = "experimental")
S2c.cross <- ggarrange(b1both.S2c, t1both.S2c,  ncol = 2, labels = c('A','B'), 
                       common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2c.cross", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))
#non adherent participants receive inferior treatment 
b1both.S2d <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, nIterations = nIterations, interval = interval, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.2, ymax = 0.3)
t1both.S2d <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = nImargin,  p.alt = p.alt, cross.over = F, nIterations = 3000, interval = interval, nonadhere.pop = "experimental")
S2d.inf <- ggarrange(b1both.S2d, t1both.S2d,  ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2d.inf", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))

#control group non adherent 
#cross over 
b1both.S2e <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, nIterations = nIterations, interval = interval, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.2, ymax = 0.3)
t1both.S2e <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = nImargin,  p.alt = p.alt, cross.over = T, nIterations = 3000, interval = interval, nonadhere.pop = "stdcare")
S2e.cross <- ggarrange(b1both.S2e, t1both.S2e,  ncol = 2, labels = c('A','B'), 
                   common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2e.cross", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))
#non adherent participants receive inferior treatment 
b1both.S2f <- bias.nonconfounding(n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, nIterations = nIterations, interval = interval, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.2, ymax = 0.3)
t1both.S2f <- type1.nonconfounding(n = n, p.experiment = p.experiment, NImargin = nImargin,  p.alt = p.alt, cross.over = F, nIterations = 3000, interval = interval, nonadhere.pop = "stdcare")
S2f.inf <- ggarrange(b1both.S2f, t1both.S2f,  ncol = 2, labels = c('A','B'), 
                   common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2f.inf", Sys.Date(), ".pdf"), width = width, height = height, units = c("in"))

###confounding factors 
## higher value of confounder makes intervention less likely, outcome more likely
#both groups 
# non adherent participants receive inferior treatment 
b2both1.inf.2h <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "both", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.inf.2h <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin,  confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "both")
case2.inf.2h <- ggarrange(b2both1.inf.2h, t2both1.inf.2h, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case2.inf.2h", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#experimental group 
# cross over
b2both1.S2i <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2i <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "experimental")
S2i.cross <- ggarrange(b2both1.S2i, t2both1.S2i, ncol = 2, labels = c('A','B'), 
                       common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2i.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")
#non adherent participants receive inferior treatment 
b2both1.S2j <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval, confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2j <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin, confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "experimental")
S2j.inf <- ggarrange(b2both1.S2j, t2both1.S2j, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2j.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#control group 
# cross over
b2both1.S2k <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2k <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "stdcare")
S2k.cross <- ggarrange(b2both1.S2k, t2both1.S2k, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2k.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")
#non adherent participants receive inferior treatment 
b2both1.S2l <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2l <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin,  confounder.intervention = "Decrease likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "stdcare")
S2l.inf <- ggarrange(b2both1.S2l, t2both1.S2l, ncol = 2, labels = c('A','B'), 
                   common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2l.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

## higher value of confounder makes intervention and outcome more likely
#both groups
#non adherent participants receive inferior treatment 
b2both2.inf.2n <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "both", true.effect=true.effect,ymin = -0.2, ymax = 0.4)
t2both2.inf.2n <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "both")
case3.inf.2n <- ggarrange(b2both2.inf.2n, t2both2.inf.2n, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("case3.inf.2n", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#experimental group 
# cross over
b2both1.S2o <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2o  <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "experimental")
S2o.cross <- ggarrange(b2both1.S2o, t2both1.S2o, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2o.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")
#non adherent participants receive inferior treatment 
b2both1.S2f <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "experimental", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2f <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "experimental")
S2p.inf <- ggarrange(b2both1.S2p, t2both1.S2p, ncol = 2, labels = c('A','B'), 
                   common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2p.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")

#control group 
# cross over
b2both1.S2q <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = T, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2q <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = T, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "stdcare")
S2q.cross <- ggarrange(b2both1.S2q, t2both1.S2q, ncol = 2, labels = c('A','B'), 
                     common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2q.cross", Sys.Date(), ".pdf"), width = width, height = height, units = "in")
#non adherent participants receive inferior treatment 
b2both1.S2r <- bias.confounding (n = n, p.experiment = p.experiment, p.stdcare = p.stdcare, p.alt = p.alt, cross.over = F, interval = interval,confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations, nonadhere.pop = "stdcare", true.effect=true.effect, ymin = -0.15, ymax = 0.3)
t2both1.S2r <- type1.confounding(n = n, p.experiment = p.experiment, p.alt = p.alt, cross.over = F, NImargin = nImargin,  confounder.intervention = "Increase likelihood", confounder.outcome = "Increase likelihood", nIterations = nIterations*3, interval = interval, nonadhere.pop = "stdcare")
S2r.inf <- ggarrange(b2both1.S2r, t2both1.S2r, ncol = 2, labels = c('A','B'), 
                   common.legend = TRUE, legend = "bottom" )
ggsave(filename = paste("S2r.inf", Sys.Date(), ".pdf"), width = width, height = height, units = "in")
