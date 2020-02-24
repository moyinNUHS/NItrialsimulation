rm(list = ls())
library('gmm'); library('speedglm')

####Deploy app
# setwd('/Users/moyin/Documents/nBox/git_projects/NItrialsimulation/shiny/samplesize_nonadherence/')
# require('rsconnect')
# rsconnect::setAccountInfo(name='moru',
#                           token='30819BAEDD492333CE7CD293F3B08D42',
#                           secret='j24TwHndsiybKqQZERFOt0s3L3ro8IG47/OgpLk/')
# deployApp(account="moru",appName="samplesize_nonadherence")

####parameters####
p.experiment = p.stdcare = 0.4
NImargin = 0.1
n = 505
adhere.stdcare = 1
adhere.experiment = 1
confounder.intervention = 'Decrease likelihood'
confounder.outcome = 'Increase likelihood'

###functions###
analysis.ub <- function(simdata) {
  ## intention to treat
  pz1.value = sum(simdata[which(simdata[, 2] == 1),][, 5])/n
  pz0.value = sum(simdata[which(simdata[, 2] == 0),][, 5])/n
  eff.itt = pz1.value - pz0.value
  #sd.itt = sqrt(pz1.value * (1 - pz1.value) / n + pz0.value * (1 - pz0.value) / n)
  
  ## per protocol
  pp = simdata[which(simdata[, 2] == simdata[, 4]),] # per protocol population
  p.experiment.vector = pp[which(pp[, 2] == 1), ][, 5]
  p.experiment.value = sum(p.experiment.vector)/length(p.experiment.vector)
  p.stdcare.vector = pp[which(pp[, 2] == 0),][, 5]
  p.stdcare.value = sum(p.stdcare.vector)/length(p.stdcare.vector)
  eff.pp = p.experiment.value - p.stdcare.value
  #sd.pp = sqrt(
  #  p.experiment.value * (1 - p.experiment.value) / 
  #    length(p.experiment.vector) + p.stdcare.value * (1 - p.stdcare.value) / length(p.stdcare.vector)
  #)
  
  ## inverse probability weights on per protocol patients
  pp = as.data.frame(pp)
  colnames(pp) = c('id',
                   'randomisation',
                   'confounder',
                   'intervention',
                   'outcome')
  ipwmodel = glm(intervention ~ confounder,
                 family = binomial(link = "logit"),
                 data = pp) #calculate denominators used in inverse probability weights
  score = predict(ipwmodel, type = "response")
  mean = (sum(pp$intervention)/length(pp$intervention))
  weight = pp$intervention * mean / score + (1 - pp$intervention) *
    (1 - mean) / (1 - score) #create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel = speedglm(
    outcome ~ intervention,
    family = binomial(link = "identity"),
    weights = weight,
    data = pp
  ) #identity link for risk difference
  eff.mpp = coef(outcomemodel)[2]
  #sd.mpp = sqrt(diag(vcovHC(outcomemodel)))[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[, 5] ~ simdata[, 4], x = simdata[, 2], vcov = "iid")
  eff.iv = summary(asmm)$coefficients[2, 1]
  #sd.iv = summary(asmm)$coefficients [2, 2]
  
  return(list (
    c(eff.itt, eff.pp, eff.mpp, eff.iv)
    #, c(sd.itt, sd.pp, eff.mpp, eff.iv)
  ))
}

getoutcome <- function(outcome0, outcome1, intervention) {
  outcome = c()
  outcome.matrix = matrix(c(outcome0, outcome1), ncol = 2)
  for (i in 1:length(intervention)) {
    outcome[i] = outcome.matrix[i, (intervention[i] + 1)]
  }
  return(unlist(outcome))
}

########Non-confounding##########
nc.power <- function (n,
                      p.experiment,
                      p.stdcare,
                      adhere.experiment,
                      adhere.stdcare,
                      NImargin) {
  nIterations = 1000
  
  id = seq(1, (2 * n), by = 1) #create participant id
  randomisation = sample(rep(0:1, n)) #randomisation
  confounder = rbeta(n = 2 * n, shape1 = 2, shape2 = 2) #confounder beta distribution ranging 0-1
  
  ###ALTERNATIVE HYPOTHESIS for power
  means.alt.nc = sds.alt.nc = c()
  
  for (i in 1:nIterations) {
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    outcome0 = rbinom(2 * n, prob = p.stdcare, size = 1)
    outcome1 = rbinom(2 * n, prob = p.experiment, size = 1)
    
    #INTERVENTION dependent on adherence
    intervention = rep(NA, 2 * n)
    intervention[sample(which(randomisation == 1), size = adhere.experiment * n)] = 1
    intervention[which(randomisation == 1) %in% which(is.na(intervention))] = 0
    intervention[sample(which(randomisation == 0), size = adhere.stdcare * n)] = 0
    intervention[which(randomisation == 0) %in% which(is.na(intervention))] = 1
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0 = outcome0,
                         outcome1 = outcome1,
                         intervention = intervention)
    simdata = matrix(
      data = c(id, randomisation, confounder, intervention, outcome),
      nrow = (2 * n)
    )
    
    ub = analysis.ub(simdata = simdata)
    means.alt.nc[[i]] = ub[[1]]
  }
  
  means.df.alt.nc = matrix(unlist(means.alt.nc), ncol = 4, byrow = T)
  sds.alt.nc = unlist(apply(means.df.alt.nc, 2, sd))
  ub.df.alt.nc = means.df.alt.nc + rep(qnorm(0.975) * sds.alt.nc, each = nIterations)
  
  p.itt.nc = sum(ub.df.alt.nc[, 1] < NImargin)
  p.pp.nc = sum(ub.df.alt.nc[, 2] < NImargin)
  p.mpp.nc = sum(ub.df.alt.nc[, 3] < NImargin)
  p.iv.nc = sum(ub.df.alt.nc[, 4] < NImargin)
  
  return (list(
    c(p.itt.nc, p.pp.nc, p.mpp.nc, p.iv.nc)/nrow(ub.df.alt.nc),
    ub.df.alt.nc
  ))
}

# p.505.nc = nc.power (n = 505, p.experiment, p.stdcare, adhere.experiment, adhere.stdcare, NImargin)

###########Confounding############
c.power <- function (n = 505,
                     p.experiment,
                     p.stdcare,
                     adhere.experiment = 1,
                     adhere.stdcare = 1,
                     confounder.intervention = 'Increase likelihood',
                     confounder.outcome = 'Increase likelihood',
                     NImargin = 0.1, 
                     shape2.min, 
                     shape2.max) {
  nIterations = 1000
  
  id = seq(1, (2 * n), by = 1) #create participant id
  randomisation = sample(rep(0:1, n)) #randomisation
  
  shape2.outcome = runif(1, min = shape2.min, max = shape2.max)
  shape2.intervention = runif(1, min = 2, max = 10)
  
  ###ALTERNATIVE HYPOTHESIS for power
  means.alt.c = sds.alt.c = eff.conf.out = c()
  
  for (i in 1:nIterations) {
    confounder = rbeta(n = (2 * n),
                       shape1 = 2,
                       shape2 = 2) #confounder beta distribution ranging 0-1
    
    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
    if (confounder.outcome == "Increase likelihood") {
      #probability of outcome is drawn from beta distribution shape1<1 and shape2<1 (U shaped) such that confounder correlates with outcome
      shape1 =  shape2.outcome * p.stdcare / (1 - p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind.gen = rbeta(
        n = (2 * n),
        shape1 = shape1,
        shape2 = shape2.outcome
      )
      p.stdcare.ind = p.stdcare.ind.gen[order(p.stdcare.ind.gen)]#individual probability with mean of p.stdcare, in increasing order
      outcome0 =  rbinom(2 * n, 1, prob = p.stdcare.ind) #increasing confounder value will have increasing probability for outcome
      
      shape1 =  shape2.outcome * p.experiment / (1 - p.experiment)
      p.experiment.ind.gen = rbeta(
        n = (2 * n),
        shape1 = shape1,
        shape2 = shape2.outcome
      )
      p.experiment.ind = p.experiment.ind.gen[order(p.experiment.ind.gen)] #individual probability with mean of p.experiment, in increasing order
      outcome1 =  rbinom(2 * n, 1, prob = p.experiment.ind) #increasing confounder value will have increasing probability for outcome
    } else {
      shape1 =  shape2.outcome * p.stdcare / (1 - p.stdcare) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      p.stdcare.ind.gen = rbeta(
        n = (2 * n),
        shape1 = shape1,
        shape2 = shape2.outcome
      )
      p.stdcare.ind = p.stdcare.ind.gen[order(p.stdcare.ind.gen, decreasing = TRUE)] #individual probability with mean of p.stdcare, in decreasing order
      outcome0 =  rbinom(2 * n, 1, prob = p.stdcare.ind) #increasing confounder value will have decreasing probability for outcome
      
      shape1 =  shape2.outcome * p.experiment / (1 - p.experiment)
      p.experiment.ind.gen = rbeta(
        n = (2 * n),
        shape1 = shape1,
        shape2 = shape2.outcome
      )
      p.experiment.ind =  p.experiment.ind.gen[order(p.experiment.ind.gen, decreasing = TRUE)] #individual probability with mean of p.experiment, in decreasing order
      outcome1 =  rbinom(2 * n, 1, prob = p.experiment.ind) #increasing confounder value will have decreasing probability for outcome
    }
    
    d = matrix(data = c(id, randomisation, confounder),
               ncol = 3)
    d.ordered = matrix(cbind(d[order(d[, 3]), ], outcome0, outcome1), ncol = 5) #order confounder in ascending order
    d.grouped = rbind(d.ordered[which(d.ordered[, 2] == 1),], d.ordered[which(d.ordered[, 2] == 0),]) #group according to randomisation
    
    #INTERVENTION dependent on randomisation and confounders
    if (confounder.intervention == "Increase likelihood") {
      shape1 =  shape2.intervention * adhere.experiment / (1 - adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind.gen = rbeta(
        n = n,
        shape1 = shape1,
        shape2 = shape2.intervention
      )
      adhere.experiment.ind = adhere.experiment.ind.gen[order(adhere.experiment.ind.gen)]#individual probability with mean of p.experiment, in increasing order
      int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for intervention
      
      shape1 =  shape2.intervention * (1 - adhere.stdcare) / (1 - (1 - adhere.stdcare))
      adhere.stdcare.ind.gen = rbeta(
        n = n,
        shape1 = shape1,
        shape2 = shape2.intervention
      )
      adhere.stdcare.ind = adhere.stdcare.ind.gen[order(adhere.stdcare.ind.gen)]#individual probability with mean of p.stdcare, in increasing order
      cont.intervention =  rbinom(n, 1, prob = adhere.stdcare.ind)
    } else {
      shape1 =  shape2.intervention * adhere.experiment / (1 - adhere.experiment) #mean of beta distribution is a/(a+b), a is shape1, b is shape2
      adhere.experiment.ind.gen = rbeta(
        n = n,
        shape1 = shape1,
        shape2 = shape2.intervention
      )
      adhere.experiment.ind = adhere.experiment.ind.gen[order(adhere.experiment.ind.gen, decreasing = TRUE)] #individual probability with mean of p.experiment, in decreasing order
      int.intervention =  rbinom(n, 1, prob = adhere.experiment.ind) #increasing confounder value will have decreasing probability for outcome
      
      shape1 =  shape2.intervention * adhere.stdcare / (1 - adhere.stdcare)
      adhere.stdcare.ind.gen = rbeta(
        n = n,
        shape1 = shape1,
        shape2 = shape2.intervention
      ) 
      adhere.stdcare.ind =adhere.stdcare.ind.gen[order(adhere.stdcare.ind.gen)]#individual probability with mean of p.stdcare, in decreasing order
      cont.intervention =  rbinom(n, 1, prob = 1 - adhere.stdcare.ind)
    }
    
    intervention = c(int.intervention, cont.intervention)
    
    #ACTUAL OUTCOMES depend on intervention
    outcome = getoutcome(outcome0 = d.grouped[, 4], outcome1 = d.grouped[, 5], intervention)
    simdata = matrix(cbind(d.grouped[,-5:-4], intervention, outcome), ncol = 5)
    
    # report effect of confounder on outcome 
    df.ordered = as.data.frame(d.ordered)
    model = speedglm(V4 ~ V3, family = binomial(link = "logit"), data = df.ordered)
    
    eff.conf.out[[i]] = coefficients(model)[2]
    
    ub = analysis.ub(simdata = simdata)
    means.alt.c[[i]] = ub[[1]]
  }
  
  means.df.alt.c = matrix(unlist(means.alt.c), ncol = 4, byrow = T)
  sds.alt.c = unlist(apply(means.df.alt.c, 2, sd))
  ub.df.alt.c = means.df.alt.c + rep(qnorm(0.975) * sds.alt.c, each = nIterations)
  
  p.itt.alt.c = sum(ub.df.alt.c[, 1] < NImargin)
  p.pp.alt.c = sum(ub.df.alt.c[, 2] < NImargin)
  p.mpp.alt.c = sum(ub.df.alt.c[, 3] < NImargin)
  p.iv.alt.c = sum(ub.df.alt.c[, 4] < NImargin)
  
  return (list(
    c(p.itt.alt.c, p.pp.alt.c, p.mpp.alt.c, p.iv.alt.c)/nIterations,
    ub.df.alt.c, 
    eff.conf.out
  ))
}

#p.505.c = c.power (n = 505, p.experiment, p.stdcare, adhere.experiment, adhere.stdcare, 
                 #confounder.intervention, confounder.outcome, NImargin, shape2.min = 1, shape2.max = 2)


# par(mfrow = c(1,1))
# ub.df.alt.nc = p.505.nc [[2]]
# ub.df.alt.c = p.505.c[[2]]
# ##compare alt distributions
# plot(density(ub.df.alt.nc[, 1]), xlim = c(-0.1, 0.15), main = paste('Effect of confounder on intervention =', 
#                                                                     round(p.505.c[[3]], 3)))
# lines(density(ub.df.alt.c[, 1]), col = 'red')
# abline(v = mean(ub.df.alt.nc[, 1]),
#        lwd = 1,
#        lty = 2)
# abline(
#   v = mean(ub.df.alt.c[, 1]),
#   col = 'red',
#   lwd = 1,
#   lty = 2
# )
# lines(
#   x = c(mean(ub.df.alt.nc[, 1]), mean(ub.df.alt.nc[, 1]) + sd(ub.df.alt.nc[, 1])),
#   col = 'grey',
#   y = c(3, 3),
#   lwd = 1
# )
# lines(
#   x = c(mean(ub.df.alt.c[, 1]), mean(ub.df.alt.c[, 1]) + sd(ub.df.alt.c[, 1])),
#   y = c(3.1, 3.1),
#   col = 'pink',
#   lwd = 1
# )
# text(round(sd(ub.df.alt.c[, 1]), digits = 4),
#      x = mean(ub.df.alt.c[, 1]) + 0.02,
#      y = 3.4)
# text(round(sd(ub.df.alt.nc[, 1]), digits = 4),
#      x = mean(ub.df.alt.nc[, 1]) + 0.02,
#      y = 2.6)
# ##compare z values
# abline(v = 0.1, col = 'blue', lwd = 0.5)
# ##compare power
# text(paste('power.nc =', p.505.nc[[1]][1]), 
#      x = -0.075, 
#      y = 10)
# text(paste('power.c =', p.505.c[[1]][1]), 
#      x = -0.075, 
#      y = 9, 
#      col = 'red')
# 
# #distribution of coefficient of outcome0 ~ confounder 
# plot(density(p.505.c[[3]]))
# mean(p.505.c[[3]])
# median((p.505.c[[3]]))
