analysis.estimate <- function(simdata) {
  
  n = nrow(simdata)/2
  
  ## intention to treat
  pz1.value = sum(simdata[simdata[, 2] == 1, ][, 5])/n
  pz0.value = sum(simdata[simdata[, 2] == 0, ][, 5])/n
  eff.itt = pz1.value - pz0.value
  
  ## per protocol
  pp = simdata[simdata[, 2] == simdata[, 4], ] # per protocol population
  p.experiment.value = sum(pp[pp[, 2] == 1,][, 5])/length(pp[pp[, 2] == 1,][, 5])
  p.stdcare.value = sum(pp[pp[, 2] == 0, ][, 5])/length(pp[pp[, 2] == 0, ][, 5])
  eff.pp = p.experiment.value - p.stdcare.value
  
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
  mean = sum(pp$intervention) / nrow(pp)
  weight = pp$intervention * mean / score + (1 - pp$intervention) *
    (1 - mean) / (1 - score) #create stabilized weights, using a null model with intervention as the dependent variable
  outcomemodel = speedglm(
    outcome ~ intervention,
    family = binomial(link = "identity"),
    weights = weight,
    data = pp
  ) #identity link for risk difference
  eff.mpp = coef(outcomemodel)[2]
  
  # iv with 2 stage regression
  asmm = gmm(simdata[, 5] ~ simdata[, 4], x = simdata[, 2], vcov = "iid")
  eff.iv = summary(asmm)$coefficients[2, 1]
  
  return(c(eff.itt, eff.pp, eff.mpp, eff.iv))
}