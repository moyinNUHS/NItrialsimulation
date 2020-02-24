analysis.power <- function(estimate.iter, z, NImargin, nIterations){
  
  estimates.df = matrix(unlist(estimate.iter), ncol = 4, byrow = T)
  sds = unlist(apply(estimates.df, 2, sd))
  upper.bounds.df = estimates.df + rep(z * sds, each = nIterations)
  
  p.itt = sum(upper.bounds.df[, 1] < NImargin)
  p.pp = sum(upper.bounds.df[, 2] < NImargin)
  p.mpp = sum(upper.bounds.df[, 3] < NImargin)
  p.iv = sum(upper.bounds.df[, 4] < NImargin)
  
  return(c(p.itt, p.pp, p.mpp, p.iv)/nIterations)
}