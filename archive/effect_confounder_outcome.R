################ Determination of effect of confounders on outcome #############
################################################################################
rm(list = ls())
library('parallel')
source('archive/powercalculator.R')

#read previous csv
prev.df = read.csv('eff.conf.outcome.csv')[,-1]
dim(prev.df)

p.outcome = seq(0.05, 0.05, by = 0.05)
shape2.interval = 0.1
shape2.min = seq(4, 8, by = shape2.interval)
shape2.max = shape2.min + shape2.interval
run.df = matrix(c(rep(rep(shape2.min, each = 1), length(p.outcome)), 
                  rep(rep(shape2.max, each = 1), length(p.outcome)), 
                  rep(p.outcome, each = 1 * length(shape2.min))), 
                ncol = 3)
dim(run.df)
head(run.df)
#vectorise
cl = makeCluster(detectCores()-1)
clusterExport(cl, list('c.power','getoutcome', 'analysis.ub', 'speedglm', 'gmm', 'n'))
list = parApply(cl = cl, run.df, 1, function(x) {
  eff.conf.outcome = c.power (n = 505, p.experiment = x[3], p.stdcare =  x[3], adhere.experiment = 1, adhere.stdcare = 1, 
                              confounder.intervention = 'Increase likelihood', confounder.outcome = 'Increase likelihood', NImargin = 0.1, 
                              shape2.min = x[1], 
                              shape2.max = x[2])
  
  return(list(mean (eff.conf.outcome[[3]]), 
              median (eff.conf.outcome[[3]])))
})
stopCluster(cl)

out.df = cbind(run.df, matrix(unlist(list), byrow = T, ncol = 2))
colnames(out.df) = c('shape2.min', 'shape2.max', 'p.outcome', 'coef.mean', 'coef.median')
new.df = rbind(prev.df, out.df)
dim(new.df) > dim(prev.df)
write.csv(new.df,'eff.conf.outcome.csv')

##review results 
out = read.csv('eff.conf.outcome.csv')
out10 = out[which(out[,5] <=10.5),] #remove rows with coeff values more than 10
plot(x=out10[,4], y=out10[,5], axes=FALSE, xlab = 'outcome prob', ylab = 'confounder coeff', pch = 4) #view view values need to be added 
axis(side=1, at=seq(0, 0.5, by = 0.05))
axis(side=2, at=seq(0, 10.5, by = 1))

##make table to be used by shiny 
df = read.csv(file = 'eff.conf.outcome.csv')[,-1]
coef.values = seq(0.25, 10.25, by = 0.5)
p.values = seq(0.05, 0.5, by = 0.05)
out = data.frame(coef.lower = rep(coef.values, each =length(p.values)),
                 coef.upper = rep(coef.values + 0.5, each = length(p.values)),
                 p = rep(p.values, length(coef.values)),
                 shape.min = NA, 
                 shape.max = NA)
index = 1
for (i in 1:length(coef.values)) {
  df.loop.i = df[which(df$coef.mean >= coef.values[i] & df$coef.mean < coef.values[i] + 0.5),]
  for (j in 1:length(p.values)) {
    df.loop.j = df.loop.i[which(as.character(df.loop.i$p.outcome) == as.character(p.values[j])),]
    out$shape.min[index] = min(df.loop.j$shape2.min)
    out$shape.max[index] = max(df.loop.j$shape2.max)
    index = index + 1
  }
}
out = out[order(out$p),]
write.csv(out,'eff.conf.outcome.tab.csv')
View(out)

#.df = df[which(df$coef.mean <10.5),]
#plot(x = .df$shape2.min, y = .df$coef.mean)
