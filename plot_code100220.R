################################################################################################################
###################Using causal inference to address non-adherence in non inferiority trials####################
#################################################Plot functions#################################################
################################################################################################################

plot.eff <- function(df, method, nIterations, true.effect, ymin, ymax){
  
  #true effect 
  true.effect=p.experiment-p.stdcare
  
  k = which(method==analysis.method)
  x = c()
  y = c()
  
  for (i in 1:length(interval)) {
    x[[i]] = df[[i]][,k]
    y[[i]] = rep(interval[i],dim(df[[i]])[1])
  }
  
  plotdata = matrix(c(unlist(x),unlist(y)), ncol = 2)
  
  plot = ggplot(as.data.frame(plotdata), aes(x=V2, y=V1))+ 
    geom_point(aes(y=V1, colour=method, alpha=0.0001), size=0.2) + 
    geom_smooth(aes(y=V1, colour=method), method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
    guides(alpha=FALSE)+
    xlab("Proportion of adherent participants in each arm")+
    ylab("Effect estimate")+
    theme_minimal()+
    scale_colour_manual(values=cbPalette[k])+
    theme(legend.title=element_blank(), legend.position="none", legend.text=element_text(size=legendfontsize))+
    scale_x_continuous(limits=c(0.58, 1), breaks = seq(0.6, 1, by=0.1))+
    scale_y_continuous(limits=c(ymin, ymax))+
    geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = 0.25),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
    geom_segment(aes(x = 0.6, y = 0.1, xend = 0.6, yend = -0.05),arrow = arrow(length = unit(0.07, "inches")), colour = 'grey40') + 
    annotate(geom="text", x=0.58, y=0.11, angle = 90, label="Favour control", color="grey40", hjust=0) +
    annotate(geom="text", x=0.58, y=0.09, angle = 90, label="Favour experiment", color="grey40", hjust=1) +
    geom_hline(yintercept=true.effect, linetype='dashed', color='red', size=0.5)
  
  return(plot)
}

plot.eff.multi <- function(estimate, cross.over, ymin, ymax){
  
  analysis.method = c("Instrumental variable","Intention to treat","Inverse probability weighting","Per protocol")
  
  iv = plot.eff(df=estimate, analysis.method[1], nIterations=nIterations, true.effect = true.effect , ymin=ymin, ymax=ymax)
  itt = plot.eff(df=estimate, analysis.method[2], nIterations=nIterations, true.effect = true.effect, ymin=ymin, ymax=ymax)
  mpp = plot.eff(df=estimate, analysis.method[3], nIterations=nIterations, true.effect = true.effect, ymin=ymin, ymax=ymax)
  pp = plot.eff(df=estimate, analysis.method[4], nIterations=nIterations, true.effect = true.effect,  ymin=ymin, ymax=ymax)
  
  if (cross.over == T){
    bias.plot = ggarrange(itt, pp, 
                          mpp, iv,
                          ncol = 2, nrow = 2,
                          common.legend = FALSE )
  } else {
    bias.plot = ggarrange(itt, pp, 
                          mpp,
                          ncol = 2, nrow = 2,
                          common.legend = FALSE )
  }
}

plot.t1 <- function(df, cross.over){
  
  ymin = min(df)
  ymax = max(df)
  if (ymax < 0.1) {ymax=0.1}
  
  if (cross.over==T) {
    plot = ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V1, colour=analysis.method[1], alpha=0.0001)) + 
      geom_point(aes(y=df$V2, colour=analysis.method[2], alpha=0.0001)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
      stat_smooth(aes(y=df$V1, colour=analysis.method[1]), se = FALSE)+
      stat_smooth(aes(y=df$V2, colour=analysis.method[2]),  se = FALSE)+
      stat_smooth(aes(y=df$V3, colour=analysis.method[3]),  se = FALSE)+
      stat_smooth(aes(y=df$V4, colour=analysis.method[4]),  se = FALSE)+
      guides(alpha=FALSE)+
      xlab("Proportion of adherent participants in each arm")+
      ylab("Type 1 error")+
      theme_minimal()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom", legend.text=element_text(size=legendfontsize))+
      theme(legend.title=element_blank())+
      scale_x_continuous(limits=c(start.interval, 1))+
      scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(ymin, ymax))+
      geom_hline(yintercept=0.025, linetype='dashed', color='red', size=0.5)
  } else {
    plot = ggplot(df, aes(interval))+ 
      geom_point(aes(y=df$V2, colour=analysis.method[2], alpha=0.0001)) + 
      geom_point(aes(y=df$V3, colour=analysis.method[3], alpha=0.0001)) + 
      geom_point(aes(y=df$V4, colour=analysis.method[4], alpha=0.0001)) + 
      stat_smooth(aes(y=df$V2, colour=analysis.method[2]), se = FALSE)+
      stat_smooth(aes(y=df$V3, colour=analysis.method[3]), se = FALSE)+
      stat_smooth(aes(y=df$V4, colour=analysis.method[4]), se = FALSE)+
      guides(alpha=FALSE)+
      xlab("Proportion of adherent participants in each arm")+
      ylab("Type 1 error")+
      theme_minimal()+
      scale_colour_manual(values=cbPalette)+
      theme(legend.position="bottom", legend.text=element_text(size=legendfontsize), legend.title=element_blank())+
      scale_color_manual(values = cbPalette[2:4], breaks = analysis.method[2:4]) +
      scale_x_continuous(limits=c(start.interval, 1))+
      scale_y_continuous(breaks=seq(0,1,0.1), limits = c(ymin, ymax))+
      geom_hline(yintercept=0.025, linetype='dashed', color='red', size=0.5)
  }
  return(plot)
}