##NI margin graph 
library(ggplot2)

d= data.frame(type=c('Intention to treat', 'Per protocol', 'Inverse probability weighting', 'Instrumental variable estimation'), 
              estimates=c(0.025, -0.024, 0.114, 0.061), 
              ci.lower=c(-0.044, -0.106, 0.032, -0.107),
              ci.upper=c(0.094, 0.058, 0.195, 0.228))
d.melt=data.table::melt(d, id.var='type')

ggplot()+ 
  geom_point(size=2, aes(x=value, y=type, colour=type), data=d.melt) +
  geom_segment(aes(x = ci.lower, y = type, 
                   xend =ci.upper, yend = type, colour=type), data=d)+
  geom_segment(aes(x = -0.05, y = 4.4, 
                   xend = 0.25, yend = 4.4), arrow = arrow(length = unit(0.1, "inches")))+
  geom_segment(aes(x = 0.25, y = 4.4, 
                   xend = -0.05, yend = 4.4), arrow = arrow(length = unit(0.1, "inches")))+
  annotate('text',x=0.04,y=4.51,label='Favours experimental treatment')+
  annotate('text',x=0.15,y=4.51,label='Favours control treatment')+
  scale_linetype_identity()+
  scale_colour_manual(values =  c("#009E73", "#CC79A7", "#E69F00", "#0072B2"), name='95% Confidence intervals')+
  geom_vline(xintercept=0.1, linetype='dashed', colour='red', size=0.4)+
  labs(y='', x= 'Estimate of treatment difference')+
  scale_y_discrete(limits = c('Instrumental variable estimation', 'Inverse probability weighting', 'Per protocol', 'Intention to treat'))+
  theme_minimal()+
  theme(legend.position = 'bottom')
