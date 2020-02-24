library(ggplot2)
library(ggpubr)

df = read.csv('powertab.csv')[,1:6]
df.c.cross = as.data.frame(df[df$Confounding == 'confounding' & df$actualtx == 'cross',])
df.nc.cross = as.data.frame(df[df$Confounding == 'nonconfounding ' & df$actualtx == 'cross',])
df.c.inferior = as.data.frame(df[df$Confounding == 'confounding' & df$actualtx == 'inferior',])
df.nc.inferior = as.data.frame(df[df$Confounding == 'nonconfounding ' & df$actualtx == 'inferior',])

legendfontsize=12
alpha = 0.1
width = 20
height = 15

c.cross = ggplot(df.c.cross, aes(x=interval, y=samplesize, colour= Method, linetype=Power)) + 
  geom_point(aes(color = Method), alpha= alpha)+
  scale_alpha_manual(values = c("Inverse probability weighting" = alpha,"Instrumental variable" = alpha), guide = 'none')+
  geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  scale_color_manual(values=c("#00AFB5", "#E69F00"))+
  ylab('')+
  scale_y_continuous(breaks = seq(400, 2000, by = 200), limits = c(350, 2050)) +
  xlab('Proportion of adherent participants in each arm (%)')+
  labs(title = "Non-adherence driven by confounding factors (cross-over)")+
  theme_bw()+
  theme(legend.position = 'bottom', legend.text=element_text(size=legendfontsize))

nc.cross = ggplot(df.nc.cross, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
  geom_point(aes(color = Method), alpha= alpha )+
  geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  scale_color_manual(values=c("#00AFB5", "#E69F00"))+
  ylab('Sample size per group')+
  scale_y_continuous(breaks = seq(400, 2000, by = 200), limits = c(350, 2050)) +
  xlab('Proportion of adherent participants in each arm (%)')+
  labs(title = "Non-adherence driven by non-confounding factors (cross-over)")+
  theme_bw()+
  theme(legend.position = 'bottom',  legend.text=element_text(size=legendfontsize))

c.inferior = ggplot(df.c.inferior, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
  geom_point(aes(color = Method), alpha= alpha)+
  geom_smooth(aes(colour = Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  scale_color_manual(values = "#E69F00")+
  ylab('')+
  scale_y_continuous(breaks = seq(350, 700, by = 50), limits = c(350, 700)) +
  xlab('Proportion of adherent participants in each arm (%)')+
  labs(title = "Non-adherence driven by confounding factors (non-adherent participants take up alternative inferior treatment)")+
  theme_bw()+
  theme(legend.position = 'bottom',  legend.text=element_text(size=legendfontsize))

nc.inferior = ggplot(df.nc.inferior, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
  geom_point(aes(color = Method), alpha= alpha)+
  geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  scale_color_manual(values="#E69F00")+
  scale_y_continuous(breaks = seq(350, 700, by = 50), limits = c(350, 700)) +
  ylab('Sample size per group')+
  xlab('Proportion of adherent participants in each arm (%)')+
  labs(title = "Non-adherence driven by non-confounding factors (non-adherent participants take up alternative inferior treatment)")+
  theme_bw()+
  theme(legend.position = 'bottom',  legend.text=element_text(size=legendfontsize))

plot = ggarrange(nc.cross, c.cross, nc.inferior, c.inferior, 
               labels=c('A', 'B', 'C', 'D'),
               common.legend = T,legend = 'bottom')
ggsave(filename = paste("samplesize", Sys.Date(), ".jpeg"), width = width, height = height*3/4, units = "in")


################################################
###simplified for presentation in conferences###
################################################
Powerdata.nc=data.frame(
  Method=rep(c('Inverse probability weighting','Instrumental variable'), each=length(interval)),
  Power= rep(c('90%'),each=length(interval)),
  interval=interval,
  samplesize=c(630, 605, 580, 560, 525, 505,
               2045, 1410, 1010, 780, 625, 505
  )
)

ggplot(Powerdata.nc, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
  geom_point(aes(color = Method), alpha=0.01)+
  scale_alpha_manual(values = c("Inverse probability weighting"=0.00001,"Instrumental variable" = 0.00001), guide = 'none')+
  geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  scale_color_manual(values=c("#00AFB5", "#E69F00"))+
  ylab('Sample size per group')+
  xlab('Proportion of adherent participants (%)')+
  labs(title = "Non-adherence driven by non-confounding factors")+
  theme_bw()+
  theme(legend.position = 'bottom', legend.text=element_text(size=legendfontsize))

