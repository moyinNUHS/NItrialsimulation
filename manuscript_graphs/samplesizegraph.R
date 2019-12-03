library(ggplot2)
library(ggpubr)

interval=seq(75,100,by=5)
legendfontsize=12

Powerdata.c=data.frame(
    Method=rep(c('Inverse probability weighting','Instrumental variable'), each=2*length(interval)),
    Power= rep(c('90%', "80%"),each=length(interval)),
    interval=rep(interval,2),
    samplesize=c(760, 710, 660, 620, 550, 500, 
                 590, 550, 515, 470, 420, 370, 
                 1600, 1155, 830, 700, 550, 500, 
                 1300, 880, 720, 525, 420, 370)
)

Powerdata.nc=data.frame(
    Method=rep(c('Inverse probability weighting','Instrumental variable'), each=2*length(interval)),
    Power= rep(c('90%', "80%"),each=length(interval)),
    interval=rep(interval,2),
    samplesize=c(630, 605, 580, 560, 525, 505, 
                 495, 470, 450, 430, 410, 370, 
                 2045, 1410, 1010, 780, 625, 505, 
                 1500, 1070, 785, 500, 465, 370
    )
)

nc=ggplot(Powerdata.nc, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
    geom_point(aes(color = Method), alpha=0.01)+
    scale_alpha_manual(values = c("Inverse probability weighting"=0.00001,"Instrumental variable" = 0.00001), guide = 'none')+
    geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    scale_color_manual(values=c("#00AFB5", "#E69F00"))+
    ylab('Sample size per group')+
    xlab('Proportion of adherent participants (%)')+
    labs(title = "Non-adherence driven by non-confounding factors")+
    theme_bw()+
    theme(legend.position = 'bottom', legend.text=element_text(size=legendfontsize))

c=ggplot(Powerdata.c, aes(x=interval, y=samplesize, colour=Method, linetype=Power)) +
    geom_point(aes(color = Method), alpha=0.01)+
    geom_smooth(aes(colour=Method), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    scale_color_manual(values=c("#00AFB5", "#E69F00"))+
    ylab('')+
    xlab('Proportion of adherent participants (%)')+
    labs(title = "Non-adherence driven by confounding factors")+
    theme_bw()+
    theme(legend.position = 'bottom',  legend.text=element_text(size=legendfontsize))

plot=ggarrange(nc,c,labels=c('A','B'),
          common.legend = T,legend = 'bottom')
ggsave(filename = paste("samplesize", Sys.Date(), ".pdf"), width = width, height = height*3/4, units = "in")

###simplified for presentation
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

