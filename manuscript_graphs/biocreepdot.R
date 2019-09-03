######Biocreep dotplot#######

setwd("/Users/moyin/Desktop/NItrialsimulation-master/") #set working directory 
source('/Users/moyin/Desktop/NItrialsimulation-master/causalinference 11 March 2019 (no AT and matching).R')

#set number of data points at which simulated data is analysed  
start.interval=0.7
interval <- seq(from=start.interval, to=1,by=0.1)

type1.nonconfounding.dotplot<- function(n, nIterations, interval,NImargin, noncomply){  
    
    p.experiment=seq(0.3, 0.7, by=0.1)
    
    #alpha error and critical value 
    z = qnorm(0.975) #alpha=0.025 (one sided)
    
    #make up vectors for simulations 
    itt<-c()  #for saving output from each interval
    estimate<-c()
    type1=list()
    
    #simulate data
    for (j in 1:length(p.experiment)){
        for(i in 1:length(interval)) { print(paste ("interval value",i, "out of", length(interval)))
            for(l in 1:nIterations) {  
                tryCatch({
                    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
                    p.stdcare=p.experiment[j]-NImargin 
                    
                    comply.experiment<-interval
                    comply.stdcare<-interval 
                    
                    id=seq(1,(2*n), by=1) #create participant id  
                    randomisation=c(rep(1,n), rep(0,n)) #randomisation
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    outcome1 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.experiment[j],1-p.experiment[j]))  #probability of outcome if intervention = 1
                    outcome0 = sample(c(1,0), size=2*n, replace=TRUE, prob=c(p.stdcare, 1-p.stdcare))       #probability of outcome if intervention = 0
                    
                    #INTERVENTION dependent on compliance 
                    experiment.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(1-comply.experiment[i],comply.experiment[i]))
                    stdcare.intervention = sample(c(0,1), size=n, replace=TRUE, prob=c(comply.stdcare[i],1-comply.stdcare[i]))
                    intervention = c(experiment.intervention,stdcare.intervention)
                    
                    #ACTUAL OUTCOMES depend on intervention
                    outcome<-getoutcome(vector.outcome1=outcome1, vector.outcome0=outcome0, intervention=intervention)
                    
                    simdata=matrix(data=c(id,randomisation,intervention,outcome), nrow=(2*n))
                    
                    ## intention to treat 
                    pz1.value = mean(simdata[which(simdata[,2]==1),][,4])
                    pz0.value = mean(simdata[which(simdata[,2]==0),][,4])
                    eff.itt = pz1.value-pz0.value
                    
                    ## intention to treat 
                    var.eff.itt<- pz1.value*(1-pz1.value)/n + pz0.value*(1-pz0.value)/n
                    CI.itt<- eff.itt + z*sqrt(var.eff.itt)
                    itt[[l]]<-CI.itt<NImargin
                    
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
            }
            # mean of type 1 error from iterated data 
            estimate[[i]]<-mean(itt)
        }
        
        type1[[j]]=estimate
        
    }
    
    raw=matrix(unlist(type1), ncol=length(interval), byrow=T)
    multi= matrix(round(c(rep(1,4), raw[1,], raw[1,]*raw[2,], raw[1,]*raw[2,]*raw[3,], raw[1,]*raw[2,]*raw[3,]*raw[4,])*100), 
                   ncol=length(interval), byrow=T)
    
    df= data.frame(
        adherence= rep(c("70% Adherence","80% Adherence","90% Adherence","100% Adherence"), each=5),
        trial.label=rep(c("SOC vs A", "A vs B", "B vs C", "C vs D", 'D vs E'), 4),
        Trials=rep(1:5,4),
        y= rep(c(0,-10,-20,-30,-40), 4),
        `Type 1 error`= as.vector(multi)
    )
    
    return(df)
    
}

foo <- data.frame(name = c("SOC vs A", "A vs B", "B vs C", "C vs D"),
                  X = seq(1.5,4.5, by=1))

plotdf=type1.nonconfounding.dotplot(n=500,  nIterations=10000, interval=interval,noncomply='both',NImargin=0.1)

a=ggplot(plotdf[1:5,], aes(Trials, y)) + 
    labs(title=plotdf[1:5,]$adherence) + 
    geom_text(aes(label=paste0(Type.1.error,'%')),hjust=-0.5, vjust=0.5)+
    scale_y_continuous(limits = c(-40, 0),breaks = seq(-40, 0, by = 10))+
    scale_x_continuous(breaks=seq(1.5,4.5, by=1), labels = foo$name, limits = c(1,6))+
    scale_colour_manual(values=c('#3DA09B', '#16BFB6','#5AE2DC', '#35FFF4', '#215B58'))+
    geom_vline(xintercept = c(1,2,3,4,5), color = "grey", size=0.2)+
    geom_hline(yintercept = c(-40,-30, -20, -10, 0), color = "grey", size=0.2)+
    geom_point(aes(col=trial.label, size=Type.1.error))+
    ylab('')+
    xlab('')+
    theme_minimal()+theme(text = element_text(size=15),
                          legend.position="none", panel.grid.major = element_blank(),panel.grid.minor= element_blank(),
                          axis.text.x=element_text(angle = 90, hjust = 1,face="bold"))

b=ggplot(plotdf[6:9,], aes(Trials, y)) + 
    labs(title=plotdf[6:10,]$adherence) + 
    geom_text(aes(label=paste0(Type.1.error,'%')),hjust=-0.5, vjust=0.5)+
    scale_y_continuous(limits = c(-40, 0),breaks = seq(-40, 0, by = 10))+
    scale_x_continuous(breaks=seq(1.5,3.5, by=1), labels = foo$name[1:3], limits = c(1,5))+
    scale_colour_manual(values=c('#3DA09B', '#16BFB6','#5AE2DC', '#215B58'))+
    geom_vline(xintercept = c(1,2,3,4), color = "grey", size=0.2)+
    geom_hline(yintercept = c(-40,-30, -20, -10, 0), color = "grey", size=0.2)+
    geom_point(aes(col=trial.label, size=Type.1.error))+
    ylab('')+
    xlab('')+
    theme_minimal()+theme(text = element_text(size=15),
                          legend.position="none", panel.grid.major = element_blank(),panel.grid.minor= element_blank(),
                          axis.text.x=element_text(angle = 90, hjust = 1,face="bold"))

c=ggplot(plotdf[11:13,], aes(Trials, y)) + 
    labs(title=plotdf[11:15,]$adherence) + 
    geom_text(aes(label=paste0(Type.1.error,'%')),hjust=-0.5, vjust=0.5)+
    scale_y_continuous(limits = c(-40, 0),breaks = seq(-40, 0, by = 10))+
    scale_x_continuous(breaks=seq(1.5,2.5, by=1), labels = foo$name[1:2], limits = c(1,4))+
    scale_colour_manual(values=c('#3DA09B', '#16BFB6', '#215B58'))+
    geom_vline(xintercept = c(1,2,3), color = "grey", size=0.2)+
    geom_hline(yintercept = c(-40,-30, -20, -10, 0), color = "grey", size=0.2)+
    geom_point(aes(col=trial.label, size=Type.1.error))+
    ylab('')+
    xlab('')+
    theme_minimal()+theme(text = element_text(size=15),
                          legend.position="none", panel.grid.major = element_blank(),panel.grid.minor= element_blank(),
                          axis.text.x=element_text(angle = 90, hjust = 1, face="bold"))

d=ggplot(plotdf[16:17,], aes(Trials, y)) + 
    labs(title=plotdf[16:17,]$adherence) + 
    geom_text(aes(label=paste0(Type.1.error,'%')),hjust=-0.5, vjust=0.5)+
    scale_y_continuous(limits = c(-40, 0),breaks = seq(-40, 0, by = 10))+
    xlab('')+
    scale_x_continuous(breaks=seq(1.5,1.5, by=1), labels = foo$name[1], limits = c(1,3))+
    scale_colour_manual(values=c('#3DA09B', '#215B58'))+
    geom_vline(xintercept = c(1,2), color = "grey", size=0.2)+
    geom_hline(yintercept = c(-40,-30, -20, -10, 0), color = "grey", size=0.2)+
    geom_point(aes(col=trial.label, size=Type.1.error))+
    ylab('Decrease in treatment efficacy compared to initial standard-of-care (%)')+
    theme_minimal()+theme(text = element_text(size=15),
                          legend.position="none", panel.grid.major = element_blank(),panel.grid.minor= element_blank(),
                          axis.text.x=element_text(angle = 90, hjust = 1, face="bold"))

p=ggarrange(d,c,b,a, ncol=4)

annotate_figure(p,
                bottom = text_grob("Consecutive non-inferiority trials")
)
