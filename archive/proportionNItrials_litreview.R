setwd('/Users/moyin/Desktop//')
d<-as.data.frame(read.csv('Book2.csv',stringsAsFactors=FALSE, header = FALSE))
tot<-as.data.frame(table(d[,1]))
ni<-as.data.frame(table(d[,2]))
comb<- merge(tot,ni, by ="Var1")
dates<- as.Date(c("2017-04-01",
                "2017-08-01",
                "2017-02-01",
                "2017-01-01",
                "2017-07-01",
                "2017-06-01",
                "2017-03-01",
                "2017-05-01",
                "2017-11-01",
                "2017-10-01",
                "2017-09-01",
                "2018-04-01",
                "2018-08-01",
                "2018-02-01",
                "2018-01-01",
                "2018-07-01",
                "2018-03-01",
                "2018-05-01",
                "2018-11-01",
                "2018-10-01",
                "2018-09-01",
                "2019-01-01"))
per<-comb$Freq.y/comb$Freq.x*100
comb<-cbind(comb,dates,per)
comb<-comb[order(comb$dates),]
comb.add<- data.frame(Var1=c("2017 Dec", "2018 Jun", "2018 Dec"), 
                      Freq.x=c(0,0,0),
                      Freq.y=c(0,0,0),
                      dates=as.Date(c("2017-12-01","2018-06-01","2018-12-01")),
                                    per=c(0,0,1/22*100))
comb<-rbind(comb,comb.add)
comb<-comb[order(comb$dates),]
plot(x=comb$dates, y=comb$per, xlab='months', ylab='Proportion of NI trials published per month')
abline(lm(comb$per~comb$dates))
