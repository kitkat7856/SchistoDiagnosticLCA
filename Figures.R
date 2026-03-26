
### FIGURES ####

#### RAW DATA ####
kkmean<- apply(kk, c(1,2), mean, na.rm=T)
plot(kkmean, ct, xlab="mean egg counts", ylab="stool qpcr CT")


#### Boxplot for infection intensities #####
boxplot(as.numeric(dbspcr$CT), main="DBS qPCR CT")
dbs.df<- data.frame(-dbs)
colnames(dbs.df)<- c("PreTx", "3 weeks", "9 weeks", "6 months")
boxplot(dbs.df, xlab="timepoints", main="DBS qPCR CT")

ggplot(dbspcr, aes(x=factor(tmp, levels=c("preT", "3wk", "9wk", "6mo")), y=-as.numeric(CT))) +
  geom_boxplot() +
  theme_bw()

ct.df<- data.frame(-ct)
colnames(ct.df)<- c("PreTx", "3 weeks", "9 weeks", "6 months")
boxplot(ct.df, xlab="timepoints", main="Stool qPCR CT")

g.df<- data.frame(gscore+1)
colnames(g.df)<- c("PreTx", "3 weeks", "9 weeks", "6 months")
boxplot(g.df, xlab="timepoints", main="G-score")

kk.df<- data.frame(sqrt(kkmean))
colnames(kk.df)<- c("PreTx", "3 weeks", "9 weeks", "6 months")
boxplot(kk.df, xlab="timepoints", main="KK infection intensity", ylab="Square root egg counts", ylim=c(0,8))

### raw prev ####
KK.prev<- data.frame(PreTx = mean(kkmean[,1]>0, na.rm=T), 
                     Wk3 = mean(kkmean[,2]>0, na.rm=T),
                     Wk9 = mean(kkmean[,3]>0, na.rm=T),
                     Mnth6 = mean(kkmean[,4]>0, na.rm=T))

cca.prev<- data.frame(PreTx = mean(gscore[,1]>=2, na.rm=T), 
                      Wk3 = mean(gscore[,2]>=2, na.rm=T),
                      Wk9 = mean(gscore[,3]>=2, na.rm=T),
                      Mnth6 = mean(gscore[,4]>=2, na.rm=T))

qpcr.prev<- data.frame(PreTx = mean(ct[,1]<=38, na.rm=T), 
                       Wk3 = mean(ct[,2]<=38, na.rm=T),
                       Wk9 = mean(ct[,3]<=38, na.rm=T),
                       Mnth6 = mean(ct[,4]<=38, na.rm=T))

obsstatus<- matrix(NA,210,4)
obsstatus[kkmean==0 | ct>38 | gscore<2]<-0
obsstatus[kkmean>0 | ct<=38 | gscore>=2]<-1
colMeans(obsstatus, na.rm = T)

obs.prev<- rbind(KK=KK.prev,CCA=cca.prev,qPCR=qpcr.prev, dbs.prev= colMeans(dbs.ct<38,na.rm=T), all=colMeans(obsstatus, na.rm = T))
obs.prev<-round(obs.prev*100,2)

### raw sums pos ####
KK.sum<- data.frame(PreTx = sum(kkmean[,1]>0, na.rm=T), 
                    Wk3 = sum(kkmean[,2]>0, na.rm=T),
                    Wk9 = sum(kkmean[,3]>0, na.rm=T),
                    Mnth6 = sum(kkmean[,4]>0, na.rm=T))

cca.sum<- data.frame(PreTx = sum(gscore[,1]>=2, na.rm=T), 
                     Wk3 = sum(gscore[,2]>=2, na.rm=T),
                     Wk9 = sum(gscore[,3]>=2, na.rm=T),
                     Mnth6 = sum(gscore[,4]>=2, na.rm=T))

qpcr.sum<- data.frame(PreTx = sum(ct[,1]<=38, na.rm=T), 
                      Wk3 = sum(ct[,2]<=38, na.rm=T),
                      Wk9 = sum(ct[,3]<=38, na.rm=T),
                      Mnth6 = sum(ct[,4]<=38, na.rm=T))

dbs.sum<- data.frame(PreTx = sum(ct[,1]<=38, na.rm=T), 
                     Wk3 = sum(ct[,2]<=38, na.rm=T),
                     Wk9 = sum(ct[,3]<=38, na.rm=T),
                     Mnth6 = sum(ct[,4]<=38, na.rm=T))

obsstatus<- matrix(NA,210,4)
obsstatus[kkmean==0 | ct>38 | gscore<2]<-0
obsstatus[kkmean>0 | ct<=38 | gscore>=2]<-1

obs.prev1<- rbind(KK=KK.sum,CCA=cca.sum,qPCR=qpcr.sum, all=colSums(obsstatus, na.rm = T))

###### logistic curves #####
eggs<- seq(0,100, by=0.1)

#gscore mean logistic curve
k<-as.numeric(point_estimate(out$k, "median"))
int<-as.numeric(point_estimate(out$int, "median"))
d<- data.frame(x=eggs,y=9/(1+exp(-k*(eggs-int))) +1) #+1 for  gscore
d$ymin<- d$y- 2*(sqrt(1/1.09)) #-2sd
d$ymax<- d$y+ 2*(sqrt(1/1.09)) #+2sd
d$ymax[d$ymax>10]<-10
ggplot(d, aes(x=x*24,y=y)) +
  geom_line() +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.5, fill="orange") +
  scale_y_continuous(breaks = seq(1,10), limits = c(1,10)) +
  xlim(0,1500) +
  xlab("Eggs per gram") +
  ylab("G-score") +
  theme_bw()

##### KK to stool CT #####
a<- as.numeric(point_estimate(out$a_ct, "median"))
b<- as.numeric(point_estimate(out$b_ct, "median"))
d<- data.frame(x=eggs, y= -a*log(eggs+0.01)+b)
d$ymin<-  d$y- 2*(sqrt(1/0.08)) #-2sd
d$ymax<-  d$y + 2*(sqrt(1/0.08)) #+2sd
z<- data.frame(y=c(ct),x=c(kkmean))

ggplot(d, aes(x=x*24,y=y)) +
  geom_line() +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.5, fill="blue") +
  xlab("Eggs per gram") +
  ylab("Stool qPCR CT") +
  geom_point(data=z)+
  theme_bw()


##### Gscore logcurve taking all posteriors  ######
x<- sample(40000,1000)
k<- as.numeric(c(mod$mcmc[[1]][,"k"], mod$mcmc[[2]][,"k"]))[x]
int<- as.numeric(c(mod$mcmc[[1]][,"int"], mod$mcmc[[2]][,"int"]))[x]
d2<- data.frame()
for(i in 1:length(x)){
  d<- data.frame(x=eggs, y=9/(1+exp(-k[i]*(eggs-int[i]))) +1)
  d2<- rbind(d2,d)
}
d2$ymin<- d2$y- 2*(sqrt(1/1.09))
d2$ymax<- d2$y+ 2*(sqrt(1/1.09))

gscore_logcurve<- d2 %>%
  group_by(x) %>%
  summarise(meanline= quantile(y, 0.5),
            ymin=quantile(ymin, 0.5),
            ymax= quantile(ymax, 0.5))
gscore_logcurve$ymax[gscore_logcurve$ymax>10]<-10

ggplot(gscore_logcurve, aes(x=x*24,y=meanline)) +
  geom_line() +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.5, fill="orange") +
  scale_y_continuous(breaks = seq(1,10,by=1), limits=c(1,10)) +
  xlab("Eggs per gram") +
  ylab("G-score") +
  theme_bw()

##### DBS logcurve taking all posteriors 
x<- sample(40000,1000)
k2<- as.numeric(c(mod$mcmc[[1]][,"k2"], mod$mcmc[[2]][,"k2"]))[x]
int2<- as.numeric(c(mod$mcmc[[1]][,"int2"], mod$mcmc[[2]][,"int2"]))[x]
d2<- data.frame()
for(i in 1:length(x)){
  d<- data.frame(x=eggs, y=1/(1+exp(-k2[i]*(eggs-int2[i]))))
  d2<- rbind(d2,d)
}
dbs_logcurve<- d2 %>%
  group_by(x) %>%
  summarise(meanline= quantile(y, 0.5),
            ymin=quantile(y, 0.025),
            ymax= quantile(y, 0.975))

ggplot(dbs_logcurve, aes(x=x*24,y=meanline)) +
  geom_line() +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.5, fill="darkred") +
  scale_y_continuous(breaks = seq(0,1,by=0.2)) +
  xlab("Eggs per gram") +
  ylab("Probability of DBS") +
  theme_bw()

### Status at each timepoint #####
timepoints<- c("PreT", "3wk", "9wk", "6mo")
prob<- select(out, starts_with("prob["))
prob.samples<- prob[sample(40000, 500),]
x<- melt(prob.samples)
x$tmp<- NA
x$tmp[str_detect(x$variable, "1,")]<- "PreT"
x$tmp[str_detect(x$variable, "2,")]<- "3wk"
x$tmp[str_detect(x$variable, "3,")]<- "9wk"
x$tmp[str_detect(x$variable, "4,")]<- "6mo"

x$status<-NA
x$status[str_detect(x$variable, ",1")]<- "1"
x$status[str_detect(x$variable, ",2")]<- "2"
x$status[str_detect(x$variable, ",3")]<- "3"

prev<- x%>%
  filter(status == "1") %>%
  group_by(tmp) %>%
  summarise(prev= (1-median(value)))

ggplot(x, aes(x=factor(tmp, level=timepoints), y=value, color=status)) +
  geom_jitter(alpha=1/5) +
  stat_summary(geom = "point", fun = "median", size=3, aes(fill=status), colour="black", shape=21) +
  geom_point(data=prev, aes(x=tmp, y=prev), color="red", shape= 4, size=5) +
  xlab("Timepoint") +
  ylab("Proportion of population")+
  theme_bw()

z<-x[x$status==3,]
ggplot(z, aes(x=factor(tmp, level=timepoints), y=value))+
  geom_jitter(alpha=1/2, color="#619CFF") +
  stat_summary(geom = "point", fun = "median", size=3, colour="black", shape=4) +
  xlab("Timepoint") +
  ylab("Proportion of population shedding")+
  theme_bw()

#######
status<- select(out, starts_with("Status"))
status.samples<- status[sample(40000,1000),]
status1<- status.samples[,1:210]
status2<- status.samples[,211:420]
status3<- status.samples[,421:630]
status4<- status.samples[,631:840]

prev1<- data.frame(notshed=rowSums(status1==2)/rowSums(status1==2 | status1==3), shed=rowSums(status1==3)/rowSums(status1==2 | status1==3), shed.t=rowSums(status1==3)/210, total=rowSums(status1==2 | status1==3)/210, tmp="PreT")
prev2<- data.frame(notshed=rowSums(status2==2)/rowSums(status2==2 | status2==3), shed=rowSums(status2==3)/rowSums(status2==2 | status2==3), shed.t=rowSums(status2==3)/210,total=rowSums(status2==2 | status2==3)/210, tmp="3wk")
prev3<- data.frame(notshed=rowSums(status3==2)/rowSums(status3==2 | status3==3), shed=rowSums(status3==3)/rowSums(status3==2 | status3==3), shed.t=rowSums(status3==3)/210,total=rowSums(status3==2 | status3==3)/210, tmp="9wk")
prev4<- data.frame(notshed=rowSums(status4==2)/rowSums(status4==2 | status4==3), shed=rowSums(status4==3)/rowSums(status4==2 | status4==3), shed.t=rowSums(status4==3)/210,total=rowSums(status4==2 | status4==3)/210, tmp="6mo")
prev.df<- rbind(prev1,prev2,prev3,prev4)

prev<- matrix(NA, 3,4)
prev[1,1]<- paste0(round(median(prev1$total)*100,1),"%"," ","(",round(hdi(prev1$total)$CI_low*100,1),",",round(hdi(prev1$total)$CI_high*100,1),")")
prev[1,2]<- paste0(round(median(prev2$total)*100,1),"%"," ","(",round(hdi(prev2$total)$CI_low*100,1),",",round(hdi(prev2$total)$CI_high*100,1),")")
prev[1,3]<- paste0(round(median(prev3$total)*100,1),"%"," ","(",round(hdi(prev3$total)$CI_low*100,1),",",round(hdi(prev3$total)$CI_high*100,1),")")
prev[1,4]<- paste0(round(median(prev4$total)*100,1),"%"," ","(",round(hdi(prev4$total)$CI_low*100,1),",",round(hdi(prev4$total)$CI_high*100,1),")")

prev[2,1]<- paste0(round(median(prev1$shed.t)*100,1),"%"," ","(",round(hdi(prev1$shed.t)$CI_low*100,1),",",round(hdi(prev1$shed.t)$CI_high*100,1),")")
prev[2,2]<- paste0(round(median(prev2$shed.t)*100,1),"%"," ","(",round(hdi(prev2$shed.t)$CI_low*100,1),",",round(hdi(prev2$shed.t)$CI_high*100,1),")")
prev[2,3]<- paste0(round(median(prev3$shed.t)*100,1),"%"," ","(",round(hdi(prev3$shed.t)$CI_low*100,1),",",round(hdi(prev3$shed.t)$CI_high*100,1),")")
prev[2,4]<- paste0(round(median(prev4$shed.t)*100,1),"%"," ","(",round(hdi(prev4$shed.t)$CI_low*100,1),",",round(hdi(prev4$shed.t)$CI_high*100,1),")")

prev[3,1]<- paste0(round(median(prev1$shed)*100,1),"%"," ","(",round(hdi(prev1$shed)$CI_low*100,1),",",round(hdi(prev1$shed)$CI_high*100,1),")")
prev[3,2]<- paste0(round(median(prev2$shed)*100,1),"%"," ","(",round(hdi(prev2$shed)$CI_low*100,1),",",round(hdi(prev2$shed)$CI_high*100,1),")")
prev[3,3]<- paste0(round(median(prev3$shed)*100,1),"%"," ","(",round(hdi(prev3$shed)$CI_low*100,1),",",round(hdi(prev3$shed)$CI_high*100,1),")")
prev[3,4]<- paste0(round(median(prev4$shed)*100,1),"%"," ","(",round(hdi(prev4$shed)$CI_low*100,1),",",round(hdi(prev4$shed)$CI_high*100,1),")")


ggplot(prev.df, aes(x=factor(tmp, level=timepoints), y=shed)) +
  geom_point(position = "jitter", alpha=0.5, col="#619CFF") +
  stat_summary(geom = "point", fun="median", col="black", shape = 4) +
  scale_y_continuous(breaks = seq(0,1,by=0.2), limits = c(0,1)) +
  xlab("Timepoints") +
  ylab("Proportion of infected that are shedding eggs") +
  theme_bw()

ggplot(prev.df, aes(x=factor(tmp, level=timepoints), y=shed.t)) +
  geom_point(position = "jitter", alpha=0.5, col="#619CFF") +
  stat_summary(geom = "point", fun="median", col="black", shape = 4) +
  scale_y_continuous(breaks = seq(0,1,by=0.2), limits = c(0,1)) +
  xlab("Timepoints") +
  ylab("Proportion of population shedding eggs") +
  theme_bw()

ggplot(prev.df, aes(x=factor(tmp, level=timepoints), y=total)) +
  geom_point(position = "jitter", alpha=0.5, col="red") +
  stat_summary(geom = "point", fun="median", col="black", shape = 4) +
  scale_y_continuous(breaks = seq(0,1,by=0.2), limits = c(0,1)) +
  xlab("Timepoints") +
  ylab("Estimated prevalence") +
  theme_bw()

#######

#### ROC Curves ####
#### ROC CT ####
#status 1 = uninfected(1), status 2 and 3 = infected (2)
status.samples[status.samples==3]<-2

#split status by timepoint, each column is a posterior sample, pair with corresponding CT values
status.1<- data.frame(CT=-ct[,1], t(status.samples[,1:210])) #timepoint 1
status.2<- data.frame(CT=-ct[,2], t(status.samples[,211:420]))#timepoint 2
status.3<- data.frame(CT=-ct[,3], t(status.samples[,421:630])) #timepoint 3
status.4<- data.frame(CT=-ct[,4], t(status.samples[,631:840]))#timepoint 4


#### pre T ROC 
#removeNA
status.1<- status.1[!is.na(status.1$CT),]
#predictions for ROCR- predictions= CT data, list of 500 for each sample
preTpred<- rep(list(status.1$CT),500)
#labels for ROCR- labels = status
preTlabels<- as.list(status.1[,2:501])

pred.1<- prediction(preTpred, preTlabels)
roc.ct.1<-performance(pred.1, "tpr", "fpr")
ss.ct.1<-performance(pred.1, measure="sens", x.measure = "spec")

#### 3 week ROC
#removeNA
status.2<- status.2[!is.na(status.2$CT),]
#predictions for ROCR- predictions= CT data, list of 500 for each sample
wk3pred<- rep(list(status.2$CT),500)
#labels for ROCR- labels = status
wk3labels<- as.list(status.2[,2:501])

pred.2<- prediction(wk3pred, wk3labels)
roc.ct.2<-performance(pred.2, "tpr","fpr")
ss.ct.2<-performance(pred.2, "sens", "spec")
plot(ss.ct.2)
plot(roc.ct.2)

#### 9 week ROC
#removeNA
status.3<- status.3[!is.na(status.3$CT),]
#predictions for ROCR- predictions= CT data, list of 500 for each sample
wk9pred<- rep(list(status.3$CT),500)
#labels for ROCR- labels = status
wk9labels<- as.list(status.3[,2:501])

pred.3<- prediction(wk9pred, wk9labels)
roc.ct.3<-performance(pred.3, "tpr",  "fpr")
ss.ct.3<-performance(pred.3, "sens", "spec")
plot(ss.ct.3)
plot(roc.ct.3)

### 6 month ROC ####
status.4<- status.4[!is.na(status.4$CT),]
#predictions for ROCR- predictions= CT data, list of 500 for each sample
mnth6pred<- rep(list(status.4$CT),500)
#labels for ROCR- labels = status
mnth6labels<- as.list(status.4[,2:501])

#roc
pred.4<- prediction(mnth6pred, mnth6labels)
roc.ct.4<-performance(pred.4, "tpr", "fpr")
ss.ct.4<-performance(pred.4, "sens", "spec")
plot(ss.ct.4)
plot(roc.ct.4)

### plot and save roc final figures ####
cut1<- getoptcut(ss.ct.1,500)
cut2<- getoptcut(ss.ct.2,500)
cut3<- getoptcut(ss.ct.3,500)
cut4<- getoptcut(ss.ct.4,500)

ctstatus<- data.frame(ct=c(-ct), t(status.samples))
ctstatus<- ctstatus[!is.na(ctstatus$ct),]
pred<- rep(list(ctstatus$ct),500)
label<- as.list(ctstatus[,2:501])
pred.ct<- prediction(pred,label)
perf<- performance(pred.ct,"sens", "spec")


pdf("ROC CT curves.pdf")
par(mfrow=c(2,2))
plot(roc.ct.1,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-30,-35,-40,-45),
     print.cutoffs.at = c(-30,-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(-0.2, 1),
     lwd=3,
     col= "blue",
     main = "Pre-Treatment, Stool qPCR")
abline(0,1)
plot(roc.ct.2,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-30,-35,-40,-45),
     print.cutoffs.at = c(-30,-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(-0.2, 1.2),
     lwd=3,
     col= "blue",
     main = "Three weeks, Stool qPCR")
abline(0,1)
plot(roc.ct.3,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-30,-35,-40,-45),
     print.cutoffs.at = c(-30,-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "blue",
     main = "nine weeks, Stool qPCR")
abline(0,1)
plot(roc.ct.4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-30,-40,-45, round(cut4$cut2)),
     print.cutoffs.at = c(-30,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "blue",
     main = "six months, Stool qPCR")
abline(0,1)

dev.off()


#### ROC KK ####
#data frames separating gscore values by each timepoint, status for each tmp
status.1<- data.frame(kk=kkmean[,1], t(status.samples[,1:210])) #timepoint 1
status.2<- data.frame(kk=kkmean[,2], t(status.samples[,211:420]))#timepoint 2
status.3<- data.frame(kk=kkmean[,3], t(status.samples[,421:630])) #timepoint 3
status.4<- data.frame(kk=kkmean[,4], t(status.samples[,631:840]))#timepoint 4

#removeNA
status.1<- status.1[!is.na(status.1$kk),]
#predictions for ROCR- predictions= gscore data, list of 500 for each sample
preTpred<- rep(list(status.1$kk),500)
#labels for ROCR- labels = statu s
preTlabels<- as.list(status.1[,2:501])

pred.1<- prediction((preTpred),(preTlabels))
roc.kk.1<-performance(pred.1, measure="tpr", x.measure = "fpr")
ss.kk.1<-performance(pred.1, measure="sens", x.measure = "spec")

### 3 week ROC
status.2<- status.2[!is.na(status.2$kk),]
#predictions for ROCR- predictions= kkscore data, list of 500 for each sample
wk3pred<- rep(list(status.2$kk),500)
#labels for ROCR- labels = status
wk3labels<- as.list(status.2[,2:501])

pred.2<- prediction(wk3pred, wk3labels)
roc.kk.2<-performance(pred.2, measure="tpr", x.measure = "fpr")
ss.kk.2<-performance(pred.2, measure="sens", x.measure = "spec")
plot(ss.kk.2)
plot(roc.kk.2)

### 9 week ROC
status.3<- status.3[!is.na(status.3$kk),]
#predictions for ROCR- predictions= kkscore data, list of 500 for each sample
wk9pred<- rep(list(status.3$kk),500)
#labels for ROCR- labels = status
wk9labels<- as.list(status.3[,2:501])

pred.3<- prediction(wk9pred, wk9labels)
roc.kk.3<-performance(pred.3, measure="tpr", x.measure = "fpr")
ss.kk.3<-performance(pred.3, measure="sens", x.measure = "spec")
plot(ss.kk.3)
plot(roc.kk.3)

### 6 month ROC
status.4<- status.4[!is.na(status.4$kk),]
#predictions for ROCR- predictions= kkscore data, list of 500 for each sample
mnth6pred<- rep(list(status.4$kk),500)
#labels for ROCR- labels = status
mnth6labels<- as.list(status.4[,2:501])

pred.4<- prediction(mnth6pred, mnth6labels)
roc.kk.4<-performance(pred.4, measure="tpr", x.measure = "fpr")
ss.kk.4<-performance(pred.4, measure="sens", x.measure = "spec")
plot(ss.kk.4)
plot(roc.kk.4)

### plot and save roc final figures ####

pdf("ROC KK curves.pdf")
par(mfrow=c(2,2))
plot(roc.kk.1,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     text.adj = c(-0.5, 0.5),
     lwd=3,
     col= "brown",
     main = "Pre-Treatment, KK")
abline(0,1)
plot(roc.kk.2,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "brown",
     main = "three weeks, KK") 
abline(0,1)
plot(roc.kk.3,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "brown",
     main = "nine weeks, KK")
abline(0,1)
plot(roc.kk.4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "brown",
     main = "six months, KK")
abline(0,1)

dev.off()

##### ROC GSCORE ######
#data frames separating gscore values by each timepoint, status for each tmp
status.1<- data.frame(G=gscore[,1]+1, t(status.samples[,1:210])) #timepoint 1
status.2<- data.frame(G=gscore[,2]+1, t(status.samples[,211:420]))#timepoint 2
status.3<- data.frame(G=gscore[,3]+1, t(status.samples[,421:630])) #timepoint 3
status.4<- data.frame(G=gscore[,4]+1, t(status.samples[,631:840]))#timepoint 4

#removeNA
status.1<- status.1[!is.na(status.1$G),]
#predictions for ROCR- predictions= gscore data, list of 500 for each sample
preTpred<- rep(list(status.1$G),500)
#labels for ROCR- labels = status
preTlabels<- as.list(status.1[,2:501])

pred.1<- prediction((preTpred),(preTlabels))
roc.G.1<-performance(pred.1, measure="tpr", x.measure = "fpr")
ss.G.1<-performance(pred.1, measure="sens", x.measure = "spec")
plot(ss.G.1)
plot(roc.G.1)


### 3 week ROC
status.2<- status.2[!is.na(status.2$G),]
#predictions for ROCR- predictions= gscore data, list of 500 for each sample
wk3pred<- rep(list(status.2$G),500)
#labels for ROCR- labels = status
wk3labels<- as.list(status.2[,2:501])

pred.2<- prediction(wk3pred, wk3labels)
roc.G.2<-performance(pred.2, measure="tpr", x.measure = "fpr")
ss.G.2<-performance(pred.2, measure="sens", x.measure = "spec")
plot(ss.G.2)
plot(roc.G.2)

### 9 week ROC
status.3<- status.3[!is.na(status.3$G),]
#predictions for ROCR- predictions= gscore data, list of 500 for each sample
wk9pred<- rep(list(status.3$G),500)
#labels for ROCR- labels = status
wk9labels<- as.list(status.3[,2:501])

pred.3<- prediction(wk9pred, wk9labels)
roc.G.3<-performance(pred.3, measure="tpr", x.measure = "fpr")
ss.G.3<-performance(pred.3, measure="sens", x.measure = "spec")
plot(ss.G.3)
plot(roc.G.3)

### 6 month ROC
status.4<- status.4[!is.na(status.4$G),]
#predictions for ROCR- predictions= gscore data, list of 500 for each sample
mnth6pred<- rep(list(status.4$G),500)
#labels for ROCR- labels = status
mnth6labels<- as.list(status.4[,2:501])

pred.4<- prediction(mnth6pred, mnth6labels)
roc.G.4<-performance(pred.4, measure="tpr", x.measure = "fpr")
ss.G.4<-performance(pred.4, measure="sens", x.measure = "spec")
plot(ss.G.4)
plot(roc.G.4)

### plot and save roc final figures ####
cut1<- getoptcut(ss.G.1,1)
cut2<- getoptcut(ss.G.2,500)
cut3<- getoptcut(ss.G.3,500)
cut4<- getoptcut(ss.G.4,500)

pdf("ROC Gscore curves.pdf")
par(mfrow=c(2,2))
plot(roc.G.1,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(1,2,3,4,5,6,7,8,9,10),
     print.cutoffs.at = c(1,2,3,4,5,6,7,8,9,10),
     text.adj = c(-0.5, 0.5),
     lwd=3,
     col= "orange",
     main = "Pre-Treatment,Gscore")
abline(0,1)
plot(roc.G.2,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     how.spread.at = c(1,2,3,4,5,6,7,8,9,10),
     print.cutoffs.at = c(1,2,3,4,5,6,7,8,9,10),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "orange",
     main = "three weeks,Gscore") 
abline(0,1)
plot(roc.G.3,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     how.spread.at = c(1,2,3,4,5,6,7,8,9,10, round(cut3$cut2)),
     print.cutoffs.at = c(1,2,3,4,5,6,7,8,9,10, round(cut3$cut2)),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "orange",
     main = "nine weeks, Gscore")
abline(0,1)
plot(roc.G.4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at =  c(1,2,3,4,5,6,7,8,9,10, round(cut3$cut2)),
     print.cutoffs.at =  c(1,2,3,4,5,6,7,8,9,10, round(cut3$cut2)),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "orange",
     main = "six months, Gscore")
abline(0,1)

dev.off()



##### ROC DBS pcr######
#data frames separating pcr values by each timepoint, status for each tmp
status.1<- data.frame(pcr=-dbs.ct[,1], t(status.samples[,1:210])) #timepoint 1
status.2<- data.frame(pcr=-dbs.ct[,2], t(status.samples[,211:420]))#timepoint 2
status.3<- data.frame(pcr=-dbs.ct[,3], t(status.samples[,421:630])) #timepoint 3
status.4<- data.frame(pcr=-dbs.ct[,4], t(status.samples[,631:840]))#timepoint 4

#removeNA
status.1<- status.1[!is.na(status.1$pcr),]
#predictions for ROCR- predictions= pcr data, list of 500 for each sample
preTpred<- rep(list(status.1$pcr),500)
#labels for ROCR- labels = status
preTlabels<- as.list(status.1[,2:501])

pred.1<- prediction((preTpred),(preTlabels))
roc.pcr.1<-performance(pred.1, measure="tpr", x.measure = "fpr")
ss.pcr.1<-performance(pred.1, measure="sens", x.measure = "spec")
plot(ss.pcr.1)
plot(roc.pcr.1)


### 3 week ROC
status.2<- status.2[!is.na(status.2$pcr),]
#predictions for ROCR- predictions= pcr data, list of 500 for each sample
wk3pred<- rep(list(status.2$pcr),500)
#labels for ROCR- labels = status
wk3labels<- as.list(status.2[,2:501])

pred.2<- prediction(wk3pred, wk3labels)
roc.pcr.2<-performance(pred.2, measure="tpr", x.measure = "fpr")
ss.pcr.2<-performance(pred.2, measure="sens", x.measure = "spec")
plot(ss.pcr.2)
plot(roc.pcr.2)

### 9 week ROC
status.3<- status.3[!is.na(status.3$pcr),]
#predictions for ROCR- predictions= pcr data, list of 500 for each sample
wk9pred<- rep(list(status.3$pcr),500)
#labels for ROCR- labels = status
wk9labels<- as.list(status.3[,2:501])

pred.3<- prediction(wk9pred, wk9labels)
roc.pcr.3<-performance(pred.3, measure="tpr", x.measure = "fpr")
ss.pcr.3<-performance(pred.3, measure="sens", x.measure = "spec")
plot(ss.pcr.3)
plot(roc.pcr.3)

### 6 month ROC
status.4<- status.4[!is.na(status.4$pcr),]
#predictions for ROCR- predictions= pcr data, list of 500 for each sample
mnth6pred<- rep(list(status.4$pcr),500)
#labels for ROCR- labels = status
mnth6labels<- as.list(status.4[,2:501])

pred.4<- prediction(mnth6pred, mnth6labels)
roc.pcr.4<-performance(pred.4, measure="tpr", x.measure = "fpr")
ss.pcr.4<-performance(pred.4, measure="sens", x.measure = "spec")
plot(ss.pcr.4)
plot(roc.pcr.4)

### plot and save roc final figures ####

pdf("ROC DBS pcr curves.pdf")
par(mfrow=c(2,2))
plot(roc.pcr.1,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-35,-40,-45),
     print.cutoffs.at = c(-35,-40,-45),
     text.adj = c(0.5,2),
     lwd=3,
     col= "red",
     main = "Pre-Treatment, DBS qPCR")
abline(0,1)
plot(roc.pcr.2,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     how.spread.at = c(-35,-40,-45),
     print.cutoffs.at = c(-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "red",
     main = "three weeks, DBS qpcr") 
abline(0,1)
plot(roc.pcr.3,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     how.spread.at = c(-35,-40,-45),
     print.cutoffs.at = c(-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "red",
     main = "nine weeks, DBS qpcr")
abline(0,1)
plot(roc.pcr.4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at =  c(-35,-40,-45),
     print.cutoffs.at =  c(-35,-40,-45),
     cutoff.label.function = abs,
     text.adj = c(0.75, 2),
     lwd=3,
     col= "red",
     main = "six months, DBS qpcr")
abline(0,1)

dev.off()


sesp<-matrix(NA,4,8)
sesp[1,]<- as.numeric(c(getss(ss.kk.1,500,0.1)[,1:2],getss(ss.kk.2,500,0.1)[,1:2],getss(ss.kk.3,500,0.15)[,1:2],getss(ss.kk.4,500,0.1)[,1:2]))
sesp[2,]<- as.numeric(c(getss(ss.ct.1,500,-40)[,1:2],getss(ss.ct.2,500,-40)[,1:2],getss(ss.ct.3,500,-40)[,1:2],getss(ss.ct.4,500,-40)[,1:2]))
sesp[3,]<- as.numeric(c(getss(ss.G.1,500,3)[,1:2],getss(ss.G.2,500,3)[,1:2],getss(ss.G.3,500,3)[,1:2],getss(ss.G.4,500,3)[,1:2]))
sesp[4,]<- as.numeric(c(getss(ss.pcr.1,500,-39)[,1:2],getss(ss.pcr.2,500,-39)[,1:2],getss(ss.pcr.3,500,-39)[,1:2],getss(ss.pcr.4,500,-39)[,1:2]))

