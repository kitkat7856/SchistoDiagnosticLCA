
### Select status posteriors, sample 500
set.seed(123)
status<- select(out, starts_with("Status"))
x<- sample(40000,500, replace=F)
status_samples<- status[x,]

### Convert status to binary, so status 1 = 0 for uninfected, and status 2 and 3 = 1 for infected
status_samples_bin<- matrix(NA,500,840)
status_samples_bin[status_samples==1]<-0
status_samples_bin[status_samples==2 | status_samples==3]<-1

### Separate by timepoint
status1_bin<-status_samples_bin[,1:210] ### preT
status2_bin<-status_samples_bin[,211:420] ### 3wks
status3_bin<-status_samples_bin[,421:630] ### 9wks
status4_bin<-status_samples_bin[,631:840] ### 6mo

### Each row is a model draw, get prev estimate for each draw (500 samples)
prev_1<- rowSums(status1_bin)/210 
prev_2<- rowSums(status2_bin)/210 
prev_3<- rowSums(status3_bin)/210 
prev_4<- rowSums(status4_bin)/210 

median(prev_1)
quantile(prev_1, c(0.025, 0.975))


######## Alternative method, what I used to create entire table with dif proportions #####
status1<- status_samples[,1:210]
status2<- status_samples[,211:420]
status3<- status_samples[,421:630]
status4<- status_samples[,631:840]


##### df of proportion of each status at each tmp
## shed = proportion of infected individuals that are shedding (3/2+3)
## shed.t = proportion of population that are infected and shedding (3/210)
## total = infection prevalence (2+3/210) 
prev1<- data.frame(shed=rowSums(status1==3)/rowSums(status1==2 | status1==3), shed.t=rowSums(status1==3)/210, total=rowSums(status1==2 | status1==3)/210, tmp="PreT")
prev2<- data.frame(shed=rowSums(status2==3)/rowSums(status2==2 | status2==3), shed.t=rowSums(status2==3)/210, total=rowSums(status2==2 | status2==3)/210, tmp="3wk")
prev3<- data.frame(shed=rowSums(status3==3)/rowSums(status3==2 | status3==3), shed.t=rowSums(status3==3)/210, total=rowSums(status3==2 | status3==3)/210, tmp="9wk")
prev4<- data.frame(shed=rowSums(status4==3)/rowSums(status4==2 | status4==3), shed.t=rowSums(status4==3)/210, total=rowSums(status4==2 | status4==3)/210, tmp="6mo")
prev_df<- rbind(prev1,prev2,prev3,prev4)

#### summarise - median, 95% CI
prev_sum<- prev_df %>%
  group_by(tmp=factor(tmp, levels = c("PreT", "3wk", "9wk", "6mo"))) %>%
  summarise(across(where(is.numeric), 
                   list(med=median,
                        Qlow= ~quantile(.,0.025),
                        Qhi= ~quantile(.,0.975)
                   ),
                   .names="{col}_{fn}" 
  ))



####### ROCR Curves #######
library(ROCR)
### labels are binary (status), predictions are the measurements
### ROCR allows us to generate average curve for multiple runs 

## CT ## CT values must be put in as negative since negative correlation!!
### Timepoint 1 ###

predictions1<- rep(list(na.omit(-ct[,1])),10) ### must omit NA values for ROCR

labels1<-status1_bin[c(sample(500,10, replace=F)),which(!is.na(ct[,1]))] ### binary status draws for which there are no NA values in CT
labels1_list<-as.list(as.data.frame(t(labels1))) ### convert binary status matrix, so each draw is an element in list

pred1<- prediction(predictions1, labels1_list) ## create prediction object ###
perf1<- performance(pred1,"tpr","fpr") ## ROC performance object- true positive vs false positive

plot(perf1,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-38,-40),
     print.cutoffs.at = c(-38,-40),
     text.adj = c(-1, 1),
     lwd=2,
     main= "CT, PreT")
#perf1<- performance(pred, "sens", "spec") ### Sensitivty/specificity plots


#### Timepoint 2 ###
predictions2<- rep(list(na.omit(-ct[,2])),500) 

labels2<-status2_bin[,which(!is.na(ct[,2]))] 
labels2_list<-as.list(as.data.frame(t(labels2))) 

pred2<- prediction(predictions2, labels2_list) 
perf2<- performance(pred2,"tpr","fpr") 
plot(perf2,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-35,-38,-40,-45),
     print.cutoffs.at = c(-35,-38,-40,-45),
     text.adj = c(0, -1),
     lwd=2,
     main= "CT, wk3")

#### Timepoint 3 ###
predictions3<- rep(list(na.omit(-ct[,3])),500) 

labels3<-status3_bin[,which(!is.na(ct[,3]))] 
labels3_list<-as.list(as.data.frame(t(labels3))) 

pred3<- prediction(predictions3, labels3_list) 
perf3<- performance(pred3,"tpr","fpr") 
plot(perf3,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-38,-40),
     print.cutoffs.at = c(-38,-40),
     text.adj = c(0, -1),
     lwd=2,
     main= "CT, wk9")


#### Timepoint 4 ###
predictions4<- rep(list(na.omit(-ct[,4])),500) 

labels4<-status4_bin[,which(!is.na(ct[,4]))] 
labels4_list<-as.list(as.data.frame(t(labels4))) 

pred4<- prediction(predictions4, labels4_list) 
perf4<- performance(pred4,"tpr","fpr") 
plot(perf4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-35,-38,-40),
     print.cutoffs.at = c(-35,-38,-40),
     text.adj = c(0, -1),
     lwd=2,
     main= "CT, 6mo")

#### All Timepoint averaged ####
predictions_all<- rep(list(na.omit(c(-ct))),500)

labels_all<-status_samples_bin[,which(!is.na(c(ct)))] 
labels_all_list<-as.list(as.data.frame(t(labels_all))) 

pred_all<- prediction(predictions_all, labels_all_list) 
perf_all<- performance(pred_all,"tpr","fpr") 
plot(perf4,
     avg= "threshold",
     spread.estimate= "stddev",
     spread.scale=2,
     show.spread.at = c(-35,-38,-40),
     print.cutoffs.at = c(-35,-38,-40),
     text.adj = c(0, -1),
     lwd=2,
     main= "CT, all time")


getROC<- function(diagnostic, timepoint, binary_status){
  if(length(unique(as.vector(as.matrix(binary_status))))!=2){
    stop("Status must be binary")
  }
  
  a<- 1 + 210*(timepoint-1)
  b<- 210*timepoint
  
  l<- binary_status[,a:b]
  
  prediction<- rep(list(na.omit(diagnostic[,timepoint])),500)
  label<- l[,which(!is.na(diagnostic[,timepoint]))]
  label_list<- as.list(as.data.frame(t(label)))
  
  pred<-prediction(prediction,label_list)
  roc<- performance(pred,"tpr","fpr")
  
  return(roc)
}

