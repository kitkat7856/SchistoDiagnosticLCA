library(rjags)
library(runjags)
library(dplyr) 
library(reshape)
library(stringr)
library(lubridate)
library(ggplot2)
library(bayestestR)
library(ROCR)

source("Functions.R")

kkdata<- read.csv("KK.4tp.csv")
ccadata<- read.csv("cca.clean.csv")
stoolqpcr<- read.csv("Stool Samples qpcr data.csv") 
dbspcr<- read.csv("DBS_qpcr.csv")
hatchingdata<- read.csv("hatching_longform_bg_alldates.csv")

stoolqpcr$tmp<- sub("_.*", "",stoolqpcr$Sample.Name)
stoolqpcr$CID<- sub(".*_", "",stoolqpcr$Sample.Name)

dbspcr$tmp<- sub("_.*", "",dbspcr$Sample.Name)
dbspcr$CID<- sub(".*_", "",dbspcr$Sample.Name)

# Match IDS, so each matrix row corresponds to same id
kkIDs<- unique(kkdata$child_id)
ccaIDs<- unique(ccadata$CID)
ctIDs<- unique(stoolqpcr$CID) 
dbsIDs<- unique(dbspcr$CID)
hatchIDs<- unique(hatchingdata$id_full)
hatchIDs<- unique(as.numeric(str_sub(hatchIDs,6,9)))

CIDs<- kkIDs[kkIDs %in% ccaIDs] ### 210

#make data into matrix for jags, individual x tmp
kk<- loadKK(kkdata, CIDs)
cca<- loadCCA(ccadata, CIDs)
gscore<- loadGscore(ccadata, CIDs)
ct<-loadCT(stoolqpcr, CIDs)
dbs<-loadDBS(dbspcr, CIDs,2) #binary data, ct 38 cutoff
dbs.ct<-loadDBS(dbspcr, CIDs,1)
hatch<- loadHatch(hatchingdata,CIDs,1)
hatch.binary<- loadHatch(hatchingdata,CIDs,2)

### Run model ####
mod<-eggworm(210,4,kk,gscore,ct,hatch.binary) 
out<- as.data.frame(as.mcmc(mod))

### Check model output #### 
#check posteriors normally distributed, > convergence
plot(mod)

par(mfrow=c(2,2))
plot(density(out$`prob[1,1]`)); plot(density(out$`prob[2,1]`)); plot(density(out$`prob[3,1]`)); plot(density(out$`prob[4,1]`))
plot(density(out$`prob[1,2]`)); plot(density(out$`prob[2,2]`)); plot(density(out$`prob[3,2]`)); plot(density(out$`prob[4,2]`))
plot(density(out$`prob[1,3]`)); plot(density(out$`prob[2,3]`)); plot(density(out$`prob[3,3]`)); plot(density(out$`prob[4,3]`))

plot(density(out$rtnb));plot(density(out$sh));plot(density(out$rt))
plot(density(out$`sigma[2]`));plot(density(out$`sigma[3]`));plot(density(out$`sigma[4]`));

plot(density(out$k));plot(density(out$int));plot(density(out$negcca))
plot(density(out$a_ct)); plot(density(out$b_ct))

par(mfrow=c(1,1))


###### summary posteriors #####

parameter<- select(out, !starts_with("Status"))
parameter.df<- data.frame()
for(i in 1:length(parameter)){
  name<- names(parameter)[i]
  median<- point_estimate(parameter[,i], centrality="median")
  CIlow<- hdi(parameter[,i], ci=0.89)[2]
  CIhigh<- hdi(parameter[,i], ci=0.89)[3]
  
  a<- data.frame(Parameter=name, median=round(median$Median,3), CI_low= round(CIlow$CI_low,3), CI_high= round(CIhigh$CI_high,3))
  parameter.df<- rbind(parameter.df,a)
}
parameter.df$paste<- paste0(parameter.df$median," ","(", paste(parameter.df$CI_low,parameter.df$CI_high,sep=","), ")")

## prob of each status
prob.df<- matrix(parameter.df$median[1:12], nrow=4, ncol=3)

