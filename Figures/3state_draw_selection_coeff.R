library(ggplot2)
library(cowplot)
library(plyr)
library(reshape2)
#make graph of selection coefficents from data made with 3_state.py

s_wf<-NULL
s_mh<-NULL

MH.files<-list.files(pattern="MH_3_state*")
WF.files<-list.files(pattern="WF_3_state*")


#read in the data and calculate selective coefficents
for( i in 1:length(MH.files)){
WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")

#remove burn in phase
wf_score<-tail(WF$Score,2500)
mh_score<-tail(MH$Score,2500)

#calculate selective coeff
for(i in 2:2500){
   print(wf_score[i])
   print(wf_score[i-1])

   s_wf<-c(s_wf, log(wf_score[i])-log(wf_score[i-1]))
   s_mh<-c(s_mh, log(mh_score[i])-log(mh_score[i-1]))
   }

}

#run a ks, a t test and f test
ks.test(s_wf,s_mh)
t.test(s_wf,s_mh)
var.test(s_wf,s_mh)

#choose colors for ggplot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggcols=ggplotColours(n = 2)
cols <- c("Accelerated\nSampling"=ggcols[1],"Original\nSampling"=ggcols[2])


#put data into data frame to make ggplot happy
temp<-data.frame(wf=as.vector(as.numeric(s_wf)),mh=as.vector(as.numeric(s_mh)))


sc<-ggplot(temp, aes(mh,fill="Accelerated\nSampling"))+geom_density(alpha=.5)+geom_density(aes(wf,fill="Original\nSampling"),alpha=.5)+xlim(-.28,.28)+ theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.key.size = unit(2, 'lines'))+scale_fill_manual(name="",values=cols) +xlab("Selection Coefficient")+ylab("Density")

#use save_plot to make the aspect_ratio look reasonable
save_plot("Selective_Coeff_dist_3_state.pdf",sc,base_aspect_ratio=1.1)
