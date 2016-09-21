library(ggplot2)
library(cowplot)
library(gam)
require(mgcv)
#draws probability of accepting reversion mutation over next 15 steps


#color scale to match ggplots default
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
ggcols=ggplotColours(n = 2)
cbbPalette <- c(ggcols[1],ggcols[2])
cols <- c("Accelerated\nSampling"=cbbPalette[1],"Original\nSampling"=cbbPalette[2])

#read in data for MH reversion probs
MH.files<-list.files(pattern="rev_15_score_rep*")
by_15<-NULL

for( i in 1:length(MH.files)){
  MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")
  prob<-MH$Prob
  mat<-matrix(prob,nrow=15)
  by_15<-cbind(by_15,mat)
}

#normalize MH probs
by_15<-by_15 / t(replicate(nrow(by_15), colSums(by_15)))

mh_r_m<-rowMeans(by_15,na.rm = TRUE)

#read in data for WF reversion probs
WF.files<-list.files(pattern="rev_15_score_fix_rep*")
by_15_wf<-NULL

for( i in 1:length(WF.files)){
  WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
  prob<-WF$Prob
  mat<-matrix(prob,nrow=15)
  by_15_wf<-cbind(by_15_wf,mat)
}

#normalize wf probs
by_15_wf<-by_15_wf / t(replicate(nrow(by_15_wf), colSums(by_15_wf)))

wf_r_m<-rowMeans(by_15_wf)

#set up data frame mh and wf reversions
dat<-data.frame(revert=1:15,mh=as.vector(as.numeric(by_15)),wf=as.vector(as.numeric(by_15_wf)),wf_mean=as.vector(as.numeric(wf_r_m)),mh_mean=as.vector(as.numeric(mh_r_m)))


lp<-ggplot(data=dat, aes(x=revert, y=mh))
lp<-lp+geom_point(aes(x=revert,y=wf_mean,colour="Original\nSampling"),dat)
lp<-lp+geom_point(aes(x=revert,y=mh_mean,colour="Accelerated\nSampling"),dat)
lp<-lp+geom_smooth(aes(x=revert,y=mh,colour="Accelerated\nSampling"),dat,fill=cbbPalette[1],method=gam,formula=y~I(log(x))) 
lp<-lp+geom_smooth(aes(x=revert,y=wf,colour="Original\nSampling"),dat,fill=cbbPalette[2],method=gam,formula=y~I(log(x))) 
lp<-lp+labs(x="Markov Step",y="Normalized Probability of Accepting Reversion")
lp<-lp+guides(color=guide_legend(override.aes=list(fill=NA)))
lp<-lp+scale_x_continuous(breaks=seq(1,15))+theme(legend.key.size = unit(2, 'lines'), legend.key = element_rect(fill = "transparent", colour = "transparent"))+scale_colour_manual(name="",values=cols) +theme(legend.justification = c(1, 1), legend.position = c(1, 1))

pdf("stokes.pdf", onefile=FALSE)
lp
dev.off()


