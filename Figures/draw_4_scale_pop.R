library(ggplot2)
library(cowplot)
library(gridExtra)
library(gam)
require(mgcv)

#draws the plot for the scaled population size, when wf pop size is 50 and MH pop is 37


#color scale to match ggplots default
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggcols=ggplotColours(n = 2)
cbbPalette <- c(ggcols[1],ggcols[2])

base <-theme(panel.border=element_blank(), axis.line=element_line())
  base <- base + theme(axis.title.x = element_text(size=32/3, vjust=-0.25/3))
  base <- base + theme(axis.text.x = element_text(size=24/3, vjust=1.3/3))
  base <- base + theme(axis.title.y = element_text(size=32/3, vjust=-0.5/3))
  base <- base + theme(axis.text.y = element_text(size=24/3, hjust=1/3))
  base <- base + theme(axis.line = element_line(colour = 'black', size = 1.5/3))
  base <- base + theme(axis.ticks.y = element_line(colour = 'black', size = 1.5/3))
  base <- base + theme(axis.ticks.x = element_line(colour = 'black', size = 1.5/3))



data_MH <- read.csv('mh_rep1_500_N50_1.35.csv',header=TRUE)
data2_MH <- read.csv('mh_rep2_500_N50_1.35.csv',header=TRUE)
data3_MH <- read.csv('mh_rep3_500_N50_1.35.csv',header=TRUE)

data_WF <- read.csv('fix_rep1_500_N50.csv',header=TRUE)
data2_WF <- read.csv('fix_rep_2_500_N50.csv',header=TRUE)
data3_WF <- read.csv('fix_rep3_500_N50.csv',header=TRUE)


#get the data to draw the delta G, these only have the first 500 mutations in them 
gens_MH<-as.vector(data_MH$Rosetta)
gens2_MH<-as.vector(data2_MH$Rosetta)
gens3_MH<-as.vector(data3_MH$Rosetta)

gens_WF<-as.vector(data_WF$Rosetta)
gens2_WF<-as.vector(data2_WF$Rosetta)
gens3_WF<-as.vector(data3_WF$Rosetta)

#put everything into a matrix to calc means
mat_MH<-cbind(as.vector(as.numeric(gens_MH)),as.vector(as.numeric(gens2_MH)),as.vector(as.numeric(gens3_MH)))

mat_WF<-cbind(as.vector(as.numeric(gens_WF)),as.vector(as.numeric(gens2_WF)),as.vector(as.numeric(gens3_WF)))

mh_mean<-rowMeans(mat_MH)

wf_mean<-rowMeans(mat_WF)


nam<-paste(paste("delta"),"G of Mutant")
cols <- c("Accelerated\nSampling"=cbbPalette[1],"Original\nSampling"=cbbPalette[2])

temp_data<-data.frame(seq=0:500,MH=as.vector(as.numeric(mat_MH)),WF=as.vector(as.numeric(mat_WF)),MH_mean=mh_mean,WF_mean=wf_mean)

mh_all<-ggplot(temp_data,aes(x=seq,y=MH))
mh_all<-mh_all+geom_smooth(aes(x=seq,y=WF,colour= "Original\nSampling"),temp_data,method=loess,span=.5,fill=cbbPalette[2],fullrange=FALSE)
mh_all<-mh_all+geom_point(aes(x=seq,y=WF),col=cbbPalette[2],size=.2,alpha=.25)
mh_all<-mh_all+geom_point(aes(x=seq,y=MH,colour="Accelerated\nSampling"),col=cbbPalette[1],size=.2,alpha=.25)
mh_all<-mh_all+geom_smooth(aes(x=seq,y=MH,colour="Accelerated\nSampling"),fill=cbbPalette[1],method=loess,span=.5,fullrange=FALSE)
mh_all<-mh_all+ylab(expression(paste(Delta,G)))
mh_all<-mh_all+xlab("Accepted Mutations")
mh_all<-mh_all+ coord_cartesian(ylim=c(-400,-250))
mh_all<-mh_all+scale_colour_manual(values=c(cbbPalette[1],cbbPalette[2]))
mh_all<-mh_all+ guides(color=guide_legend(override.aes=list(fill=NA)))
mh_all<-mh_all+theme(legend.position="none")


#get delta delta G data, and do the same thing as before
gens_MH<-as.vector(data_MH$Score)
gens2_MH<-as.vector(data2_MH$Score)
gens3_MH<-as.vector(data3_MH$Score)


gens_WF<-as.vector(data_WF$Score)
gens2_WF<-as.vector(data2_WF$Score)
gens3_WF<-as.vector(data3_WF$Score)


mat_MH<-cbind(as.vector(as.numeric(gens_MH)),as.vector(as.numeric(gens2_MH)),as.vector(as.numeric(gens3_MH)))
mat_WF<-cbind(as.vector(as.numeric(gens_WF)),as.vector(as.numeric(gens2_WF)),as.vector(as.numeric(gens3_WF)))


mh_mean<-rowMeans(mat_MH)
wf_mean<-rowMeans(mat_WF)


temp_data<-data.frame(seq=0:500,MH=as.vector(as.numeric(mat_MH)),WF=as.vector(as.numeric(mat_WF)),MH_mean=mh_mean,WF_mean=wf_mean)

mh_all2<-ggplot(data=temp_data,aes(x=seq,y=MH,colour="Accelerated\nSampling"))
mh_all2<-mh_all2 +ylab(expression(paste(paste(Delta,Delta),G)))+xlab("Accepted Mutations")
mh_all2<-mh_all2+ geom_point(aes(x=seq,y=WF),col=cbbPalette[2],temp_data,size=.2,alpha=.25)
mh_all2<-mh_all2+ geom_point(aes(x=seq,y=MH),col=cbbPalette[1],temp_data,size=.2,alpha=.25)
mh_all2<-mh_all2+ geom_smooth(aes(x=seq,y=WF,colour="Original\nSampling"),temp_data,fill=cbbPalette[2])
mh_all2<-mh_all2+ geom_smooth(aes(x=seq,y=MH,colour="Accelerated\nSampling"),temp_data,fill=cbbPalette[1]) 
mh_all2<-mh_all2+ theme(legend.position="none")
mh_all2<-mh_all2+ coord_cartesian(ylim=c(-10,10))
mh_all2<-mh_all2+scale_colour_manual(name="",values=cols) 


#read in a different set of files to draw the selective coefficent distributions
WF.files<-c('fix_rep1_50.csv','fix_rep2_50.csv','fix_rep3_50.csv')
MH.files<-c('mh_rep1_37.csv','mh_rep2_37.csv','mh_rep_3_37.csv')

s_wf<-NULL
s_mh<-NULL

for( i in 1:3){
WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")


wf_score<-WF$Rosetta
wf_delta<-WF$Score

mh_score<-MH$Rosetta
mh_delta<-MH$Score

#calc the selective coefficents
for(i in 2:500){
   wf_org<-wf_score[i]
   mh_org<-mh_score[i]

   wf_to<-wf_score[i+1]
   mh_to<-mh_score[i+1]

   wf_fit<-1/(exp(wf_org-wf_score[1]/2)+1)
   mh_fit<-1/(exp(mh_org-mh_score[1]/2)+1)

 
   wf_fit_to<-1/(exp(wf_to-wf_score[1]/2)+1)
   mh_fit_to<-1/(exp(mh_to-mh_score[1]/2)+1)

   s_wf<-c(s_wf, log(wf_fit_to)-log(wf_fit) )
   s_mh<-c(s_mh, log(mh_fit_to)-log(mh_fit) )
}

}

#put selective coeffs into a data frame
temp<-data.frame(wf=as.vector(as.numeric(s_wf)),mh=as.vector(as.numeric(s_mh)))

selec<-ggplot(temp, aes(mh,fill="Accelerated\nSampling"))+geom_density(alpha=.5)+geom_density(aes(wf,fill="Original\nSampling"),alpha=.5)+scale_fill_manual(name="",values=cols)+ xlab("Selection Coefficient")+ylab("\nDensity")+ theme(legend.position="none")


#read in files with the version probabilities in them
MH.files<-c('rep_1_mh_37_rev_15_score.csv','rep_2_mh_37_rev_15_score.csv','rep_3_mh_37_rev_15_score.csv')

by_15<-NULL

for( i in 1:length(MH.files)){

MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")
prob<-MH$Prob
mat<-matrix(prob,nrow=15)

by_15<-cbind(by_15,mat)

}

#normalize the reversion prob
by_15<-by_15 / t(replicate(nrow(by_15), colSums(by_15)))

mh_r_m<-rowMeans(by_15,na.rm = TRUE)


#get the version probs for WF
WF.files<-c('rep_1_fix_50_rev_15_score.csv','rep_2_fix_50_rev_15_score.csv','rep_3_fix_50_rev_15_score.csv')

by_15_wf<-NULL

for( i in 1:length(WF.files)){


WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
prob<-WF$Prob
mat<-matrix(prob,nrow=15)

by_15_wf<-cbind(by_15_wf,mat)

}

#normalize reversion probs WF
by_15_wf<-by_15_wf / t(replicate(nrow(by_15_wf), colSums(by_15_wf)))


wf_r_m<-rowMeans(by_15_wf)

#put everything in dataframe for ggplot
dat<-data.frame(revert=1:15,mh=as.vector(as.numeric(by_15)),wf=as.vector(as.numeric(by_15_wf)),wf_mean=as.vector(as.numeric(wf_r_m)),mh_mean=as.vector(as.numeric(mh_r_m)))


rev<-ggplot(data=dat, aes(x=revert, y=mh))
rev<-rev+geom_point(aes(x=revert,y=wf_mean,colour="Original\nSampling"),dat)
rev<-rev+geom_point(aes(x=revert,y=mh_mean,colour="Accelerated\nSampling"),dat)
rev<-rev+geom_smooth(aes(x=revert,y=mh,colour="Accelerated\nSampling"),dat,fill=cbbPalette[1],method=gam,formula=y~I(log(x))) 
rev<-rev+geom_smooth(aes(x=revert,y=wf,colour="Original\nSampling"),dat,fill=cbbPalette[2],method=gam,formula=y~I(log(x))) 
rev<-rev+labs(x="Markov Step",y="Normalized Probability\nof Accepting Reversion")
rev<-rev+guides(color=guide_legend(override.aes=list(fill=NA)))
rev<-rev+scale_x_continuous(breaks=c(1,5,10,15))+theme(legend.key.size = unit(2, 'lines'), legend.key = element_rect(fill = "transparent", colour = "transparent"))+scale_colour_manual(name="",values=cols) +theme(legend.justification = c(1, 1), legend.position = c(1, 1))


pdf("all4_rescale_pop.pdf", onefile=FALSE)
plot_grid(mh_all,selec,mh_all2,rev,labels=c("A","C","B","D"),ncol=2,align='v')
dev.off()


