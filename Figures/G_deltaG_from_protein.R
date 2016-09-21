library(ggplot2)
library(cowplot)
library(gridExtra)
#draws G and delta G for protein simulations. This program assumes that you are reading files that are only 500 long


#color scale to match ggplots default
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggcols=ggplotColours(n = 2)
cbbPalette <- c(ggcols[1],ggcols[2])

#make a base theme for ggplot
base <-theme(panel.border=element_blank(), axis.line=element_line())
  base <- base + theme(axis.title.x = element_text(size=32/3, vjust=-0.25/3))
  base <- base + theme(axis.text.x = element_text(size=24/3, vjust=1.3/3))
  base <- base + theme(axis.title.y = element_text(size=32/3, vjust=-0.5/3))
  base <- base + theme(axis.text.y = element_text(size=24/3, hjust=1/3))
  base <- base + theme(axis.line = element_line(colour = 'black', size = 1.5/3))
  base <- base + theme(axis.ticks.y = element_line(colour = 'black', size = 1.5/3))
  base <- base + theme(axis.ticks.x = element_line(colour = 'black', size = 1.5/3))

#labels and colors
nam<-paste(paste("delta"),"G of Mutant")
cols <- c("Accelerated\nSampling"=cbbPalette[1],"Original\nSampling"=cbbPalette[2])


data_MH <- read.csv('MH_run1.csv',header=TRUE)
data2_MH <- read.csv('MH_run2.csv',header=TRUE)
data3_MH <- read.csv('MH_run3.csv',header=TRUE)

data_WF <- read.csv('WF_run1.csv',header=TRUE)
data2_WF <- read.csv('WF_run2.csv',header=TRUE)
data3_WF <- read.csv('WF_run3.csv',header=TRUE)


gens_MH<-as.vector(data_MH$Rosetta)
gens2_MH<-as.vector(data2_MH$Rosetta)
gens3_MH<-as.vector(data3_MH$Rosetta)


gens_WF<-as.vector(data_WF$Rosetta)
gens2_WF<-as.vector(data2_WF$Rosetta)
gens3_WF<-as.vector(data3_WF$Rosetta)


#put everything into a big matrix so we can calc means
mat_MH<-cbind(as.vector(as.numeric(gens_MH)),as.vector(as.numeric(gens2_MH)),as.vector(as.numeric(gens3_MH)))
mat_WF<-cbind(as.vector(as.numeric(gens_WF)),as.vector(as.numeric(gens2_WF)),as.vector(as.numeric(gens3_WF)))


mh_mean<-rowMeans(mat_MH)

wf_mean<-rowMeans(mat_WF)

#t test on row means
t.test(mh_mean,wf_mean)


#collect data from delta G plot
temp_data<-data.frame(seq=0:500,MH=as.vector(as.numeric(mat_MH)),WF=as.vector(as.numeric(mat_WF)),MH_mean=mh_mean,WF_mean=wf_mean)

mh_all<-ggplot(temp_data,aes(x=seq,y=MH))
mh_all<-mh_all+geom_smooth(aes(x=seq,y=WF,colour= "Original\nSampling"),temp_data,method=loess,span=.5,fill=cbbPalette[2],fullrange=FALSE)
mh_all<-mh_all+geom_point(aes(x=seq,y=WF),col=cbbPalette[2],size=.2,alpha=.25)
mh_all<-mh_all+geom_point(aes(x=seq,y=MH,colour="Accelerated\nSampling"),col=cbbPalette[1],size=.2,alpha=.25)
mh_all<-mh_all+geom_smooth(aes(x=seq,y=MH,colour="Accelerated\nSampling"),fill=cbbPalette[1],method=loess,span=.5,fullrange=FALSE)
mh_all<-mh_all+ylab(expression(paste(Delta,G)))
mh_all<-mh_all+xlab("Accepted Mutations")
mh_all<-mh_all+coord_cartesian(ylim=c(-400,-250))
mh_all<-mh_all+scale_colour_manual(values=c(cbbPalette[1],cbbPalette[2]))
mh_all<-mh_all+theme(legend.justification = c(1, 1), legend.position = c(1, .5))
mh_all<-mh_all+ guides(color=guide_legend(override.aes=list(fill=NA)))
mh_all<-mh_all+theme(legend.key.size = unit(2, 'lines'), legend.key = element_rect(fill = "transparent", colour = "transparent"))+scale_colour_manual(name="",values=cols) 

#repeat the process for delta delta G
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

#put delta delta G into data frame
temp_data<-data.frame(seq=0:500,MH=as.vector(as.numeric(mat_MH)),WF=as.vector(as.numeric(mat_WF)),MH_mean=mh_mean,WF_mean=wf_mean)

mh_all2<-ggplot(data=temp_data,aes(x=seq,y=MH,colour="Accelerated\nSampling"))
mh_all2<-mh_all2 +ylab(expression(paste(paste(Delta,Delta),G)))+xlab("Accepted Mutations")
mh_all2<-mh_all2+ geom_point(aes(x=seq,y=WF),col=cbbPalette[2],temp_data,size=.2,alpha=.25)
mh_all2<-mh_all2+ geom_point(aes(x=seq,y=MH),col=cbbPalette[1],temp_data,size=.2,alpha=.25)
mh_all2<-mh_all2+ geom_smooth(aes(x=seq,y=WF,colour="Original\nSampling"),temp_data,fill=cbbPalette[2])
mh_all2<-mh_all2+geom_smooth(aes(x=seq,y=MH,colour="Accelerated\nSampling"),temp_data,fill=cbbPalette[1]) 
mh_all2<-mh_all2+theme(legend.position="none")
mh_all2<-mh_all2+ coord_cartesian(ylim=c(-10,10))

pdf("Engergies_first500accepts.pdf", onefile=FALSE)
plot_grid(mh_all,mh_all2,labels=c("A","B"),ncol=1,align='v')
dev.off()


