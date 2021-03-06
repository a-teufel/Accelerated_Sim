library(ggplot2)
library(cowplot)
library(plyr)
library(reshape2)

#draws distribution of selective coefficents for protein simulations


#read in data
WF.files<-list.files(pattern="fix_rep*")
MH.files<-list.files(pattern="mh_rep*")
s_wf<-NULL
s_mh<-NULL

for( i in 1:3){
print(WF.files[i])
print(MH.files[i])
WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")


wf_score<-WF$Rosetta
wf_delta<-WF$Score

mh_score<-MH$Rosetta
mh_delta<-MH$Score

#calculate selection coeff
for(i in 2:500){
   wf_org<-wf_score[i]
   mh_org<-mh_score[i]

   wf_to<-wf_score[i+1]
   mh_to<-mh_score[i+1]

   #[1] is the inital minimized score
   wf_fit<-1/(exp(wf_org-wf_score[1]/2)+1)
   mh_fit<-1/(exp(mh_org-mh_score[1]/2)+1)

   print(wf_fit)

   wf_fit_to<-1/(exp(wf_to-wf_score[1]/2)+1)
   mh_fit_to<-1/(exp(mh_to-mh_score[1]/2)+1)

   s_wf<-c(s_wf, log(wf_fit_to)-log(wf_fit) )
   s_mh<-c(s_mh, log(mh_fit_to)-log(mh_fit) )
   }


}

#set up colors
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
ggcols=ggplotColours(n = 2)
cols <- c("Accelerated\nSampling"=ggcols[1],"Original\nSampling"=ggcols[2])


#put everything in a data frame
temp<-data.frame(wf=as.vector(as.numeric(s_wf)),mh=as.vector(as.numeric(s_mh)))


#there is an NA at the end because I go too far in the loop...it doesnt really matter
w<-na.exclude(temp)$wf
m<-na.exclude(temp)$mh

#figure out which ones are neutral or nearly neutral, changes with population size
neu_w<- w[w<.01 & w>-.01]
neu_m<- m[m<.01 & m>-.01]

print(length(neu_w))
print(length(neu_m))

print(length(w))
print(length(m))

#fraction of changes that are neutral or nearly neutral
print(length(neu_w)/length(w))
print(length(neu_m)/length(m))

#run a ks, t and f test
ks.test(w,m)
t.test(w,m)
var.test(w,m)


#function to estimate what the scaler on mh method should be so that the selective coefficient distributions are the same. You should change that 100 part if you are using a differnt population size
fr<-function(s){
	scale_mh<-s*m*2*100
        scale_wf<-w*2*100
        ans<-1-ks.test(scale_mh,scale_wf)$p.value
	return(ans)

}

#run opt on the function to see what it comes back with, try a few differnt values to make sure it always comes back with the same thing and doesnt get stuck. hill climbers can be kinda crap. 
optim(c(1.5),fr)

#make density plotof selective coeffients
dens_sc<-ggplot(temp, aes(mh,fill="Accelerated\nSampling"))+geom_density(alpha=.5)+geom_density(aes(wf,fill="Original\nSampling"),alpha=.5)+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.key.size = unit(2, 'lines'))+scale_fill_manual(name="",values=cols) +xlab("Selection Coefficient")+ylab("\nDensity")


#sort the data into a new data frame to draw a qq plot
temp<-data.frame(wf=sort(as.vector(as.numeric(s_wf))),mh=sort(as.vector(as.numeric(s_mh))))
qq<- ggplot(data=temp, aes(x=(wf),y=(mh))) + geom_point()+geom_abline(colour = "red", size = 1)+ylim(-.026,.026)+xlim(-.026,.026)
qq<-qq+ylab("Selection Coefficient\nAccelerated Sampling")+xlab("Selection Coefficient Original Sampling")


pdf("protein_sc_qq.pdf", onefile=FALSE)
final_plot<-plot_grid(dens_sc,q,labels=c("A","B"),ncol=1,align='v')
final_plot
dev.off()

