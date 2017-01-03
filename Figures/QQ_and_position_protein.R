library(ggplot2)
library(cowplot)
library(plyr)
library(reshape2)
#draws a QQ plot and a position specific relation for the number of subsitutions at each site. Makes an all vs. all plots for within replicate experiments

#color scale to match ggplots default
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#correlation coeff helper function
lm_eqn <- function(df){

    eq <- substitute(~R~"="~r2, 
         list(
             r2 = format(cor(df$wf,df$mh), digits = 3)))
   
    as.character(as.expression(eq));                 
}


wf_sites<-NULL
mh_sites<-NULL

WF.files<-list.files(pattern="fix_rep*")
MH.files<-list.files(pattern="mh_rep*")

for( i in 1:3){
print(WF.files[i])
print(MH.files[i])
WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")

#get the positions an remove the amino acids associated with the change
sites<-WF$Variant
s<-gsub("[^0-9]", "", sites) 
wf_sites<-c(wf_sites,s)


sitesmh<-MH$Variant
smh<-gsub("[^0-9]", "", sitesmh) 
mh_sites<-c(mh_sites,smh)

}


#put sites into a numeric vector
sit<-as.vector(as.numeric(wf_sites))
sit2<-as.vector(as.numeric(mh_sites))


#there are 300 amino acids in the protein, so break it into factors so we can use tabulate. Tabulate is needed because some sites dont have any changes
sit<-factor(sit,levels=c(1:300))
sit2<-factor(sit2,levels=c(1:300))

d<-tabulate(as.numeric(sit))
d2<-tabulate(as.numeric(sit2))


#sort the data and put it in a data frame for the qq plot
temp<-data.frame(wf=sort(as.vector(as.numeric(d))),mh=sort(as.vector(as.numeric(d2))))

qq<- ggplot(data=temp, aes(x=(wf),y=(mh))) + geom_point()+geom_abline(colour = "red", size = 1)+xlim(1,20)+ylim(1,20)
qq<-qq+ylab("# of Subsitutions\nwith Accelerated Sampling")+xlab("# of Subsitutions with Original Sampling")
qq

#unsorted into a data frame for the position based number of substitutions
temp<-data.frame(wf=d,mh=d2)
pos<-ggplot(temp,aes(wf,mh))+geom_point()+ ylim(0, 20)+ xlim(0,20)+xlab("# of Substitutions per Site with Original Sampling")+ylab("# of Substitutions per Site\nwith Accelerated Sampling")
pos<-pos +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
pos<-pos + annotate("text",x = 15, y = 19, label = lm_eqn(temp), fontface=1,parse = TRUE,size=5)
pos

pdf("QQ_site.pdf", onefile=FALSE)
plot_grid(qq,g,labels=c("A","B"),ncol=1)
dev.off()

#all vs all replicates for position based substitutions, yes the naming convention is not great. 
wf_sites<-NULL
mh_sites<-NULL

WF_1<-read.csv(file=WF.files[1],head=TRUE,sep=",")
WF_2<-read.csv(file=WF.files[2],head=TRUE,sep=",")
WF_3<-read.csv(file=WF.files[3],head=TRUE,sep=",")

MH_1<-read.csv(file=MH.files[1],head=TRUE,sep=",")
MH_2<-read.csv(file=MH.files[2],head=TRUE,sep=",")
MH_3<-read.csv(file=MH.files[3],head=TRUE,sep=",")


wf_sites_1<-NULL
wf_sites_2<-NULL
wf_sites_3<-NULL

mh_sites_1<-NULL
mh_sites_2<-NULL
mh_sites_3<-NULL

sites_1<-WF_1$Variant
s_1<-gsub("[^0-9]", "", sites_1) 
wf_sites_1<-c(wf_sites_1,s_1)

sitesmh_1<-MH_1$Variant
smh_1<-gsub("[^0-9]", "", sitesmh_1) 
mh_sites_1<-c(mh_sites_1,smh_1)

sit_1<-as.vector(as.numeric(wf_sites_1))
sit2_1<-as.vector(as.numeric(mh_sites_1))

sit_1<-factor(sit_1,levels=c(1:300))
sit2_1<-factor(sit2_1,levels=c(1:300))

d_1<-tabulate(as.numeric(sit_1))
d2_1<-tabulate(as.numeric(sit2_1))


sites_2<-WF_2$Variant
s_2<-gsub("[^0-9]", "", sites_2) 
wf_sites_2<-c(wf_sites_2,s_2)

sitesmh_2<-MH_2$Variant
smh_2<-gsub("[^0-9]", "", sitesmh_2) 
mh_sites_2<-c(mh_sites_2,smh_2)

sit_2<-as.vector(as.numeric(wf_sites_2))
sit2_2<-as.vector(as.numeric(mh_sites_2))

sit_2<-factor(sit_2,levels=c(1:300))
sit2_2<-factor(sit2_2,levels=c(1:300))

d_2<-tabulate(as.numeric(sit_2))
d2_2<-tabulate(as.numeric(sit2_2))


sites_3<-WF_3$Variant
s_3<-gsub("[^0-9]", "", sites_3) 
wf_sites_3<-c(wf_sites_3,s_3)

sitesmh_3<-MH_3$Variant
smh_3<-gsub("[^0-9]", "", sitesmh_3) 
mh_sites_3<-c(mh_sites_3,smh_3)

sit_3<-as.vector(as.numeric(wf_sites_3))
sit2_3<-as.vector(as.numeric(mh_sites_3))

sit_3<-factor(sit_3,levels=c(1:300))
sit2_3<-factor(sit2_3,levels=c(1:300))

d_3<-tabulate(as.numeric(sit_3))
d2_3<-tabulate(as.numeric(sit2_3))


temp_wf_1<-data.frame(wf=d_1,mh=d_2)
g1<-ggplot(temp_wf_1,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub. Orig. Rep 1")+ylab("# Sub. Orig. Rep 2")
g1<-g1 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g1<-g1 + annotate("text",x = 4, y = 10, label = lm_eqn(temp_wf_1), fontface=1,parse = TRUE,size=5)
g1

temp_wf_2<-data.frame(wf=d_1,mh=d_3)
g2<-ggplot(temp_wf_2,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub. Orig. Rep 1")+ylab("# Sub. Orig. Rep 3")
g2<-g2 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g2<-g2 + annotate("text",x = 4, y = 10, label = lm_eqn(temp_wf_2), fontface=1,parse = TRUE,size=5)
g2

temp_wf_3<-data.frame(wf=d_2,mh=d_3)
g3<-ggplot(temp_wf_3,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub. Orig. Rep 2")+ylab("# Sub. Orig. Rep 3")
g3<-g3 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g3<-g3 + annotate("text",x = 4, y = 10, label = lm_eqn(temp_wf_3), fontface=1,parse = TRUE,size=5)
g3


temp_mh_1<-data.frame(wf=d2_1,mh=d2_2)
g4<-ggplot(temp_mh_1,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub. Accel. Rep 1")+ylab("# Sub. Accel. Rep 2")
g4<-g4 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g4<-g4 + annotate("text",x = 4, y = 10, label = lm_eqn(temp_mh_1), fontface=1,parse = TRUE,size=5)
g4


temp_mh_2<-data.frame(wf=d2_1,mh=d2_3)
g5<-ggplot(temp_mh_2,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub. Accel. Rep 1")+ylab("# Sub. Accel. Rep 3")
g5<-g5 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g5<-g5 + annotate("text",x = 4, y = 10, label = lm_eqn(temp_mh_2), fontface=1,parse = TRUE,size=5)
g5


temp_mh_3<-data.frame(wf=d2_2,mh=d2_3)
g6<-ggplot(temp_mh_3,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub Accel. Rep 2")+ylab("# Sub. Accel. Rep 3")
g6<-g6 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g6<-g6 + annotate("text",x = 4, y = 10,  fontface=1,label = lm_eqn(temp_mh_3),parse = TRUE,size=5)
g6



temp_mh_vs<-data.frame(wf=d_1,mh=d2_1)
g7<-ggplot(temp_mh_vs,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub.\nOrig. Rep 1")+ylab("# Sub.\nAccel. Rep 1")
g7<-g7 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g7<-g7 + annotate("text",x = 4, y = 10,  fontface=1,label = lm_eqn(temp_mh_vs),parse = TRUE,size=3)
g7


temp_mh_vs_2<-data.frame(wf=d_2,mh=d2_2)
g8<-ggplot(temp_mh_vs_2,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub.\nOrig. Rep 2")+ylab("# Sub.\nAccel. Rep 2")
g8<-g8 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g8<-g8 + annotate("text",x = 4, y = 10,  fontface=1,label = lm_eqn(temp_mh_vs_2),parse = TRUE,size=3)
g8


temp_mh_vs_3<-data.frame(wf=d_3,mh=d2_3)
g9<-ggplot(temp_mh_vs_3,aes(wf,mh))+geom_point()+ ylim(0, 10)+ xlim(0,10)+xlab("# Sub.\nOrig. Rep 3")+ylab("# Sub.\nAccel. Rep 3")
g9<-g9 +geom_smooth(method="lm",formula = y ~ x,se=FALSE,colour="red")
g9<-g9 + annotate("text",x = 4, y = 10,  fontface=1,label = lm_eqn(temp_mh_vs_3),parse = TRUE,size=3)
g9


plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,labels=c("A","B","C","D","E","F","G","H","I"))

pdf("by_site_all_V_all.pdf", onefile=FALSE)
plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,labels=c("A","B","C","D","E","F","G","H","I"),align='h')
dev.off()





