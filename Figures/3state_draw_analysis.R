library(ggplot2)
library(cowplot)
library(grid)
#draws figures for 3 state simulations, calculates transitions

##set of helper functions

#std function for counts
std <- function(x) sd(x)/sqrt(length(x))

#color scale to match ggplots default
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#counts number of transition types from a data file
CountTrans<-function(method){
AB=0
AC=0
BA=0
BC=0
CA=0
CB=0

#removes the burn in phase
vars<-tail(method$Variant,2500)

for(i in 1:(length(vars)-1)){
	
	if(vars[i] == "A" & vars[i+1] == "B"){
		AB<-AB+1
	}
	if(vars[i] == "A" & vars[i+1] == "C"){
		AC<-AC+1
	}
	
	if(vars[i] == "B" & vars[i+1] == "A"){
		BA<-BA+1
	}
	if(vars[i] == "B" & vars[i+1] == "C"){
		BC<-BC+1
	}
	
	if(vars[i] == "C" & vars[i+1] == "A"){
		CA<-CA+1
	}
	if(vars[i] == "C" & vars[i+1] == "B"){
		CB<-CB+1
	}

}

return(c(AB,AC,BA,BC,CA,CB))
}


##draws the images from data that 3_state.py makes
#get the full sequences
MH.files<-list.files(pattern="MH_3_state*")
WF.files<-list.files(pattern="WF_3_state*")
all_MH_cnt<-NULL
all_WF_cnt<-NULL


for( i in 1:length(MH.files)){
  MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")
  WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")
  
  #take just the end to remove burn in phase
  cnt_MH<-table(tail(MH$Variant,2500))
  all_MH_cnt<-rbind(all_MH_cnt,cnt_MH)

  cnt_WF<-table(tail(WF$Variant,2500))
  all_WF_cnt<-rbind(all_WF_cnt,cnt_WF)

}

#run t and f test
print(all_MH_cnt)
print(all_WF_cnt)
print(t.test(all_MH_cnt[,1],all_WF_cnt[,1]))
print(t.test(all_MH_cnt[,2],all_WF_cnt[,2]))
print(t.test(all_MH_cnt[,3],all_WF_cnt[,3]))

print(var.test(all_MH_cnt[,1],all_WF_cnt[,1]))
print(var.test(all_MH_cnt[,2],all_WF_cnt[,2]))
print(var.test(all_MH_cnt[,3],all_WF_cnt[,3]))


#jams data into a frame to make ggplot happy
dat1 <- data.frame(
    Method = factor(c("Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling")),
    State = factor(c("A","A","B","B","C","C"), levels=c("A","B","C")),
    total = c(mean(all_MH_cnt[,1]), mean(all_WF_cnt[,1]), mean(all_MH_cnt[,2]), mean(all_WF_cnt[,2]), mean(all_MH_cnt[,3]), mean(all_WF_cnt[,3])),
    se=c(std(all_MH_cnt[,1]), std(all_WF_cnt[,1]), std(all_MH_cnt[,2]), std(all_WF_cnt[,2]), std(all_MH_cnt[,3]), std(all_WF_cnt[,3]))
)


bar<-ggplot(data=dat1, aes(x=State, y=total, fill=Method)) +
geom_bar(stat="identity", position="dodge")+ 
geom_errorbar(aes(ymin=total-se, ymax=total+se),width=.2,           
position=position_dodge(.9))
bar<-bar+labs(y="# of Times Sampled",x="State",fill="")
bar<-bar+ theme(legend.position = "right",legend.key.size = unit(2, 'lines')) +scale_colour_manual(name="",values=cols)


#makes types of transitions plot
#read in data again, cuz I'm lazy and just cut and pasted this
MH.files<-list.files(pattern="MH_3_state*")
WF.files<-list.files(pattern="WF_3_state*")


all_MH_trans<-NULL
all_WF_trans<-NULL
for( i in 1:length(MH.files)){
  MH<-read.csv(file=MH.files[i],head=TRUE,sep=",")
  WF<-read.csv(file=WF.files[i],head=TRUE,sep=",")

  #get counts and normalize
  transMH<-CountTrans(MH)
  all_MH_trans<-rbind(all_MH_trans,transMH/sum(transMH))

  transWF<-CountTrans(WF)
  all_WF_trans<-rbind(all_WF_trans,transWF/sum(transWF))

}


print(all_MH_trans)
print(all_WF_trans)

#check to see if the means of the normalized counts are different
print(t.test(all_MH_trans[,1],all_WF_trans[,1]))
print(t.test(all_MH_trans[,2],all_WF_trans[,2]))
print(t.test(all_MH_trans[,3],all_WF_trans[,3]))
print(t.test(all_MH_trans[,4],all_WF_trans[,4]))
print(t.test(all_MH_trans[,5],all_WF_trans[,5]))
print(t.test(all_MH_trans[,6],all_WF_trans[,6]))



dat2 <- data.frame(
    Method = factor(c("Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling","Accelerated\nSampling","Original\nSampling")),
    State = factor(c("A->B","A->B","A->C","A->C","B->A","B->A","B->C","B->C","C->A","C->A","C->B","C->B"), levels=c("A->B","A->C","B->A","B->C","C->A","C->B")),
    total = c(mean(all_MH_trans[,1]), mean(all_WF_trans[,1]), mean(all_MH_trans[,2]), mean(all_WF_trans[,2]), mean(all_MH_trans[,3]), 
mean(all_WF_trans[,3]), mean(all_MH_trans[,4]), mean(all_WF_trans[,4]), mean(all_MH_trans[,5]), mean(all_WF_trans[,5]), mean(all_MH_trans[,6]), 
mean(all_WF_trans[,6])),
se=c(std(all_MH_trans[,1]), std(all_WF_trans[,1]), std(all_MH_trans[,2]), std(all_WF_trans[,2]), std(all_MH_trans[,3]), std(all_WF_trans[,3]),
std(all_MH_trans[,4]), std(all_WF_trans[,4]), std(all_MH_trans[,5]), std(all_WF_trans[,5]), std(all_MH_trans[,6]), std(all_WF_trans[,6])
)
)

labs=expression(A%->%B,A%->%C,B%->%A,B%->%C,C%->%A,C%->%B)

# figure out where to draw the lines on the sig diff bars, you have to look at the data and figure out which one is the higher one
ABhigh=mean(all_MH_trans[,1])+std(all_MH_trans[,1])
ABlow=mean(all_WF_trans[,1])+std(all_WF_trans[,1])

AChigh=mean(all_WF_trans[,2])+std(all_WF_trans[,2])
AClow=mean(all_MH_trans[,2])+std(all_MH_trans[,2])


BAhigh=mean(all_MH_trans[,3])+std(all_MH_trans[,3])
BAlow=mean(all_WF_trans[,3])+std(all_WF_trans[,3])

BChigh=mean(all_MH_trans[,4])+std(all_MH_trans[,4])
BClow=mean(all_WF_trans[,4])+std(all_WF_trans[,4])

CAhigh=mean(all_WF_trans[,5])+std(all_WF_trans[,5])
CAlow=mean(all_MH_trans[,5])+std(all_MH_trans[,5])

CBhigh=mean(all_MH_trans[,6])+std(all_MH_trans[,6])
CBlow=mean(all_WF_trans[,6])+std(all_WF_trans[,6])


bar2<-ggplot(data=dat2, aes(x=State, y=total, fill=Method)) +
geom_bar(stat="identity", position="dodge")+ 
geom_errorbar(aes(ymin=total-se, ymax=total+se),width=.2,           
position=position_dodge(.9))
bar2<-bar2+labs(y="Fraction of Transitions",x="Transition",fill="")
bar2<-bar2 + scale_x_discrete(labels = labs);
bar2<-bar2+ylim(0,0.5)
bar2<-bar2+ theme(legend.position="none", axis.text.x = element_text(size=9))
bar2<-bar2+ annotate(x=c(.75,.75,1.25,1.25),y=c(ABhigh+.01,ABhigh+.02,ABhigh+.02,ABlow+.01),"path")
bar2<-bar2+ annotate(x=c(1.75,1.75,2.25,2.25),y=c(AClow+.01,AChigh+.02,AChigh+.02,AChigh+.01),"path")
bar2<-bar2+ annotate(x=c(3.75,3.75,4.25,4.25),y=c(BChigh+.01,BChigh+.02,BChigh+.02,BClow+.01),"path")
bar2<-bar2+ annotate(x=c(4.75,4.75,5.25,5.25),y=c(CAlow+.01,CAhigh+.02,CAhigh+.02,CAhigh+.01),"path")
bar2<-bar2+ annotate(x=c(5.75,5.75,6.25,6.25),y=c(CBhigh+.01,CBhigh+.02,CBhigh+.02,CBlow+.01),"path")
bar2<-bar2+ annotate("text",x=1,y=ABhigh+.025,label="**")
bar2<-bar2+ annotate("text",x=2,y=AChigh+.025,label="***")
bar2<-bar2+ annotate("text",x=4,y=BChigh+.025,label="***")
bar2<-bar2+ annotate("text",x=5,y=CAhigh+.025,label="***")
bar2<-bar2+ annotate("text",x=6,y=CBhigh+.025,label="***")
MH.files<-list.files(pattern="MH_3_state*")
WF.files<-list.files(pattern="WF_3_state*")


#change the widths so the two images line up
g.bar <- ggplotGrob(bar) # convert to gtable
g.bar2 <- ggplotGrob(bar2) # convert to gtable

if(getRversion() < "3.3.0"){
  # We need to convert the widths to unit lists for the subsequent
  # manipulations to be possible.
  # Once R 3.3.0 is released, this will not be necessary anymore.
  g.bar$widths <- grid:::unit.list(g.bar$widths)   
  g.bar2$widths <-  grid:::unit.list(g.bar2$widths)
}

bar.widths <- g.bar$widths[1:3] # extract the first three widths, 
                                  # corresponding to left margin, y lab, and y axis
bar2.widths <- g.bar2$widths[1:3] # same for mpg plot
max.widths <- unit.pmax(bar.widths, bar2.widths) # calculate maximum widths
g.bar$widths[1:3] <- max.widths # assign max widths to iris gtable
g.bar2$widths[1:3] <- max.widths # assign max widths to mpg gtable

# plot_grid() can work directly with gtables, so this works
pdf("3_state_analysis.pdf", onefile=FALSE)
plot_grid(g.bar, g.bar2, labels = "AUTO", ncol = 1)
dev.off()


