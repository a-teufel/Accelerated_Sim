library(shape)
library(cowplot)
library(grid)
require(gridExtra)
rm(list=ls())

f<-function(b, thresh,x){
	return( 1/(exp(b*(x-thresh))+1))
}



a<-function(i,j,N,beta){

thresh<--17

xi<-f(beta,thresh,i)
xj<-f(beta,thresh,j)

if(xi==xj){
	return(1/N)
}


top<-(xi/xj)^2
bottom<-(xi/xj)^(2*N)

return((1-top)/(1-bottom))

}


fit1<-seq(from=-50,to=-0,by=.1)
fit2<--90

#differnt beta terms
row_N1_B1<-NULL
for(i in 1:length(fit1)){
  row_N1_B1<-rbind(row_N1_B1,a(fit2,fit1[i],4,.15))
  
}


row_N1_B5<-NULL
for(i in 1:length(fit1)){

  row_N1_B5<-rbind(row_N1_B5,a(fit2,fit1[i],4,1))
}

row_N1_B2<-NULL
for(i in 1:length(fit1)){


  row_N1_B2<-rbind(row_N1_B2,a(fit2,fit1[i],4,.07))
}


te<-expression(symbol('\256'))
nam<-paste(paste("delta"),"G of Mutant")

temp<-data.frame(x=fit1,y=row_N1_B1)

p <- ggplot(temp,aes(x=fit1,y=row_N1_B1))+geom_line()+xlim(c(fit1[1],fit1[length(fit1)]))
p<-p+geom_line(aes(x=fit1,y=row_N1_B5),col=c("red"))
p<-p+geom_line(aes(x=fit1,y=row_N1_B2),col=c("red"))+
geom_segment(aes(x = fit1[length(fit1)/2]+3.2, y = .1413, xend =  fit1[length(fit1)/2]+3.2, yend =.1285), arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+
geom_segment(aes(x = fit1[length(fit1)/2]-6.45, y = .1285, xend =  fit1[length(fit1)/2]-6.45, yend =.1413), arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+

geom_segment(aes(x = fit1[length(fit1)/2]-1.35, y = .125, xend = fit1[length(fit1)/2]+6 , yend = .125),size=.05, arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+
geom_segment(aes(x = fit1[length(fit1)/2]-2.5, y = .125, xend = fit1[length(fit1)/2]-12 , yend = .125),size=.05, arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+

annotate("text",x=fit1[length(fit1)/2]+2.38,y=.135,label=paste("beta"),parse=TRUE,col=c("black"))+annotate("text",x=fit1[length(fit1)/2]-7.25,y=.135,label=paste("beta"),parse=TRUE,col=c("black"))+
xlab(label=expression(paste(Delta,"G of Mutant")))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
scale_y_continuous(name="Probability of\nAccepting Mutation",limit=c(0,.25))

base <-theme(panel.border=element_blank(), axis.line=element_line())

p<-p+base


#data for differnt Ns
row_N10_B1<-NULL
for(i in 1:length(fit1)){

  row_N10_B1<-rbind(row_N10_B1,a(fit2,fit1[i],10,.15))
}


row_N100_B1<-NULL
for(i in 1:length(fit1)){

  row_N100_B1<-rbind(row_N100_B1,a(fit2,fit1[i],3,.15))
}


temp<-data.frame(x=fit1,y=row_N1_B1)
p2 <- ggplot(temp,aes(x=fit1,y=row_N1_B1))+geom_line()
p2<-p2+geom_line(aes(x=fit1,y=row_N10_B1),col=c("blue"))
p2<-p2+geom_line(aes(x=fit1,y=row_N100_B1),col=c("blue"))+
geom_segment(aes(x = fit1[length(fit1)/2]-3.45, y = .095, xend =  fit1[length(fit1)/2]-3.45, yend =.1078), arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+
geom_segment(aes(x = fit1[length(fit1)/2]-3.45, y = .2212, xend =  fit1[length(fit1)/2]-3.45, yend =.2082), arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+

geom_segment(aes(x = fit1[length(fit1)/2]-5.75, y = .175, xend = fit1[length(fit1)/2]-5.75 , yend = .25),size=.05, arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+
geom_segment(aes(x = fit1[length(fit1)/2]-5.75, y = .165, xend = fit1[length(fit1)/2]-5.75 , yend = .035),size=.05, arrow = arrow(length = unit(0.1, "cm"),type = "closed"))+

annotate("text",x=fit1[length(fit1)/2]-4.69,y=.2125,label=paste("N[e]"),parse=TRUE,col=c("black"))+
annotate("text",x=fit1[length(fit1)/2]-4.69,y=.1,label=paste("N[e]"),parse=TRUE,col=c("black"))+
xlab(label=expression(paste(Delta,"G of Mutant")))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
scale_y_continuous(name="Probability of\nAccepting Mutation",limit=c(0,.35),breaks=c(0,.05,.10,.15,.20,.25,.30,.35))
p2<-p2+base


t<-list()
t[[1]]<-p
t[[2]]<-p2

s<-plot_grid(plotlist=t,labels=c("A","B"),nrow=2)
s

