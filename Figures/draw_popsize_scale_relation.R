library(ggplot2)
library(cowplot)
library(plyr)
library(reshape2)
#draws relationship between modified population size for the accelrated method and the true population, the conversion factors came from running the draw_selection_qq_protein.R on a bunch of runs with different population sizes


x<-c(100,50,25)
y<-c(100/1.5375,50/1.35,25/.99)

temp<-data.frame(x=x,y=y)

last<-ggplot(data=temp, aes(x=x, y=y)) +geom_point()+xlab(expression(paste(N[e]," Original Sampling")))+ylab(expression(paste("Rescaled ", N[e]," Accelerated Sampling")))+xlim(0,100)+ylim(0,100)+geom_smooth(method="lm",se=FALSE,colour="red",fullrange=TRUE)

#fits linear model to data
smooth_vals = lm(y~x)
print(smooth_vals)

save_plot("convertN.pdf",last,base_aspect_ratio=1.1)
