library(ggplot2)
library(cowplot)
require(gridExtra)
library(gridGraphics)
library("qgraph")
library("igraph")


layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(oma=c(2,2,1,1)) 
par(mar=c(5.1,4.1,2.1,2.1) ) 

# An adjacency matrix:
A <- matrix(1,3,3)
# igraph graph and layout:
Graph <- graph.adjacency(A)
Layout <- matrix(c(0,.5, -.5,-.5, .5,-.5),byrow=TRUE,nrow=3)


#barplot of landscape
barplot(c(1,1.15,1.25), type="o",ylab=c("Fitness"),yaxt="n", pch=19,xlab=c("State"),ylim=c(0,1.25),names=c("A","B","C"))
axis(side=1, at=c(0.7,1.9,3.1), labels = c("","",""))
axis(side=2, at=c(0,.25,.5,.75,1,1.25), labels = c("0.0","0.25","0.50","0.75","1.00","1.25"))
mtext("A", side=3, line=-2, adj=0.0, cex=1, outer=TRUE,font=2)  


#directed graphs
qgraph(get.adjacency(Graph,sparse=FALSE),layout=Layout,diag=TRUE,directed=TRUE,labels=c("A","B","C"),edge.color=c("black"),vsize=15,rescale=FALSE,esize=1,edge.labels=c("0.0","0.04","0.009","0.5","0.46","0.111","0.5","0.5","0.88"),edge.label.cex=1.5,edge.label.position=c(.5,.6,.6,.6,.5,.6,.6,.6,.5))
mtext("B", side=3, line=-22, adj=0.0, cex=1, outer=TRUE,font=2)  

Graph <- graph.adjacency(A)
Layout <- matrix(c(0,.5, -.5,-.5, .5,-.5),byrow=TRUE,nrow=3)

qgraph(get.adjacency(Graph,sparse=FALSE),layout=Layout,diag=TRUE,directed=TRUE,labels=c("A","B","C"),edge.color=c("black"),vsize=15,rescale=FALSE,esize=1,edge.labels=c("0.688","0.01","0.003","0.131","0.895","0.021","0.182","0.095","0.976"),edge.label.cex=1.5,edge.label.position=c(.5,.6,.6,.6,.5,.6,.6,.6,.5))
mtext("C", side=3, line=-22, adj=0.55, cex=1, outer=TRUE,font=2)  






