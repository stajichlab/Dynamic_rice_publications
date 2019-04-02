#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def plottrait(trait, parent):

    Rcmd='''
pdf("DrawQTLtraitALL.pdf")
par(mfrow=c(4,3))

error.bar <- function(x, y, upper,color, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x-upper,y, x+upper, y, angle=90, code=3,col=color, length=length, ...)
}

read.table("''' + trait + '''",sep="\\t",header=TRUE) -> trait
read.table("''' + parent + '''",sep="\\t") -> parent
traits <- colnames(trait)
for (i in 2:ncol(trait)){
    brk = 10
    if (traits[i] == "Unique.mPing"){
       brk = c(0,2,5,10,20,40,80,120,140,180)
    } 
    hist(trait[,i],xlim=c(min(as.numeric(na.omit(trait[,2]))),max(as.numeric(na.omit(trait[,2])))),breaks=brk,plot=FALSE) -> xh
    step=10;
    yadj=max(xh$counts)*1.4/step;
    barplot(xh$counts,axes=FALSE,border=F,ylim=c(0,max(xh$counts)*1.4),col=c("cornflowerblue"),xlab=traits[i],ylab="# of RILs") -> xx
    for (j in 1:length(xx)) { # adj can not take vector, so we use loops to add text
      text(xx[j],-yadj,labels=paste(xh$breaks[j],"-",xh$breaks[j+1],sep=""),cex=0.7,srt=65,adj=c(1,1),xpd=TRUE)
    }
    axis(1,c(0,max(xx)+0.5),line=0.2,labels=c("",""))
    axis(2,seq(0,max(xh$counts)*1.4,by=step),cex.axis=0.7)
    index <- i-1
    pix <- (max(xx)-min(xx)+1)/(max(xh$mids)-min(xh$mids)+1)
    x1 <- (parent[index,2]-min(xh$mids)+1)*pix
    x2 <- (parent[index,5]-min(xh$mids)+1)*pix
    points(x1,max(xh$counts)*1.1,col="red")
    points(x2,max(xh$counts)*1.1,col="blue") 
    
    ex <- c(x1,x2)
    ey <- c(max(xh$counts)*1.1,max(xh$counts)*1.1)
    sd <- c(parent[index,4]*pix,parent[index,7]*pix)
    ec <- c("red","blue") 
    error.bar(ex,ey,sd,ec)
    
    legend(max(xx)-5,max(xh$counts)*1.2,bty="n",lty=c(1,1),cex=0.7,c("NB","HEG4"),col=c("red","blue")) 
}
dev.off()
'''
    with open ("DrawQTLtraitALL.R",'w') as filefh:
        print >> filefh, Rcmd
    os.system('cat DrawQTLtraitALL.R | R --slave')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--trait')
    parser.add_argument('-p', '--parent')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.trait) > 0
    except:
        usage()
        sys.exit(2)
   
    plottrait(args.trait,args.parent)

if __name__ == '__main__':
    main()

