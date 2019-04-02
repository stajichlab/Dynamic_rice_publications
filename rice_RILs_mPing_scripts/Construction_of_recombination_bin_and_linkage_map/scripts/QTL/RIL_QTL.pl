#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"csv:s","qtlcart:s","project:s","help");


my $help=<<USAGE;
perl $0 --csv test.csv 
perl $0 --qtlcart test
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

if ($opt{qtlcart}){
   $opt{project} ||= $opt{qtlcart};
   Rqtl($opt{qtlcart},"qtlcart");
}else{
   $opt{project} ||= $1.".Rqtl" if ($opt{csv}=~/(.*)\.csv/);
   Rqtl($opt{csv},"csv");
}


sub Rqtl
{
my ($qtl,$type)=@_;
my $cross="cross";
my $tranperant=90;
my $flag= $type eq "qtlcart" ? 1 : 0;
my $cmd =<<R;
library(qtl)
pdf ("$opt{project}.QTL.pdf")
nph <- 24 ## number of phenotype
## step0. read and write the data of cross
if ($flag){
   read.cross("qtlcart",dir="./",file="$qtl.cro",mapfile="$qtl.map") -> $cross
}else{
   read.cross("csv",dir="./",file="$qtl") -> $cross
}

myMaq <- pull.map($cross)
plot.map(myMaq)

if (0){
write.cross($cross,"qtlcart","$opt{project}.originalmap")
convert2riself($cross) -> $cross.riself
write.cross($cross.riself,"qtlcart","$opt{project}.riself")
est.map($cross) -> reestimate.map
replace.map($cross,reestimate.map) -> $cross.reestimate
write.cross($cross.reestimate,"qtlcart","$opt{project}.reestimatemap")
}

## step1. do single QTL analysis
## 1.1 marker regression
QTL.mr <- scanone($cross,pheno.col=1:nph,method="mr")

## 1.2 do interval method stardand invertal mapping(em), Haley-Knott regression (hk), Extended HK regression (ehk), and Multiple imputation (imp)
$cross <- calc.genoprob($cross,step=1,error.prob=0.001)

## 1.2.1 estimate the threshold of LOD using permutation test
operm <- scanone($cross,pheno.col=1:nph,n.perm=1000, verbose=FALSE)
LOD <- summary(operm,alpha=0.05) ## get LOD threshold
QTL.mr.test <- summary(QTL.mr,perms=operm,alpha=0.05,format="allpeak",pvalues=TRUE)

## 1.2.2 interval mapping
QTL.em <- scanone($cross, pheno.col=1:nph,method="em")
if (1){
QTL.hk <- scanone($cross, pheno.col=1:nph,method="hk")
QTL.ehk <- scanone($cross, pheno.col=1:nph,method="ehk")
$cross.imp <- sim.geno($cross,step=1,n.draws=64,error.prob=0.001)
QTL.imp <- scanone($cross.imp,pheno.col=1:nph, method="imp")
}
## 1.2.3 compostion interval mapping
#QTL.cim.20 <- cim($cross,pheno.col=1:nph,n.marcovar=3,window=20)
#QTL.cim.20.test <- summary(QTL.cim.20,perms=operm,pvalues=TRUE)

if (TRUE){ ###QTL effect switch
## 1.2.4 estimate QTL effect

cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
cross.imp.jit <- cross.imp
#cross.imp.jit <- jittermap(cross.imp)
#cross.imp.jit <- sim.geno(cross.imp.jit,step=1,n.draws=64,error.prob=0.001)

write.table("QTL summary by fitqtl","$opt{project}.QTL.fit.summary",row.names=FALSE, col.names=FALSE, sep="\t")
for (i in seq(4,length(QTL.mr.test),by=3)){
   print (names(QTL.mr.test)[i-1])
   write.table(names(QTL.mr.test)[i-1],"$opt{project}.QTL.fit.summary",append=TRUE,row.names=FALSE, col.names=FALSE ,sep="\t")
   select <- (QTL.mr.test[,i] <= 0.05) 
   chr <- QTL.mr.test[,1][select]
   pos <- QTL.mr.test[,i-2][select] 
   pheno <- floor(i/3)

   ###confident interval
   int_mk1 <- c()
   int_pos1 <- c()
   int_mk2 <- c()
   int_pos2 <- c()
   n <- (i-4)/3+1
   for (c in chr){ 
      interval <- lodint(QTL.mr, c, 1.5, lodcolumn=n, expandtomarkers=TRUE)
      int_mk1 <- append(int_mk1, paste(rownames(interval)[1],rownames(interval)[3], sep="-"))
      int_pos1 <- append(int_pos1, paste(interval\$pos[1],interval\$pos[3],sep="-"))
  
      interval_bayes <- bayesint(QTL.mr, c, 0.95, lodcolumn=n, expandtomarkers=TRUE)
      int_mk2 <- append(int_mk2, paste(rownames(interval_bayes)[1],rownames(interval_bayes)[3], sep="-"))
      int_pos2 <- append(int_pos2, paste(interval_bayes\$pos[1],interval_bayes\$pos[3],sep="-"))

   }

   if (length (chr) < 1){
      next
   }else if(length (chr) == 1){
      qtl <- makeqtl(cross.imp.jit,chr=chr,pos=pos)
      fit <- summary(fitqtl(cross.imp.jit,pheno.col=pheno,qtl=qtl,get.ests=TRUE))
      fitest <- fit\$ests[2:(length(chr)+1)]
      fitse  <- fit\$ests[(length(chr)+1+2):((length(chr)+1)*2)]
      fitlod <- fit\$result.full[10]
      fitvar <- fit\$result.full[13]
      fitpval <- fit\$result.full[16]
      
   }else if (length (chr) > 1){
      qtl <- makeqtl(cross.imp.jit,chr=chr,pos=pos)
      fit <- summary(fitqtl(cross.imp.jit,pheno.col=pheno,qtl=qtl,get.ests=TRUE))
      fitest <- fit\$ests[2:(length(chr)+1)]
      fitse  <- fit\$ests[(length(chr)+1+2):((length(chr)+1)*2)]
      fitlod <- fit\$result.drop[(length(chr)*2+1):(length(chr)*2+length(chr))]
      fitvar <- fit\$result.drop[(length(chr)*3+1):(length(chr)*3+length(chr))]
      fitpval <- fit\$result.drop[(length(chr)*5+1):(length(chr)*5+length(chr))]
   }
   fitsum <- cbind(chr,pos,int_pos1,int_mk1,int_pos2,int_mk2,fitlod,fitvar,fitpval,fitest,fitse)
   print (fitsum)
   write.table(fitsum,"$opt{project}.QTL.fit.summary",append=TRUE,row.names=FALSE, col.names=TRUE,sep="\t") 
}
}### QTL effect switch

## 2 write result into table and pdf
## 2.1 write QTL loci and LOD into table
write.table(QTL.mr,"$opt{project}.QTL.mr.table",sep="\t",col.names=NA)
write.table(LOD,"$opt{project}.QTL.mr.table.LOD_threshold",sep="\t",col.names=NA)
write.table(QTL.mr.test,"$opt{project}.QTL.mr.table.test",sep="\t",col.names=NA)
write.table(QTL.em,"$opt{project}.QTL.em.table",sep="\t",col.names=NA)
if(0){
write.table(QTL.hk,"$opt{project}.QTL.hk.table",sep="\t",col.names=NA)
write.table(QTL.ehk,"$opt{project}.QTL.ehk.table",sep="\t",col.names=NA)
write.table(QTL.imp,"$opt{project}.QTL.imp.table",sep="\t",col.names=NA)
}
#write.table(QTL.cim.20,"$opt{project}.QTL.cim.20",sep="\t",col.names=NA)
#write.table(QTL.cim.20.test,"$opt{project}.QTL.cim.20.test",sep="\t",col.names=NA)
#write.table(qtlsummary,"$opt{project}.QTL.cim.20.summary",sep="\t",col.names=NA)

## 2.2 write QTL curves into pdf
#pdf ("$opt{project}.QTL.pdf")
for (trait in 1:nph){
     traitname <- colnames(cross\$pheno)[trait]
     ## draw all chr together
     chr <- 1:12 
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (MR)", col="#000000$tranperant")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EM)", col="#FF0000$tranperant")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (HK)", col="#808000$tranperant")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EHK)", col="#008000$tranperant")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (IMP)", col="#0000FF$tranperant")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (CIM)", col="#800000$tranperant")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,lwd=0.2,ylab="LOD score",col="#000000$tranperant") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF0000$tranperant","#808000$tranperant","#008000$tranperant"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF$tranperant",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#800000$tranperant",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)
    
    ## draw chr by chr
    for (chr in 1:12){
        xtitle <- paste("Chr",chr," Map position (cM)",sep="")
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (MR)", col="#000000$tranperant")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EM)", col="#FF0000$tranperant")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (HK)", col="#808000$tranperant")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EHK)", col="#008000$tranperant")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (IMP)", col="#0000FF$tranperant")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (CIM)", col="#800000$tranperant")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,lwd=0.2,ylab="LOD score",col="#000000$tranperant") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF0000$tranperant","#808000$tranperant","#008000$tranperant"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF$tranperant",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#800000$tranperant",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)}
    }
}

dev.off()


R

open OUT, ">$opt{project}.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.R | /opt/linux/centos/7.x/x86_64/pkgs/R/3.2.0/bin/R --vanilla --slave`;

}
 
