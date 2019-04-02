#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"trait:s","help");


my $help=<<USAGE;
perl $0 --trait ../input/trait/May28_2013.RIL.trait.table
Parse the trait table and generate trait matrix for parent and RILs. Draw trait distribution for each trait.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable($opt{trait});

#####################
sub readtable
{
my ($file)=@_;
my %hash;
my @sample;
my %traits;
open IN, "$file" or die "$!";
my @head=split("\t",<IN>);
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $id = $unit[1];
    push @sample, $id;
    for (my $i=4;$i<@unit;$i=$i+3){
        $trait="$1" if $head[$i]=~/Mean\s*(.*)\-All Rep/;
        $traits{$trait}=$i;
        $hash{$id}{$trait}=[$unit[$i],$unit[$i+1],$unit[$i+2]];
    }
}
close IN;

##store traits
my @traits=sort {$traits{$a} <=> $traits{$b}} keys %traits;
my @temp;
push @temp, @traits;
for(my $i=0;$i<@temp;$i++){
   $temp[$i]="\"$temp[$i]\"";
}
my $traits4r=join(",",@temp);

##store parents trait value
my @parents;
for(my $i=0;$i<@traits;$i++){
   $parents[$i]=[@{$hash{"Nipponbare"}{$traits[$i]}},@{$hash{"HEG4"}{$traits[$i]}}];
}

open OUT, ">$file.QTL.parents.txt" or die "$!";
for (my $i=0;$i<@parents;$i++){
   print OUT "$traits[$i]\t$parents[$i]->[0]\t$parents[$i]->[1]\t$parents[$i]->[2]\t$parents[$i]->[3]\t$parents[$i]->[4]\t$parents[$i]->[5]\n";
}
close OUT;

##write out QTL trait table
open OUT, ">$file.QTL.trait.txt" or die "$!";
     ##head line
     print OUT "Sample";
     for(my $i=0;$i<@traits;$i++){
        print OUT "\t$traits[$i]";
     }
     print OUT "\n";
     ##trait line
     for(my $i=0;$i<@sample;$i++){
        next unless ($sample[$i]=~/^GN/); ## RILs, not Nipporbare, HEG4 or A1XXX
        print OUT "$sample[$i]";
        for(my $j=0;$j<@traits;$j++){
           print OUT "\t$hash{$sample[$i]}{$traits[$j]}->[0]";
        }
        print OUT "\n";
     }
close OUT;


##draw QTL trait distribution
my $Rcmd=<<Rscript;
pdf("DrawQTLtrait.pdf")
par(mfrow=c(3,2))

error.bar <- function(x, y, upper,color, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x-upper,y, x+upper, y, angle=90, code=3,col=color, length=length, ...)
}


traits <- c($traits4r)
read.table("$file.QTL.trait.txt",sep="\\t",skip=1) -> trait
read.table("$file.QTL.parents.txt",sep="\\t") -> parent
for (i in 2:ncol(trait)){

    hist(trait[,i],breaks=10,plot=FALSE) -> xh
    step=20;
    yadj=max(xh\$counts)*1.2/20;
    barplot(xh\$counts,axes=FALSE,border=F,ylim=c(0,max(xh\$counts)*1.2),col=c("cornflowerblue"),xlab=traits[i-1],ylab="# of RILs") -> xx
    for (j in 1:length(xx)) { # adj can not take vector, so we use loops to add text
      text(xx[j],-yadj,labels=floor(xh\$mids[j]),cex=1,srt=65,adj=c(1,1),xpd=TRUE)
    }
    axis(1,c(0,max(xx)+0.5),line=0.2,labels=c("",""))
    axis(2,seq(0,max(xh\$counts)*1.4,by=step),cex=1.2)
    index <- i-1
    pix <- (max(xx)-min(xx)+1)/(max(xh\$mids)-min(xh\$mids)+1)
    x1 <- (parent[index,2]-min(xh\$mids)+1)*pix
    x2 <- (parent[index,5]-min(xh\$mids)+1)*pix
    points(x1,max(xh\$counts)*1.1,col="red")
    points(x2,max(xh\$counts)*1.1,col="blue") 
    
    ex <- c(x1,x2)
    ey <- c(max(xh\$counts)*1.1,max(xh\$counts)*1.1)
    sd <- c(parent[index,4]*pix,parent[index,7]*pix)
    ec <- c("red","blue") 
    error.bar(ex,ey,sd,ec)
    
    legend(max(xx)-4,max(xh\$counts)*1,bty="n",lty=c(1,1),cex=0.7,c("Nipponbare","HEG4"),col=c("red","blue")) 
}
dev.off()
Rscript

open OUT, ">DrawQTLtrait.R" or die "$!";
     print OUT "$Rcmd\n";
close OUT;
`cat DrawQTLtrait.R | R --vanilla --slave`;

return \%hash;
}
 
