#!/usr/bin/perl
=header
The scripts is designed to generate SNPs for RILs using Maq pileup command. The results will be input for MPR package developed by Qifa Zhang's lab.
--ref:   reference sequence
--fastq: dir of fastq
--split: split chromosome for pileup
--parents: parents SNPs, which could be predicted in this script by parsepileup or provide by user. for example, 0500000526A A G.
--project: project name that used for result file 
=cut

use Getopt::Long;
use FindBin qw($Bin $Script);

my %opt;
GetOptions(\%opt,"ref:s","fastq:s","trait:s","split","parents:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --ref ../input/reference/Pseudomolecules_Nipponbare.Build4.0.fasta --fastq ../input/testfastq\n";
   exit();
}

$opt{project} ||= "NB.RILs.dbSNP";

my $maq="/opt/tyler/bin/maq";
my @map=glob("$opt{fastq}/GN*.Maq.p1.map");
my $map=join(" ",@map);
my $cmd1="$maq mapmerge $opt{project}.Maq.map $map";
my $cmd2;
unless ($opt{split}){
   $cmd2="$maq pileup -v -q 40 $opt{ref}.bfa $opt{project}.Maq.map | awk '\$4 > 0' > $opt{project}.Maq.map.pileup";
}else{
   $cmd2="$maq pileup -v -q 40 $opt{ref}.bfa $opt{project}.Maq.map | awk '\$4 > 20 && \$4 < 100' | perl -n -e '/^(chromosome\\d+).*?\$/ and open FH, \">>$opt{project}.Maq.map.\$1.pileup\";print FH;'";
}
`$cmd1` unless (-e "$opt{project}.Maq.map" or $opt{parents});
`$cmd2` unless (-e "$opt{project}.Maq.map.chromosome12.pileup" or $opt{parents}); ## clean chromosome*.pileup file if want to generate new one
#print "$cmd2\n";
if ($opt{split}){
   for(my $i=1; $i<=12;$i++){
      my $chr=sprintf("%02d");
      my $pileup="$opt{project}.Maq.map.chromosome$chr.pileup";
      parsepileup($pileup) unless ($opt{parents}); ## generate parents SNPs from pileup file of Maq, skip if parents provide by user
   }
}else{
   my $pileup="$opt{project}.Maq.map.pileup";
   parsepileup($pileup) unless ($opt{parents}); ## generate parents SNPs from pileup file of Maq, skip if parents provide by user
}
my $parents=$opt{parents} ? $opt{parents} : "$opt{project}.SNPs.parents";
#my $parents="$opt{project}.SNPs.parents";
#if ($opt{parents} and $opt{parents} ne $parents){
#   `ln -s $opt{parents} $opt{project}.SNPs.parents`;
#}
#RILgenotype(\@map,$parents); ## call SNPs for RILs based on SNPs of parents, which means only the SNPs exists parents will be called
RILgenotype_effecient(\@map,$parents);  ### use this effeciency subfunction instead, which read pileup.SNP if exists.


####
sub parents
{
my ($file)=@_;
my %hash;
open OUT, ">$opt{project}.SNPs.Markers" or die "$!"; #0500007370C     T       C
print OUT "SNP_id\tAllele\n";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/V/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2]];
    my $snp=substr($unit[0],-1,1);
    print OUT "$unit[0]\t$snp\n";
}
close IN;
close OUT;
return \%hash;
}


sub RILgenotype
{
my ($map,$parents)=@_;
my $refparents=parents($parents);
my @RIL;
my %hash;
open TEMP, ">temp.pileup.list" or die "$!"; 
for(my $i=0; $i< @$map; $i++){
my $rils=$1 if ($map->[$i]=~/\/(GN\d+)\.Maq\.p1\.map/); #../input/fastq/012/GN1.Maq.p1.map
push @RIL, $rils;
print TEMP "$rils\n";
`$maq pileup -vP -q 40 $opt{ref}.bfa $map->[$i] | awk '\$4 > 0' > $map->[$i].pileup` unless (-e "$map->[$i].pileup");
#`grep "chromosome05" $map->[$i].pileup | awk '\$4 > 0' > $map->[$i].pileup.chr05` unless (-e "$map->[$i].pileup.chr05");
my $maxdepth=20; #max depth of one SNP to avoid repetitive sequence regions
my $minbaseq=20; #min base quality for at least one base of one allele in SNP site
my $minreads=1;  #min of reads support a allele in SNP site
my $minsumbq=$minreads*20; #min of sum base quality for each allele in SNP site
my ($snpID);
open IN, "$map->[$i].pileup" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] >= $maxdepth or $unit[3] == 0); ## possible repetitive regions or uncovered
    $unit[0]=~s/\D+//ig;
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[2];
    #print "$snpID\t$_\n";
    next unless (exists $refparents->{$snpID}); ## SNPs need exists in parents SNPs
    my @base=split("",$unit[4]);
    my @qual=split("",$unit[5]);
    my %allele;
    for(my $i=1;$i<@base;$i++){
       if ($base[$i]=~/\,/ or $base[$i]=~/\./){
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$unit[2]}}, $qscore ) if ($qscore >= 15); ### do not count if quality is too poor
       }elsif($base[$i]=~/[acgtACGT]/){
          $base[$i]=~tr/atcg/ATCG/;
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$base[$i]}}, $qscore ) if ($qscore >= 15); ### do not count if quality is too poor
       }
    }
    #print "$snpID\tcheck2\n";
    if (keys %allele == 1){  ##only have one alleles in SNP site, what if the coverage if deep, say ~6. how to deal with heterozygous or can not deal with. Random one?
       my @SNP=keys %allele;
       #print "$snpID\t$SNP[0]\t$refparents->{$snpID}->[0]\t$refparents->{$snpID}->[1]\n";
       next unless ($SNP[0] eq $refparents->{$snpID}->[0] or $SNP[0] eq $refparents->{$snpID}->[1]); ## SNPs need to one of the parents SNPs
       my %maxbq; my %sumbq; my %maxreads;
       foreach my $a (keys %allele){
          $maxreads{$a}=@{$allele{$a}};
          $maxbq{$a}=max(\@{$allele{$a}});
          $sumbq{$a}=sum(\@{$allele{$a}});
          #print "$a\t$maxreads{$a}\t$maxbq{$a}\t$sumbq{$a}\n";
       } 
       my $flag1;my $flag2; my $flag3;
       foreach my $v (values %maxreads){
          $flag1++ if $v >= $minreads; ## reads depth larger than $minreads==4 for each allele
       }
       foreach my $q (values %maxbq){ 
          $flag2++ if $q >= $minbaseq; ## at least have a base quality larger than $minbaseq==20 for each allele  
       }
       foreach my $q (values %sumbq){
          $flag3++ if $q >= $minsumbq; ## sum of base quality larger than $minsumbq==60 for each allele
       }
       #print "$snpID\t$rils\t$SNP[0]\t$flag1\t$flag2\t$flag3\n";
       if ($flag1 == 1 and $flag2 == 1 and $flag3 ==1){ ## need to meet all three criteria for both allele
          #print OUT "$snpID\t$SNP\n";
          #print "$snpID\t$rils\t$SNP[0]\n";
          $hash{$snpID}{$rils}=$SNP[0];
       }
    }
}#while loop
close IN;
}# for loop
close TEMP;

print "output genotype\n";
## write SNPs of RILs into file
open OUT, ">$opt{project}.SNPs.RILs" or die "$!";
open OUT1, ">$opt{project}.SNPs.sub.parents" or die "$!";
open OUT2, ">$opt{project}.SNPs.sub.Marker";
print OUT1 "V1\tV2\n";
print OUT2 "SNP_id\tAllele\n";
my $head=join("\t",@RIL);
print OUT "$head\n";
foreach my $snp (sort {$a <=> $b} keys %$refparents){
    #print OUT "$snp\t";
    my %allele;
    my @rilsgeno;
    for (my $i=0; $i<@RIL; $i++ ){
        my $allele0= $hash{$snp}{$RIL[$i]} ? $hash{$snp}{$RIL[$i]} : "NA";
        #print "$snp\t$RIL[$i]\t$hash{$snp}{$RIL[$i]}\tALLELE:$allele\n";
        push @rilsgeno, $allele0;
        $allele{$allele0}=1 if ($allele0 ne "NA");
    }
    my $temp=join("\t",@rilsgeno);
    print OUT "$snp\t$temp\n" if (keys %allele == 2); ## keep only these biallelic lines
    print OUT1 "$snp\t$refparents->{$snp}->[0]\t$refparents->{$snp}->[1]\n" if (keys %allele == 2);
    print OUT2 "$snp\t$refparents->{$snp}->[0]\n" if (keys %allele == 2);
}
close OUT;
close OUT1;
close OUT2;
}# end of sub function


#### parsing pileup file is slow. So we store the parse result in file and read it when it exists.
sub RILgenotype_effecient
{
my ($map,$parents)=@_;
my $refparents=parents($parents);
my @RIL;
my %hash;
open TEMP, ">temp.pileup.list" or die "$!"; 
for(my $i=0; $i< @$map; $i++){
my $rils=$1 if ($map->[$i]=~/\/(GN\d+)\.Maq\.p1\.map/); #../input/fastq/012/GN1.Maq.p1.map
push @RIL, $rils;
print TEMP "$rils\n";
`$maq pileup -vP -q 40 $opt{ref}.bfa $map->[$i] | awk '\$4 > 0' > $map->[$i].pileup` unless (-e "$map->[$i].pileup");
#`grep "chromosome05" $map->[$i].pileup | awk '\$4 > 0' > $map->[$i].pileup.chr05` unless (-e "$map->[$i].pileup.chr05");

unless (-e "$map->[$i].pileup.SNP"){
my $maxdepth=20; #max depth of one SNP to avoid repetitive sequence regions
my $minbaseq=20; #min base quality for at least one base of one allele in SNP site
my $minreads=1;  #min of reads support a allele in SNP site
my $minsumbq=$minreads*20; #min of sum base quality for each allele in SNP site
my ($snpID);
open IN, "$map->[$i].pileup" or die "$!";
open OUT, ">$map->[$i].pileup.SNP" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] >= $maxdepth or $unit[3] == 0); ## possible repetitive regions or uncovered
    $unit[0]=~s/\D+//ig;
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[2];
    #print "$snpID\t$_\n";
    next unless (exists $refparents->{$snpID}); ## SNPs need exists in parents SNPs
    my @base=split("",$unit[4]);
    my @qual=split("",$unit[5]);
    my %allele;
    for(my $i=1;$i<@base;$i++){
       if ($base[$i]=~/\,/ or $base[$i]=~/\./){
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$unit[2]}}, $qscore ) if ($qscore >= 15); ### do not count if quality is too poor
       }elsif($base[$i]=~/[acgtACGT]/){
          $base[$i]=~tr/atcg/ATCG/;
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$base[$i]}}, $qscore ) if ($qscore >= 15);
       }
    }
    #print "$snpID\tcheck2\n";
    if (keys %allele == 1){  ##only have one alleles in SNP site, what if the coverage if deep, say ~6. how to deal with heterozygous or can not deal with. Random one?
       my @SNP=keys %allele;
       #print "$snpID\t$SNP[0]\t$refparents->{$snpID}->[0]\t$refparents->{$snpID}->[1]\n";
       next unless ($SNP[0] eq $refparents->{$snpID}->[0] or $SNP[0] eq $refparents->{$snpID}->[1]); ## SNPs need to one of the parents SNPs
       my %maxbq; my %sumbq; my %maxreads;
       foreach my $a (keys %allele){
          $maxreads{$a}=@{$allele{$a}};
          $maxbq{$a}=max(\@{$allele{$a}});
          $sumbq{$a}=sum(\@{$allele{$a}});
          #print "$a\t$maxreads{$a}\t$maxbq{$a}\t$sumbq{$a}\n";
       } 
       my $flag1;my $flag2; my $flag3;
       foreach my $v (values %maxreads){
          $flag1++ if $v >= $minreads; ## reads depth larger than $minreads==4 for each allele
       }
       foreach my $q (values %maxbq){ 
          $flag2++ if $q >= $minbaseq; ## at least have a base quality larger than $minbaseq==20 for each allele  
       }
       foreach my $q (values %sumbq){
          $flag3++ if $q >= $minsumbq; ## sum of base quality larger than $minsumbq==60 for each allele
       }
       #print "$snpID\t$rils\t$SNP[0]\t$flag1\t$flag2\t$flag3\n";
       if ($flag1 == 1 and $flag2 == 1 and $flag3 ==1){ ## need to meet all three criteria for both allele
          #print OUT "$snpID\t$SNP\n";
          #print "$snpID\t$rils\t$SNP[0]\n";
          $hash{$snpID}{$rils}=$SNP[0];
          print OUT "$snpID\t$rils\t$SNP[0]\n"; 
       }
    }
}#while loop
close IN;
close OUT;
}else{ ###### *.pileup.SNP exists, just parse this file and store SNP in %hash
    open IN, "$map->[$i].pileup.SNP" or die "$!";
    while(<IN>){ 
        chomp $_;
        next if ($_=~/^$/);
        my @unit=split("\t",$_);
        $hash{$unit[0]}{$unit[1]}=$unit[2];
    }
    close IN;
}


}# for loop
close TEMP;



print "output genotype\n";
## write SNPs of RILs into file
open OUT, ">$opt{project}.SNPs.RILs" or die "$!";
open OUT1, ">$opt{project}.SNPs.sub.parents" or die "$!";
open OUT2, ">$opt{project}.SNPs.sub.Marker";
print OUT1 "V1\tV2\n";
print OUT2 "SNP_id\tAllele\n";
my $head=join("\t",@RIL);
print OUT "$head\n";
foreach my $snp (sort {$a <=> $b} keys %$refparents){
    #print OUT "$snp\t";
    my %allele;
    my @rilsgeno;
    for (my $i=0; $i<@RIL; $i++ ){
        my $allele0= $hash{$snp}{$RIL[$i]} ? $hash{$snp}{$RIL[$i]} : "NA";
        #print "$snp\t$RIL[$i]\t$hash{$snp}{$RIL[$i]}\tALLELE:$allele\n";
        push @rilsgeno, $allele0;
        $allele{$allele0}=1 if ($allele0 ne "NA");
    }
    my $temp=join("\t",@rilsgeno);
    print OUT "$snp\t$temp\n" if (keys %allele == 2); ## keep only these biallelic lines
    print OUT1 "$snp\t$refparents->{$snp}->[0]\t$refparents->{$snp}->[1]\n" if (keys %allele == 2);
    print OUT2 "$snp\t$refparents->{$snp}->[0]\n" if (keys %allele == 2);
}
close OUT;
close OUT1;
close OUT2;
}# end of sub function






#chromosome01    1019    A       16      @,,.,...,,g.,,,,,       @JJ,J#A:HH#(HHHDC       @~~~~z~~~~k~~~~z~       19,18,84,17,85,87,89,12,11,92,92,9,9,9,4,3,
sub parsepileup
{
my ($file)=@_;
my $maxdepth=100; #max depth of one SNP to avoid repetitive sequence regions, test value at 30
my $minbaseq=20; #min base quality for at least one base of one allele in SNP site, test value at 20
my $minreads=4;  #min of reads support a allele in SNP site, test value at 2
my $minsumbq=$minreads*20; #min of sum base quality for each allele in SNP site, should be minreads*20
my ($snpID);
open OUT, ">>$opt{project}.SNPs.parents" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] >= $maxdepth or $unit[3] == 0); ## possible repetitive regions or uncovered
    $unit[0]=~s/\D+//ig;
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[2];
    my @base=split("",$unit[4]);
    my @qual=split("",$unit[5]);
    my %allele;
    for(my $i=1;$i<@base;$i++){
       if ($base[$i]=~/\,/ or $base[$i]=~/\./){
          my $qscore=ord($qual[$i])-33;
          push @{$allele{$unit[2]}}, $qscore;
       }else{
          $base[$i]=~tr/atcg/ATCG/;
          my $qscore=ord($qual[$i])-33;
          if ($qscore >= 20){  ### only if a base have quality score >= 20, we count it as a effective base
             push @{$allele{$base[$i]}}, $qscore;
          }
       }
    }
    if (keys %allele == 2){  ##only have two effective alleles in SNP site
       #my $SNP=join("\t",keys %allele);
       my $SNP;
       my $flag0;
       #print OUT "$snpID\t$SNP\tcheck\n";
       my %maxbq; my %sumbq; my %maxreads;
       foreach my $a (keys %allele){
          if ($a eq $unit[2]){
             $flag0=1; ## one allele is reference allele
          }else{
             $SNP=$a;
          }
          $maxreads{$a}=@{$allele{$a}};
          $maxbq{$a}=max(\@{$allele{$a}});
          $sumbq{$a}=sum(\@{$allele{$a}});
          #print OUT "$a\t$maxreads{$a}\t$maxbq{$a}\t$sumbq{$a}\n";
       } 
       my $flag1; my $flag2; my $flag3;
       foreach my $v (values %maxreads){
          $flag1++ if $v >= $minreads; ## reads depth larger than $minreads==4 for each allele
       }
       foreach my $q (values %maxbq){ 
          $flag2++ if $q >= $minbaseq; ## at least have a base quality larger than $minbaseq==20 for each allele  
       }
       foreach my $q (values %sumbq){
          $flag3++ if $q >= $minsumbq; ## sum of base quality larger than $minsumbq==60 for each allele
       }
       #print OUT "$flag1\t$flag2\t$flag3\n";
       if ($flag0 == 1 and $flag1 == 2 and $flag2 == 2 and $flag3 ==2){ ## need to meet all three criteria for both allele
          print OUT "$snpID\t$unit[2]\t$SNP\n";
       }
    }
}
close IN;
close OUT;
}


sub max
{
my ($num)=@_;
my $max=0;
foreach  (@$num) {
        next if ($_ eq "NA");
        $max = $_ > $max ? $_ : $max;
}
return $max;
}


sub sum
{
my ($num)=@_;
my $loop=0;
my $total;
foreach (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
return 0 if ($loop == 0);
return $total;
}


