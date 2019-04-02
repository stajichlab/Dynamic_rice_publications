#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"pileup:s","parent:s","help");


my $help=<<USAGE;
perl $0 --pileup GN60.Maq.p1.map.pileup --parent NB.RILs.dbSNP.SNPs.parents

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

convert($opt{pileup},$opt{parent}) unless (-e "$opt{pileup}.SNP");

sub convert
{
my ($file,$parent)=@_;
my $refparents=parents($parent);
my $rils=$1 if ($file=~/\/(GN\d+)\.Maq\.p1\.map\.pileup/ or $file=~/(GN\d+)\.Maq\.p1\.map\.pileup/); #../input/fastq/012/GN1.Maq.p1.map
my $maxdepth=50; #max depth of one SNP to avoid repetitive sequence regions
my $minbaseq=15; #min base quality for at least one base of one allele in SNP site
my $minreads=1;  #min of reads support a allele in SNP site
my $minsumbq=$minreads*15; #min of sum base quality for each allele in SNP site
my ($snpID);
open IN, "$file" or die "$!";
open OUT, ">$file.SNP" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] >= $maxdepth or $unit[3] == 0); ## possible repetitive regions or uncovered
    $unit[0]=~s/\D+//ig;
    $unit[0]=sprintf("%02d",$unit[0]);
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[2];
    #print "$snpID\t$_\n";
    next unless (exists $refparents->{$snpID}); ## SNPs need exists in parents SNPs
    #remove non-base characters
    if ($unit[4] =~ /([-|+]+[0-9]+[ACGTNacgtn]+)/){
        $unit[4] =~ s/([-|+]+[0-9]+[ACGTNacgtn]+)//g;
    }
    if ($unit[4] =~ /([^.,ACGTNacgtn]+)/){
        $unit[4] =~ s/([^.,ACGTNacgtn]+)//g;
    }

    my @base=split("",$unit[4]);
    my @qual=split("",$unit[5]);
    my %allele;
    for(my $i=0;$i<@base;$i++){
       if ($base[$i]=~/\,/ or $base[$i]=~/\./){
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$unit[2]}}, $qscore ) if ($qscore >= 10); ### do not count if quality is too poor
       }elsif($base[$i]=~/[acgtACGT]/){
          $base[$i]=~tr/atcg/ATCG/;
          my $qscore=ord($qual[$i])-33;
          push (@{$allele{$base[$i]}}, $qscore ) if ($qscore >= 10); ### do not count if quality is too poor
       }
    }
   
    #remove minor allele, which probably due to sequencing error or somatic mutations due to high coverage
    my $flag = keys %allele;
    my $coverage = 0;
    for my $tempb (keys %allele){
       my $tempd = @{$allele{$tempb}};
       $coverage += $tempd;
    }
    for my $tempb (keys %allele){
       my $tempd = @{$allele{$tempb}};
       if ($tempd/$coverage <= 0.2 and $coverage >= 5){
          delete $allele{$tempb};
       }
    }

    #print "$snpID\t$flag\tcheck2\n";
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
}




####
sub parents
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/V/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
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
 
