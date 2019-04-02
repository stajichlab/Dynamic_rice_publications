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
GetOptions(\%opt,"fastq:s","trait:s","split","parents:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --fastq ../input/fastq/RILs_ALL_bam --parents NB.RILs.dbSNP.SNPs.parents\n";
   exit();
}

$opt{project} ||= "NB.RILs.dbSNP";

my @map=glob("$opt{fastq}/GN*.Maq.p1.map.pileup.SNP");
my $parents=$opt{parents} ? $opt{parents} : "$opt{project}.SNPs.parents";
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


#### parsing pileup file is slow. So we store the parse result in file and read it when it exists.
sub RILgenotype_effecient
{
my ($map,$parents)=@_;
my $refparents=parents($parents);
my @RIL;
my %hash;
open TEMP, ">temp.pileup.list" or die "$!"; 
for(my $i=0; $i< @$map; $i++){
    my $rils=$1 if ($map->[$i]=~/\/(GN\d+)\.Maq\.p1\.map\.pileup\.SNP/); #path/GN87.Maq.p1.map.pileup.SNP
    push @RIL, $rils;
    print TEMP "$rils\n";

    open IN, "$map->[$i]" or die "$!";
    while(<IN>){ 
        chomp $_;
        next if ($_=~/^$/);
        my @unit=split("\t",$_);
        $hash{$unit[0]}{$unit[1]}=$unit[2];
    }
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






