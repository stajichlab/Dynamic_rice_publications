#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"RIL:s","bam:s","help");


my $help=<<USAGE;
perl $0 --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Bam --bam /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @fq=glob("$opt{RIL}/RIL*.bam");
my %rils;
my %dup;
for (my $i=0;$i<@fq;$i++){
   #print "$fq[$i]\n"; 
   if ($fq[$i]=~/.*\/RIL(\d+)\_.*\.recal\.bam/){ #RIL100_0_TGACCA_FC0813L3.recal.bam
      if (exists $rils{$1}){ ## skip if already have fastq record (these line have multi lane sequenced)
         $dup{$1}  = 1;
         next;
      }
      $bam = "$opt{bam}/GN$1.bam";
      #print "$fq[$i]\t$bam\n";
      #next if (-e "$bam");
      #`ln -s $fq[$i] $bam`;
      $rils{$1}=1;
   }
}
foreach (keys %dup){
   print "$_\n";

}


############
#Illumina/RIL10_0/RIL10_0_GAGTGG_p1.fq
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}
 
