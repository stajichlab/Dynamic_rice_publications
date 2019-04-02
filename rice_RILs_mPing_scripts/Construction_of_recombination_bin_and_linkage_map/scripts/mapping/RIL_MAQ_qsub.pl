#!/usr/bin/perl
=header
The scripts is designed to map fastq of RILs to reference using MAQ.
Use full path of files in command line.
--ref:   reference sequence of parent1
--fastq:  dir of fastq of pair-end reads
--project: project name that used for result file 
=cut

use Getopt::Long;
use FindBin qw($Bin $Script);

my %opt;
GetOptions(\%opt,"ref:s","fastq:s","verbose","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.2X\n";
   exit();
}

$opt{project} ||= "NBxHEG4";

my $maq="/opt/tyler/bin/maq";
my @fastq=glob("$opt{fastq}/*_1.fq");

## build reference index
unless (-e "$opt{ref}.bfa"){
   `$maq fasta2bfa $opt{ref} $opt{ref}.bfa`;
}


my @cmd;
## mapping by MAQ
open OUT, ">MAQ.sampleRIL.list" or die "$!";
foreach my $fq1 (@fastq){
   if ($fq1=~/(.*)_1.fq/){
      my $line=$1;
      my $fq2="$line\_2.fq";
      if ($line=~/\/(\w+)$/){
         print OUT "$1\n";
      }
      unless (-e "$line.Maq.p1.map.pileup"){
         print "Write pair-end mapping $fq1 and $fq2 by Maq!\n";
         push @cmd, "$maq fastq2bfq $fq1 $fq1.bfq";
         push @cmd, "$maq fastq2bfq $fq2 $fq2.bfq";
         push @cmd, "$maq match -a 900 $line.Maq.p1.map $opt{ref}.bfa $fq1.bfq $fq2.bfq";
         push @cmd, "$maq pileup -vP -q 40 $opt{ref}.bfa $line.Maq.p1.map | awk '\$4 > 0' > $line.Maq.p1.map.pileup";
         print "Done!\n";
      }
   }
}
close OUT;

if(@cmd > 4){
open OUT, ">maq.sh" or die "$!";
for(my $i=0; $i<@cmd; $i++){
   print OUT "$cmd[$i]\n";
}
close OUT;
`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --lines 4 --resource nodes=1:ppn=1,mem=2G,walltime=100:00:00 maq.sh`;
}
