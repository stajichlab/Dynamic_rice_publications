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
my $bwa="/opt/tyler/bin/bwa";
my $samtools="/usr/local/bin/samtools";
my $java="/opt/java/jdk1.6.0_43/bin/java";
my $rmdup="/opt/picard/1.81/MarkDuplicates.jar";
my @fastq=glob("$opt{fastq}/*_1.fq");

## build reference index
unless (-e "$opt{ref}.sa"){
   `$bwa index $opt{ref}`;
   `$samtools faidx $opt{ref}`;
}


my @cmd;
## mapping by BWA
open OUT, ">BWA.sampleRIL.list" or die "$!";
foreach my $fq1 (@fastq){
   if ($fq1=~/(.*)_1.fq/){
      my $line=$1;
      my $fq2="$line\_2.fq";
      if ($line=~/\/(\w+)$/){
         print OUT "$1\n";
      }
      unless (-e "$line.Maq.p1.map.pileup.SNP"){
         print "Write pair-end mapping $fq1 and $fq2 by BWA!\n";
         #push @cmd, "$maq fastq2bfq $fq1 $fq1.bfq";
         #push @cmd, "$maq fastq2bfq $fq2 $fq2.bfq";
         #push @cmd, "$maq match -a 900 $line.Maq.p1.map $opt{ref}.bfa $fq1.bfq $fq2.bfq";
         #push @cmd, "$maq pileup -vP -q 40 $opt{ref}.bfa $line.Maq.p1.map | awk '\$4 > 0' > $line.Maq.p1.map.pileup";
         
         push @cmd, "$bwa aln -t 5 $opt{ref} $fq1 > $fq1.sai";
         push @cmd, "$bwa aln -t 5 $opt{ref} $fq2 > $fq2.sai";
         push @cmd, "$bwa sampe $opt{ref} $fq1.sai $fq2.sai $fq1 $fq2 > $line.sam";
         push @cmd, "$samtools view -bS -o $line.raw.bam $line.sam";
         push @cmd, "$samtools sort $line.raw.bam $line.sort";
         push @cmd, "$java -Xmx2G -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$line.sort.bam OUTPUT=$line.bam METRICS_FILE=$line.dupli";
         push @cmd, "$samtools mpileup -q 40 -Q 15 -f $opt{ref} $line.bam | awk '\$4 > 0' > $line.Maq.p1.map.pileup";
         push @cmd, "rm $line.sam $line.raw.bam $line.sort.bam $fq1.sai $fq2.sai";
         print "Done!\n";
      }
   }
}
close OUT;

if(@cmd > 8){
open OUT, ">bwa.sh" or die "$!";
for(my $i=0; $i<@cmd; $i++){
   print OUT "$cmd[$i]\n";
}
close OUT;
`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --lines 8 --resource nodes=1:ppn=5,mem=3G,walltime=100:00:00 bwa.sh`;
}
