#!/usr/bin/perl
=header
The scripts is designed to run bwa to map solexa sequencing read to reference genome.
--ref:   reference sequence
--1:     paired read with -2, SRR034638_1.fastq
--2:     paired read with -1, SRR034638_2.fastq
--tool:  mapping tools: bwa, maq, ssaha, soap
-project: project name that used for result file 
=cut

use Getopt::Long;
my %opt;
GetOptions(\%opt,"ref:s","1:s","2:s","tool:s","min:s","max:s","cpu:s","bam","verbose","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref ../input/all.con -1 1.fq -2 2.fq -min 0 -max 500 -cpu 12 -tool soap\n";
   exit();
}

$opt{project} ||= "out";
$opt{tool} ||= "bwa";
$opt{cpu} ||=12;
$opt{min} ||= 0;
$opt{max} ||= 500; 

my $bwa="/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin";
my $soap="/usr/local/bin";
my $ssaha="/home/jfchen/software/ssaha2_v2.5.5_x86_64/";
my $maq="/opt/tyler/bin/maq";
 
my $SAMtool="/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools";
#my $rmdup="/opt/picard/1.81/MarkDuplicates.jar";
my $rmdup="/opt/linux/centos/7.x/x86_64/pkgs/picard/1.130/bin/picard MarkDuplicates";

if (exists $opt{1} and exists $opt{2}){
   if ($opt{tool}=~/bwa/){
      print "Run pair-end mapping by BWA!\n";
      unless (-e "$opt{ref}.sa"){
         `$bwa/bwa index $opt{ref} > $opt{project}.index.log 2> $opt{project}.index.log2`;
      }
      print "Align Read 1!\n";
      `$bwa/bwa aln -t $opt{cpu} $opt{ref} $opt{1} > $opt{project}.1.sai 2> $opt{project}.1.bwa.log2`;
      print "Align Read 2!\n";
      `$bwa/bwa aln -t $opt{cpu} $opt{ref} $opt{2} > $opt{project}.2.sai 2> $opt{project}.2.bwa.log2`;
      print "Pairing!\n";
      `$bwa/bwa sampe -a $opt{max} $opt{ref} $opt{project}.1.sai $opt{project}.2.sai $opt{1} $opt{2} > $opt{project}.sam 2> $opt{project}.sampe.log2`;
      print "SAM 2 BAM!\n";
      `$SAMtool view -bS -@ $opt{cpu} -o $opt{project}.raw.bam $opt{project}.sam > $opt{project}.convert.log 2> $opt{project}.convert.log2`;
      print "Sort Bam!\n";
      `$SAMtool sort -@ $opt{cpu} $opt{project}.raw.bam $opt{project}.sort > $opt{project}.sort.log 2> $opt{project}.sort.log2`;
      print "Remove duplicate!\n";
      `$rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT INPUT=$opt{project}.sort.bam OUTPUT=$opt{project}.bam METRICS_FILE=$opt{project}.dupli > $opt{project}.rmdup.log 2> $opt{project}.rmdup.log2`;
      `$SAMtool index $opt{project}.bam`;
      unless ($opt{verbose}){
          `rm $opt{project}.sam $opt{project}.raw.bam $opt{project}.sort.bam`;
          `rm $opt{project}.*.log* $opt{project}.1.* $opt{project}.2.*`;
      }
      print "Done!\n";
   }
}


