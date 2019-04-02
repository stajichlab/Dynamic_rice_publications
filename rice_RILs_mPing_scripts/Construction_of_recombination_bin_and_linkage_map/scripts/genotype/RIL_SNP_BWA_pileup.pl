#!/usr/bin/perl
=header
This script is design to run pileup and SNP calling using qsub, which speed up the process in RIL_SNP_MAQ.pl.
--ref:   reference sequence
--parent: parents file
--fastq: dir of fastq
=cut

use Getopt::Long;
use FindBin qw($Bin $Script);

my %opt;
GetOptions(\%opt,"ref:s", "parent:s", "fastq:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --ref ../input/reference/Pseudomolecules_Nipponbare.Build4.0.fasta --parent NB.RILs.dbSNP.SNPs.parents --fastq ../input/testfastq\n";
   exit();
}

my @map=glob("$opt{fastq}/GN*.bam");

pileup(\@map);


sub pileup
{
my ($map)=@_;
my @cmd;
my $samtools="/opt/samtools-0.1.16/samtools"; #version used here should be the same with version used for index genome sequence
#my $convert="/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/pileup2SNP.pl";
my $convert="/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/pileup2SNP_lowquality_highcoverage.pl";
open LIST, ">BWA.sampleRIL.list" or die "$!";
for(my $i=0; $i< @$map; $i++){
   my $prefix=$1 if ($map->[$i]=~/(.*)\.bam$/);
   if ($prefix=~/\/(\w+)$/){
       print LIST "$1\n";
   }
   unless (-e "$prefix.Maq.p1.map.pileup.SNP"){
       push @cmd, "$samtools mpileup -q 40 -Q 15 -d 8000 -f $opt{ref} $map->[$i] | awk '\$4 > 0' > $prefix.Maq.p1.map.pileup";
       push @cmd, "perl $convert --pileup $prefix.Maq.p1.map.pileup --parent $opt{parent}";
       push @cmd, "rm $prefix.Maq.p1.map.pileup";
   }
   open OUT, ">pileup.sh" or die "$!";
   for(my $i=0; $i<@cmd; $i++){
      print OUT "$cmd[$i]\n";
   }
   close OUT;
}# for loop
close LIST;
#`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --maxjob 10 --resource nodes=1:ppn=1,mem=5G,walltime=100:00:00 pileup.sh`;
`perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --convert no --lines 3 --maxjob 20 --resource nodes=1:ppn=1,mem=5G,walltime=100:00:00 pileup.sh`;

}# end of sub function


