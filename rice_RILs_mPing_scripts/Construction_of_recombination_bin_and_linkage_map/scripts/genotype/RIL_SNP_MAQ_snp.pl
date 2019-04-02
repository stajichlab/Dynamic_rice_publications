#!/usr/bin/perl
=header
This script is design to run pileup and SNP calling using qsub, which speed up the process in RIL_SNP_MAQ.pl.
--ref:   reference sequence
--parents: parents file
--fastq: dir of fastq
=cut

use Getopt::Long;
use FindBin qw($Bin $Script);

my %opt;
GetOptions(\%opt,"ref:s","fastq:s","parent:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --ref ../input/reference/Pseudomolecules_Nipponbare.Build4.0.fasta --fastq ../input/testfastq\n";
   exit();
}

my $maq="/opt/tyler/bin/maq";
my @map=glob("$opt{fastq}/GN*.Maq.p1.map.pileup");

pileup(\@map);


sub pileup
{
my ($map)=@_;
my @cmd;
my $maq="/opt/tyler/bin/maq";
my $convert="/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/pileup2SNP.pl";
for(my $i=0; $i< @$map; $i++){
   #push @cmd, "$maq pileup -vP -q 40 $opt{ref}.bfa $map->[$i] | awk '\$4 > 0' > $map->[$i].pileup" unless (-e "$map->[$i].pileup");
   push @cmd, "perl $convert --pileup $map->[$i] --parent $opt{parent}" unless (-e "$map->[$i].SNP");
   open OUT, ">snp.sh" or die "$!";
   for(my $i=0; $i<@cmd; $i++){
      print OUT "$cmd[$i]\n";
   }
   close OUT;
}# for loop
#`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --maxjob 30 --resource nodes=1:ppn=1,mem=5G,walltime=100:00:00 snp.sh`;
`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --maxjob 20 --queue js --resource walltime=100:00:00 snp.sh`;
}# end of sub function


