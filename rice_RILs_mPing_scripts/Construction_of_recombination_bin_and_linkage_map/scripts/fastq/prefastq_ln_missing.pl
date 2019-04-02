#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"RIL:s","fastq:s","read:s","help");


my $help=<<USAGE;
perl $0 --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{read} ||= 950000; ##950000 reads will make 0.5X of 372 Mb of rice genome. 950000x100x2/372000000 ~= 0.5
my @fq=glob("$opt{RIL}/RIL*/*.fq");
my %rils;
#open OUT, ">subfastq.sh" or die "$!";
for (my $i=0;$i<@fq;$i++){
   #print "$fq[$i]\n";
   if ($fq[$i]=~/(.*\/RIL(\d+)\_.*)_p1\.fq/){
      my $name=$2;
      print "$name\n";
      my $fq2 =$1."_p2.fq";
      next if ($rils{$name}); ## skip if already have fastq record (these line have multi lane sequenced)
      print "CK: $name\n";
      my $file1=exists $rils{$name} ? "$opt{fastq}/"."GN".$name."a"."_1".".fq" : "$opt{fastq}/"."GN".$name."_1".".fq";
      my $file2=exists $rils{$name} ? "$opt{fastq}/"."GN".$name."a"."_2".".fq" : "$opt{fastq}/"."GN".$name."_2".".fq";
      print "$file1\n$file2\n";
      unless (-e "$file1" and -e "$file2"){
         print "$fq[$i]\t$fq2\t$file1\t$file2\n";
         #my $cmd1="/rhome/cjinfeng/software/tools/seqtk-master/seqtk sample -s100 $fq[$i] $opt{read} > $file1";
         `ln -s $fq[$i] $file1`;
         #my $cmd2="/rhome/cjinfeng/software/tools/seqtk-master/seqtk sample -s100 $fq2 $opt{read} > $file2";
         `ln -s $fq2 $file2`;
         #print OUT "$cmd1\n$cmd2\n";
      }
      $rils{$name}=1;
   }
}
#close OUT;

#`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --resource nodes=1:ppn=1,mem=2G,walltime=100:00:00 subfastq.sh`;

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
 
