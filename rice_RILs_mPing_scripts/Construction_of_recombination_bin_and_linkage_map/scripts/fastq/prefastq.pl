#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"RIL:s","help");


my $help=<<USAGE;
perl $0 --RIL ./Illumina

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @fq=glob("$opt{RIL}/RIL*/*.fq");
my %rils;
`rm GN*.fq`;
for (my $i=0;$i<@fq;$i++){
   #print "$fq[$i]\n";
   if ($fq[$i]=~/(.*\/RIL(\d+)\_\d+\_*\w*)_p1\.fq/){
      my $name=$2;
      my $fq2 =$1."_p2.fq";
      my $file1=exists $rils{$name} ? "GN".$name."a"."_1".".fq" : "GN".$name."_1".".fq";
      my $file2=exists $rils{$name} ? "GN".$name."a"."_2".".fq" : "GN".$name."_2".".fq";
      print "$fq[$i]\t$fq2\t$file1\t$file2\n";
      #`ln -s $fq[$i] $file1`;
      #`ln -s $fq2 $file2`;
      `/rhome/cjinfeng/software/tools/seqtk-master/seqtk sample -s100 $fq[$i] 350000 > $file1`;
      `/rhome/cjinfeng/software/tools/seqtk-master/seqtk sample -s100 $fq2 350000 > $file2`;
      #my $file3=exists $rils{$name} ? "GN".$name."a"."_1"."sanger.fq" : "GN".$name."_1"."sanger.fq";
      #my $file4=exists $rils{$name} ? "GN".$name."a"."_2"."sanger.fq" : "GN".$name."_2"."sanger.fq";
      #`/rhome/cjinfeng/software/tools/slx2sanger/slx2sanger $file1 $file3`;
      #`/rhome/cjinfeng/software/tools/slx2sanger/slx2sanger $file2 $file4`;
      $rils{$name}=1;
   }
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
 
