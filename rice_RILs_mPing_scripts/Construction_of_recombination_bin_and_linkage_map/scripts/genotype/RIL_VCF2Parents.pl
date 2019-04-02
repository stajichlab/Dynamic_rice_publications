#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"vcf:s","project:s","help");


my $help=<<USAGE;
perl $0 --vcf --project
Convert vcf SNP file into Parents format, which used in MPR package or RIL_SNP_MAQ.pl to predict SNP of RILs
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "NB.RILs.dbSNP.SNPs";
convert($opt{vcf});

#Chr1    31071   .       A       G       4311.41 PASS    AC=8;AF=1.00;AN
sub convert
{
my ($file)=@_;
my %hash;
open OUT, ">$opt{project}.parents" or die "$!";
print OUT "V1\tV2\n";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $unit[0]=~s/\D+//ig;
    $unit[0]=sprintf("%02d",$unit[0]);
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[3];
    next if (exists $hash{$snpID});
    print OUT "$snpID\t$unit[3]\t$unit[4]\n";    
    $hash{$snpID}=1; 
}
close IN;
close OUT;
}
 
