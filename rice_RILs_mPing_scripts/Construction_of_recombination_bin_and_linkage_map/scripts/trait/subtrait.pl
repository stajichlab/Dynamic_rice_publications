#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"trait:s","list:s","maqlist:s","help");


my $help=<<USAGE;
perl $0 --trait --list 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

if ($opt{list}){
my $list=readlist($opt{list});
readtable("$opt{trait}",$list);
}
if ($opt{maqlist}){
my $list=readmaqlist($opt{maqlist});
readtable("$opt{trait}",$list);
}

##GN-1	103	92	120.2	22	42.0
sub readtable
{
my ($file,$list)=@_;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    print "$_\n" if $_=~/^Sample/;
    my @unit=split("\t",$_);
    $unit[0]=~s/\-//;
    my $line=join("\t",@unit);
    if (exists $list->{$unit[0]}) {
       print "$line\n";
    }
}
close IN;
}

##Illumina/RIL10_0/RIL10_0_GAGTGG_p1.fq ./Illumina/RIL10_0/RIL10_0_GAGTGG_p2.fq GN10_1.fq       GN10_2.fq
sub readlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $name=$1 if ($unit[2]=~/(\w+)\_1\.fq/);
    $hash{$name}=1;
}
close IN;
return \%hash;
}

sub readmaqlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $name=$unit[0];
    $hash{$name}=1;
}
close IN;
return \%hash;
} 
