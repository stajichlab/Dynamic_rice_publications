#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"trait:s","mping:s","cold:s","help");


my $help=<<USAGE;


USAGE


if ($opt{help}){
    print "$help\n";
    exit();
}

$opt{trait} ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt";
$opt{mping} ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Aug2_2013.mPing.table";
$opt{cold}  ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Aug2_2013.Cold.table";

my $refmping=readmping($opt{mping});
my $refcold =readcold($opt{cold});

readtable($opt{trait},$refmping,$refcold);

########
#Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)
#GN-1    103     92      120.2   22      42.0

sub readtable
{
my ($file,$mping,$cold)=@_;
my %hash;
open OUT, ">$file.ALL" or die "$!";
open IN, "$file" or die "$!";
my $head=<IN>;
chomp $head;
print OUT "$head\t$mping->{Title}->[0]\t$mping->{Title}->[1]\t$cold->{Title}->[0]\t$cold->{Title}->[1]\t$cold->{Title}->[2]\n";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $mpingline= exists $mping->{$unit[0]} ? "$mping->{$unit[0]}->[0]\t$mping->{$unit[0]}->[1]" : "NA\tNA";
    my $coldline = exists $cold->{$unit[0]} ? "$cold->{$unit[0]}->[0]\t$cold->{$unit[0]}->[1]\t$cold->{$unit[0]}->[2]" : "NA\tNA\tNA";
    print OUT "$_\t";
    print OUT "$mpingline\t";
    print OUT "$coldline\n";
}
close IN;
close OUT;
return \%hash;
}
 

#RILs    Non-ref mPing   Unique mPing
#RIL39_1 181     0
sub readmping
{
my ($file)=@_;
my %hash;
$hash{"Title"}=["Non-ref mPing","Unique mPing"];
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $ril="GN-".$1 if ($unit[0]=~/RIL(\d+)\_\d+/);
    $hash{$ril}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}

#RIL ID  Cold_ratio_index        CV(cold_ratio_index)    Rank Cold Ratio Index
#GN-1    0.92    12.50   243
sub readcold
{
my ($file)=@_;
my %hash;
$hash{"Title"}=["Cold_ratio_index","CV(cold_ratio_index)","Rank Cold Ratio Index"];
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2],$unit[3]];
}
close IN;
return \%hash;
}




