#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"trait:s","mping:s","ping:s","cold:s","usda:s","seed:s","help");


my $help=<<USAGE;


USAGE


if ($opt{help}){
    print "$help\n";
    exit();
}

$opt{trait} ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt";
$opt{mping} ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Sep20_2017.mPing.table";
$opt{cold}  ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Aug2_2013.Cold.table";
$opt{ping}  ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Sep20_2017.Ping.table";
$opt{usda}  ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/Jan30_2014.USDA.table";
$opt{seed}  ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/New04_2015.SeedDimension_Salt.table";

my $refmping=readmping($opt{mping});
my $refcold =readcold($opt{cold});
my $refping =readping($opt{ping});
my $refusda =readusda($opt{usda});
my $refseed =readseed($opt{seed});

readtable($opt{trait},$refmping,$refcold,$refping,$refusda,$refseed);

########
#Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)
#GN-1    103     92      120.2   22      42.0

sub readtable
{
my ($file,$mping,$cold,$ping,$usda,$seed)=@_;
my %hash;
open OUT, ">$file.ALL" or die "$!";
open IN, "$file" or die "$!";
my $head=<IN>;
chomp $head;
print OUT "$head\t$mping->{Title}->[0]\t$mping->{Title}->[1]\t$ping->{Title}->[0]\t$cold->{Title}->[0]\t$cold->{Title}->[1]\t$cold->{Title}->[2]\t$refusda->{Title}->[0]\t$refusda->{Title}->[1]\t$refusda->{Title}->[2]\t$refusda->{Title}->[3]\t$refusda->{Title}->[4]\t$refusda->{Title}->[5]\t$refusda->{Title}->[6]\t$refseed->{Title}->[0]\t$refseed->{Title}->[1]\t$refseed->{Title}->[2]\t$refseed->{Title}->[3]\t$refseed->{Title}->[4]\t$refseed->{Title}->[5]\n";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $mpingline= exists $mping->{$unit[0]} ? "$mping->{$unit[0]}->[0]\t$mping->{$unit[0]}->[1]" : "NA\tNA";
    my $pingline = exists $ping->{$unit[0]} ? "$ping->{$unit[0]}->[0]" : "NA";
    my $coldline = exists $cold->{$unit[0]} ? "$cold->{$unit[0]}->[0]\t$cold->{$unit[0]}->[1]\t$cold->{$unit[0]}->[2]" : "NA\tNA\tNA";
    my $usdaline = exists $usda->{$unit[0]} ? "$usda->{$unit[0]}->[0]\t$usda->{$unit[0]}->[1]\t$usda->{$unit[0]}->[2]\t$usda->{$unit[0]}->[3]\t$usda->{$unit[0]}->[4]\t$usda->{$unit[0]}->[5]\t$usda->{$unit[0]}->[6]" : "NA\tNA\tNA\tNA\tNA\tNA\tNA";
    my $seedline = exists $seed->{$unit[0]} ? "$seed->{$unit[0]}->[0]\t$seed->{$unit[0]}->[1]\t$seed->{$unit[0]}->[2]\t$seed->{$unit[0]}->[3]\t$seed->{$unit[0]}->[4]\t$seed->{$unit[0]}->[5]" : "NA\tNA\tNA\tNA\tNA\tNA";
    print OUT "$_\t";
    print OUT "$mpingline\t";
    print OUT "$pingline\t";
    print OUT "$coldline\t";
    print OUT "$usdaline\t";
    print OUT "$seedline\n";
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
    my $ril="GN-".$1 if ($unit[0]=~/RIL(\d+)/);
    $hash{$ril}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}

#RIL	ping
#RIL1	1
sub readping
{
my ($file)=@_;
my %hash;
$hash{"Title"}=["Ping"];
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $ril="GN-".$1 if ($unit[0]=~/RIL(\d+)$/);
    $hash{$ril}=[$unit[1]];
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



#Type    Name    Days to Heading_TX      Height_TX       Height Single Plants    Grain yield per tiller  Harvest Index
#RIL     1       94      86      95      1.88    0.299
#RIL     2       91      81      100     3.03    0.259
#RIL     3       87      89      97      2.68    0.413
#RIL     4       90      99      112     3.06    0.384
sub readusda
{
my ($file)=@_;
my %hash;
$hash{"Title"}=["Days_to_Heading_TX", "Height_TX", "Height_Single_Plants", "Grain_yield_per_tiller", "Harvest_Index", "seedling_12C", "seedling_26C"];
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print join("\t",@unit), "END\n";
    my $ril ='GN-'.$unit[1];
    for (my $i=2;$i<@unit;$i++){
        #print "CK ",$unit[$i], "\n";
        $unit[$i] = 'NA' if ($unit[$i]!~/\d+/);
    }
    #print "RIL: ", $ril, "\n";
    $hash{$ril}=[$unit[2],$unit[3],$unit[4],$unit[5],$unit[6],$unit[7],$unit[8]];
    #print "$unit[2],$unit[3],$unit[4],$unit[5],$unit[6],$unit[7],$unit[8],End\n";
}
close IN;
return \%hash;
}


#Type    RIL     Salt16  Salt10  HT10    SeedLength      Seedwidth       L_W_Ratio
#RIL     1       7.13    4.20    181.27  7.214   3.346   2.16
sub readseed
{
my ($file)=@_;
my %hash;
$hash{"Title"}=["Salt16", "Salt10", "Height_Salt10", "AvgStraightLength", "AvgStraightWidth", "LengthWidthRatio"];
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $ril = "GN-$unit[1]";
    $hash{$ril}=[$unit[2],$unit[3],$unit[4],$unit[5],$unit[6],$unit[7]];
    #print "$unit[5],$unit[6],$unit[7]\n";
}
close IN;
return \%hash;
}






