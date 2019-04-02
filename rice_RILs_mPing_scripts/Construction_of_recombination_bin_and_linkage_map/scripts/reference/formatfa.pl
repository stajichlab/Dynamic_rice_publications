#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fa:s","project:s","help");


my $help=<<USAGE;
perl $0 --fa test.fa --project Nipponbare

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable("$opt{fa}");

sub readtable
{
my ($file)=@_;
my $out= $1.".reform.fa" if ($file=~/(.*)\.fa/);
open OUT, ">$out" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    if ($_=~/^>Chr(\d+)/){
       my $chr= length $1 > 1 ? $1 : sprintf("%02d",$1);
       print OUT ">chromosome$chr $opt{project}\n";
    }else{
       print OUT "$_\n";
    }
}
close IN;
close OUT;
}
 
