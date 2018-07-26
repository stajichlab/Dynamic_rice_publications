#!/usr/bin/perl

# Program to convert output of VCFtools' vcf-to-tab
# to FASTA alignment.

# Sample input file
#	$ head input.vcf.tab
#	chr10	94051	C	./	./	./	./	./	T/T
#	chr10	94056	T	./	./	./	./	./	C/C
#	chr10	94180	G	./	A/A	./	./	./	./

use strict;
use warnings;
use Getopt::Long;

my $exclude_het = 0;	# Should we exclude heterozygous SNPs?
my $output_ref  = 0;	# Should we output the reference calls?
my $input_tab;

my $usage = "usage: $0 [--exclude_het] [--output_ref] -i input.tab";

my $result = GetOptions (	"exclude_het"	=> \$exclude_het,
							"output_ref"	=> \$output_ref,
							"i=s"			=> \$input_tab
						) or die "Incorrect usage. $usage\n";

my $starting_col = 3;
if ($output_ref) {
	print STDERR "Including reference sequence. Remove --output_ref flag to exclude.\n";
	$starting_col = 2;
}

my %iupac = (
			'G/G' => 'G',
			'C/C' => 'C',
			'T/T' => 'T',
			'A/A' => 'A',

			'G/T' => 'K',
			'T/G' => 'K',
			'A/C' => 'M',
			'C/A' => 'M',
			'C/G' => 'S',
			'G/C' => 'S',
			'A/G' => 'R',
			'G/A' => 'R',
			'A/T' => 'W',
			'T/A' => 'W',
			'C/T' => 'Y',
			'T/C' => 'Y',

			'./.' => '.',
		);

open (TAB, "<$input_tab") or die "ERROR: Could not open input file $input_tab.\n";
my $header = <TAB>;
my @col_names = split /\t/, $header;


# Get number of columns
my $num_cols = scalar @col_names;
print STDERR "Number of columns:\t$num_cols\n";

my %hash;
my $count = 0;
my $pos = 0;
LINE: foreach my $line (<TAB>) {
	my @data = split /\t/, $line;
	$pos++;
	print "number of position parsed: ï$pos\n";
	# Skip if this is indel (Length of @data will be less than $num_cols)
	if ((scalar @data) < $num_cols) {
		print STDERR "Skipping indel.\n";
		next LINE;
	}
	
	for (my $i = $starting_col; $i < $num_cols; $i++) {
		my $bp = $data[$i]; 
		chomp $bp;
                # Skip if any basepairs are actually 2 or more together
		if ($bp =~ /\w{2,}/) {
			print STDERR "Skipping multi-basepair insertion.\n";
			next LINE;
		}
                # Exclude heterozygotes. Keep only fixed SNPs
                if ($exclude_het) {
                	if ($bp =~ /(\w)\/(\w)/) {
                        	if ($1 ne $2) {
        				print STDERR "Skipping heterozygote. ";
					print STDERR "Remove --exclude_het flag to retain.\n";
					next LINE;
                                }
                        }
                }        
	}

	# Otherwise write line to pure temporary file
	#print TEMP $line;
	$count++;
	for (my $i = $starting_col; $i < $num_cols; $i++) {
		my $nuc = $data[$i];
		chomp $nuc;
		#If we're reference, just print basepair
		if ($i == 2){
			$hash{$i} .= $nuc;
		# Haploid
		}elsif ($nuc =~ /(\w)\/$/) {
			$hash{$i} .= $1;
		# Missing data
		} elsif ($nuc eq './' or $nuc eq './.') {
			$hash{$i} .= '-';
		# Data
		} elsif ($nuc =~ /(\w)\/(\w)/) {
			my $first = $1;
			my $second = $2;
			# Homozygote
			if ($first eq $second) {
				$hash{$i} .= $first;
			# Heterozygote
			} else {
				my $gt = $first . '/' . $second;
				if ( !exists($iupac{$gt}) ) { die "ERROR: BP is $nuc\n"; }
                		$gt = $iupac{$gt};
				$hash{$i} .= $gt;
			}
		}
		if ($count == 100) {
			$count = 0;
			$hash{$i} .= "\n";
		}
	}
}
	
close TAB;

my $outfile = "$input_tab.fasta";
open (OUT, ">$outfile") or die "ERROR: Could not open output fasta file $outfile.\n";
foreach my $i (sort {$a <=> $b} keys %hash){
	my $ind = $col_names[$i];
	chomp $ind;
	print OUT ">$ind\n";
	print OUT "$hash{$i}\n";
}
close OUT;
exit;
