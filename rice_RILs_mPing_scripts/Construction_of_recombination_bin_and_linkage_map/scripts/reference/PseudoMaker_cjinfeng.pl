#perl Build_Pseudo_by_SNP.pl snp_file reference.fa parent_name
if (@ARGV < 3){
   print "perl $0 snp_file reference.fa parent2_name\n";
   exit();
}

open IN, $ARGV[0] or die "$!";
while (<IN>){
  if ($_=~/^Chr(\d+)\s+(\d+)\s+(\w+)\s+(\w+)$/){
    $chr=$1+1-1;$pos=$2;
    $v="$chr"."_"."$pos";
    $Ibase{$v}=$3;
  }
}
close (IN);
$chr=0;
open IN2, $ARGV[1] or die "$!";
open OUT, ">$ARGV[2].fa" or die "$!";
while ( defined ($a=<IN2>) ){
  if ($a=~/(chromosome\d+)/){
        print "$1\n";
  	$chr=$chr+1;$pos=0;print OUT '>';print OUT "$1\n";  	
  }
  else {
  	@temp=split //,$a;
  	foreach $temp(@temp){
  	if ($temp =~/\w/){
  		$pos=$pos+1;
  		$v="$chr"."_"."$pos";
  		$value=$Ibase{$v};
  		if (defined ($value)){
  		   $temp = $Ibase{$v};
  		   print OUT "$temp";
  		}
  		else {
  		   print OUT "$temp";
  		}
  	}
  	else {
  	 	print OUT "\n";
  	}
}
  }	
  	
}
close (IN2);close (OUT);
