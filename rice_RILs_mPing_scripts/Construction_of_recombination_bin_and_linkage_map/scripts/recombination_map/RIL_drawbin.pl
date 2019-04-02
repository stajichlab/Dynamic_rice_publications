#!/usr/bin/perl
use Getopt::Long;
use lib "/rhome/cjinfeng/BigData/software/lib";
use SVG;
use FontSize;
#use warnings;
#use strict;

my $font = FontSize->new();
my $font_family = "Arial";
my $font_size = 12;

GetOptions (\%opt,"MPR:s","chrlen:s","project:s","help");


my $help=<<USAGE;
perl $0 --MPR ./MPR_chr05 --chrlen MSU7.chr.inf
--MPR: MPR results directory
--chrlen: chromosome length file, chr01\t4300000\tcentstart\tcentend\n
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "MPR";

`mkdir $opt{MPR}/MPR_bin`;

my $width=800;
my $height=1400;
my $binh  =10;
my $chrlen  =readchrlen($opt{chrlen});
my @chrinf  =values %$chrlen;
my $maxlen  =max1(\@chrinf);
my $bp_per_pix=$maxlen/($width-150); #max chromosome length/max chromosome length in figure
my $refgeno=readgeno("$opt{MPR}/MPR.geno.data"); ##raw genotype data of RILs
my $refgenohmm=readgeno("$opt{MPR}/MPR.geno.data.HMMcr"); ##HMM corrected genotype data of RILs
my $refbin =readbin("$opt{MPR}/MPR.geno.bin");
my $refbinfill =readbin("$opt{MPR}/MPR.geno.bin.fill");
my $refbinuniq =readbin("$opt{MPR}/MPR.geno.bin.uniq");
my $refparents =readparents("$opt{MPR}/MPR.allele.MPR");
foreach my $rils (sort keys %$refgeno){
   print "$rils\n";
   $rils=~s/\"//g;
   #next;
   #next unless ($rils=~/GN15\d{1}$/ or $rils=~/MZ00.*/);
   #next unless ($rils=~/GN131/ or $rils=~/GN80/ or $rils=~/GN83/);
   next unless ($rils=~/GN131/);
   my $svg=SVG->new(width=>$width,height=>$height);
   my $xstart=100; my $ystart=0; my $count=0; 
   foreach my $chr0 (sort { $a <=> $b } keys %{$refgeno->{$rils}}){
      my $chr="Chr".$chr0;
      $ystart=110+$ystart;
      #next unless ($chr0 eq 1); ## only draw chr1
      print "$chr\t$xstart\t$ystart\n";
      my $len=$chrlen->{$chr}->[0];
      drawparents($svg,$xstart,$ystart,$binh,$refparents->{$chr0},$bp_per_pix);
      drawgeno($svg,$xstart,$ystart-($binh+5),$binh,$refgeno->{$rils}->{$chr0},$bp_per_pix);
      drawgenohmm($svg,$xstart,$ystart-2*($binh+5),$binh,$refgenohmm->{$rils}->{$chr0},$bp_per_pix);
      drawbin($svg,$xstart,$ystart-4*($binh+5),$binh,$refbin->{$rils}->{$chr0},$bp_per_pix);
      #drawbin($svg,$xstart,$ystart-5*($binh+5),$binh,$refbinfill->{$rils}->{$chr0},$bp_per_pix);
      drawbin($svg,$xstart,$ystart-5*($binh+5),$binh,$refbinuniq->{$rils}->{$chr0},$bp_per_pix); 
      drawcent($svg,$chrlen->{$chr}->[1],$chrlen->{$chr}->[2],$xstart,$ystart,$bp_per_pix);
      plot_ruler("svg",$svg,"Y",$ystart+10, "X_start",$xstart,"X_end",$xstart+$len/$bp_per_pix,"bp_start",0,"bp_end",$len,"scaletype","Mb","scaletypepos","right","scalestart","force","rulerstyle",2);
      drawtitle($svg,$xstart-50,$ystart,10,$chr);
      $count++;
   }
   drawsvg($svg,"$opt{project}.Recombination.Bin.$rils");
   `mv $opt{project}.Recombination.Bin.$rils.* $opt{MPR}/MPR_bin`;
}

if (0){
print "Draw chromosome bin:";
#draw chr bin
my $refgeno2=readgeno2("$opt{MPR}/MPR.geno.data");
my $refbin2=readbin2("$opt{MPR}/MPR.geno.bin");
my $refbinfill2 =readbin2("$opt{MPR}/MPR.geno.bin.fill");
my $refbinuniq2 =readbin2("$opt{MPR}/MPR.geno.bin.uniq");

drawchrbin($refbin2,$chrlen,$bp_per_pix,"original");
drawchrbin($refbinfill2,$chrlen,$bp_per_pix,"fill");
drawchrbin($refbinuniq2,$chrlen,$bp_per_pix,"uniq");
drawchrSNP($refgeno2,$refparents,$chrlen,$bp_per_pix,"geno");
}

###
sub drawtitle
{
my ($svg,$xstart,$ystart,$size,$t)=@_;
my $tag=$svg->text(
        x=> $xstart, y=> $ystart,
        style =>{
            'font-size' => $size
        }
   )->cdata($t);
} 


####
sub drawchrbin
{
my ($refbin2,$chrlen,$bp_per_pix,$title)=@_;
my $width2=800;
my $height2=700;
my $binh2   =2;
foreach my $chr0 (sort keys %$refbin2){
   my $chr="Chr".$chr0;
   my $len=$chrlen->{$chr}->[0];
   my $svg=SVG->new(width=>$width2,height=>$height2);
   my $xstart=100; my $ystart=0; my $count=0;
   plot_ruler("svg",$svg,"Y",$ystart+10, "X_start",$xstart,"X_end",$xstart+$len/$bp_per_pix,"bp_start",0,"bp_end",$len,"scaletype","Mb","scaletypepos", "right","scalestart","force","rulerstyle",2);
   print "$chr\tDraw this chromosome\n";
   #next; #### just test for input data
   my $max=0;
   my %linepos;##record bin line position
   foreach my $rils (sort keys %{$refbin2->{$chr0}}){
      #print "$rils\tAdd this RILs\n";
      my $yrils=$count*2+$ystart+30;
      $max = $yrils > $max ? $yrils : $max;
      drawbin($svg,$xstart,$yrils,$binh2,$refbin2->{$chr0}->{$rils},$bp_per_pix,\%linepos);
      #drawbinline($svg,$xstart,$ystart+20,$yrils+10,$refbin2->{$chr0}->{$rils},$bp_per_pix);
      drawtitle($svg,$xstart-10,$yrils+1.5,2,$rils);
      $count++;
   }
   drawbinline($svg,$xstart,$ystart+27,$max+5,\%linepos,$bp_per_pix);
   drawsvg($svg,"$opt{project}.Recombination.Bin.$chr.$title");
   `mv $opt{project}.Recombination.Bin.$chr.$title.* $opt{MPR}/MPR_bin/`;
}
}

####
sub drawchrSNP
{
my ($refgeno,$refparents,$chrlen,$bp_per_pix,$title)=@_;
my $width2=800;
my $height2=700;
my $binh2   =2;
foreach my $chr0 (sort keys %$refgeno){
   my $chr="Chr".$chr0;
   my $len=$chrlen->{$chr}->[0];
   my $svg=SVG->new(width=>$width2,height=>$height2);
   my $xstart=100; my $ystart=0; my $count=0;
   plot_ruler("svg",$svg,"Y",$ystart+10, "X_start",$xstart,"X_end",$xstart+$len/$bp_per_pix,"bp_start",0,"bp_end",$len,"scaletype","Mb","scaletypepos", "right","scalestart","force","rulerstyle",2);
   print "$chr\tDraw this chromosome\n";
   drawparents($svg,$xstart,$ystart+40,$binh2,$refparents->{$chr0},$bp_per_pix);
   drawcent($svg,$chrlen->{$chr}->[1],$chrlen->{$chr}->[2],$xstart,$ystart,$bp_per_pix); 
   foreach my $rils (sort keys %{$refgeno->{$chr0}}){
      print "$rils\tAdd this RILs\n";
      my $yrils=$count*2+$ystart+45;
      drawgeno($svg,$xstart,$yrils,$binh2,$refgeno->{$chr0}->{$rils},$bp_per_pix); 
      drawtitle($svg,$xstart-10,$yrils-0.5,2,$rils);
      $count++;
   } 
   drawsvg($svg,"$opt{project}.Recombination.Bin.$chr.$title");
   `mv $opt{project}.Recombination.Bin.$chr.$title.* $opt{MPR}/MPR_bin/`;
}
}




#####################################
sub readchrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split(/\s+/,$_);
    $hash{$unit[0]}=[$unit[1],$unit[2],$unit[3]];
}
close IN;
return \%hash;
}

sub drawcent
{
my ($svg,$cents,$cente,$xstart,$ystart,$bp_per_pix)=@_;
my $x=$xstart+$cents/$bp_per_pix;
my $y=$ystart+1;
my $width=($cente-$cents+1)/$bp_per_pix;
my $height=9;
my $tag=$svg->rectangle(
       x=>$x,y=>$y,
       width=>$width,height=>$height,
       style=>{
           'fill'=>"black"  
       }       
   );
}

sub max1
{
my ($num)=@_;
my $max=0;
foreach  (@$num) {
         next if ($_->[0] eq "NA");
         $max = $_->[0] > $max ? $_->[0] : $max;
}
return $max;
}


sub max
{
my ($num)=@_;
my $max=0;
foreach  (@$num) {
        next if ($_ eq "NA");
        $max = $_ > $max ? $_ : $max;
}
return $max;
}


###read MPR.geno.data or MPR.geno.data.HMMcr
###"0501941871G"   NA      NA      NA      NA      NA      0       NA      NA      1       1
sub readgeno
{
my ($file)=@_;
my %hash;
open IN,"$file" or die "$!";
my $header=<IN>;
chomp $header;
$header=~s/\"//g;
my @head=split("\t",$header);
#print join("\t",@head),"\n";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    next if ($unit[0] eq "NA");
    my $chr =substr($unit[0],0,2);
    my $pos =substr($unit[0],2,8);
    #print "$chr\t$pos\n";
    for(my $i=1; $i<@unit; $i++){
       push @{$hash{$head[$i]}{$chr}},[$pos,$unit[$i]]; 
    }
}
close IN;
return \%hash;
}

###
###read MPR.geno.bin or MPR.geno.bin.filled or MPR.geno.bin.uniq
###"0501941871G"   NA      NA      NA      NA      NA      0       NA      NA      1       1
sub readbin
{
my ($file)=@_;
my %hash;
open IN,"$file" or die "$!";
my $header=<IN>;
chomp $header;
$header=~s/\"//g;
my @head=split("\t",$header);
#print join("\t",@head),"\n";
my %last;
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    my $chr =substr($unit[0],0,2);
    my $pos =substr($unit[0],2,8);
    my $temp=exists $last{$chr} ? $last{$chr} : -1;
    my $start=$temp+1;
    $last{$chr}=$pos;
    #print "$chr\t$pos\n";
    for(my $i=1; $i<@unit; $i++){
       push @{$hash{$head[$i]}{$chr}},[$start,$pos,$unit[$i]]; 
    }
}
close IN;
return \%hash;
}

###read MPR.geno.data or MPR.geno.data.HMMcr
###"0501941871G"   NA      NA      NA      NA      NA      0       NA      NA      1       1
sub readgeno2
{
my ($file)=@_;
my %hash;
open IN,"$file" or die "$!";
my $header=<IN>;
chomp $header;
$header=~s/\"//g;
my @head=split("\t",$header);
#print join("\t",@head),"\n";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    next if ($unit[0] eq "NA");
    my $chr =substr($unit[0],0,2);
    my $pos =substr($unit[0],2,8);
    #print "$chr\t$pos\n";
    for(my $i=1; $i<@unit; $i++){
       push @{$hash{$chr}{$head[$i]}},[$pos,$unit[$i]]; 
    }
}
close IN;
return \%hash;
}




###
###read "0501941871G" "A" "G"
sub readparents
{
my ($file)=@_;
my %hash;
open IN,"$file" or die "$!";
my $header=<IN>;
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    my $chr =substr($unit[0],0,2);
    my $pos =substr($unit[0],2,8);
    #print "$chr\t$pos\n";
    push @{$hash{$chr}}, $pos; 
}
close IN;
return \%hash;
}



############
sub drawparents
{
my ($svg,$xstart,$ystart,$height,$pos,$bp_per_pix)=@_;
for(my $i=0;$i<@$pos;$i++){
   next if ($pos->[$i] eq "NA");
   my $x=$xstart+$pos->[$i]/$bp_per_pix;
   my $y1=$ystart-$height;
   my $y2=$ystart;
   my $color="black";
   
   my $cov=$svg->line(
      x1=>$x,y1=>$y1,
      x2=>$x,y2=>$y2,
      style=>{
           'stroke'=>$color,
           'stroke-width' => 0.1
      }
   );
}
}





sub drawgeno
{
my ($svg,$xstart,$ystart,$height,$pos,$bp_per_pix)=@_;
my $count=0;
for(my $i=0;$i<@$pos;$i++){
   next if ($pos->[$i]->[1] eq "NA");
   $count++;
   my $x=$xstart+$pos->[$i]->[0]/$bp_per_pix;
   my $y1=$ystart-$height/2;
   my $y2=$pos->[$i]->[1] > 0 ? $ystart : $ystart-$height;
   my $color=$pos->[$i]->[1] > 0 ? "blue" : "red";
   
   my $cov=$svg->line(
      x1=>$x,y1=>$y1,
      x2=>$x,y2=>$y2,
      style=>{
           'stroke'=>$color,
           'stroke-width' => 0.1
      }
   );
}
print "Draw $count points\n";
}

##
###
###read MPR.geno.bin or MPR.geno.bin.filled or MPR.geno.bin.uniq
###"0501941871G"   NA      NA      NA      NA      NA      0       NA      NA      1       1
sub readbin2
{
my ($file)=@_;
my %hash;
open IN,"$file" or die "$!";
my $header=<IN>;
chomp $header;
$header=~s/\"//g;
my @head=split("\t",$header);
#print join("\t",@head),"\n";
my %last;
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    my $chr =substr($unit[0],0,2);
    my $pos =substr($unit[0],2,8);
    my $temp=exists $last{$chr} ? $last{$chr} : -1;
    my $start=$temp+1;
    $last{$chr}=$pos;
    #print "$chr\t$pos\n";
    for(my $i=1; $i<@unit; $i++){
       push @{$hash{$chr}{$head[$i]}},[$start,$pos,$unit[$i],$unit[0]]; 
    }
}
close IN;
return \%hash;
}


###
sub drawgenohmm
{
my ($svg,$xstart,$ystart,$height,$pos,$bp_per_pix)=@_;
for(my $i=0;$i<@$pos;$i++){
   next if ($pos->[$i]->[1] eq "NA");
   my $x=$xstart+$pos->[$i]->[0]/$bp_per_pix;
   my $y1=$ystart;
   my $y2=$ystart-$height;
   my $color;
   if ($pos->[$i]->[1] == 0){
      $color="red";
   }
   if ($pos->[$i]->[1] == 1){
      $color="blue";
   }
   if ($pos->[$i]->[1] == 0.5){
      $color="green";
   }
   my $cov=$svg->line(
      x1=>$x,y1=>$y1,
      x2=>$x,y2=>$y2,
      style=>{
           'stroke'=>$color,
           'stroke-width' => 0.1
      }
   );
}
}

###
sub drawbin
{
my ($svg,$xstart,$ystart,$height,$pos,$bp_per_pix,$linepos)=@_;
for(my $i=0;$i<@$pos;$i++){
   next if ($pos->[$i]->[2] eq "NA");
   my $x=$xstart+$pos->[$i]->[0]/$bp_per_pix;
   my $y=$ystart;
   my $h=$height;
   my $w=($pos->[$i]->[1]-$pos->[$i]->[0]+1)/$bp_per_pix;
   my $temp=$xstart + $pos->[$i]->[1]/$bp_per_pix;
   $linepos->{$temp}=$pos->[$i]->[3];
   #print "$x\t$y\t$w\t$h\n";
   my $color=$pos->[$i]->[2] > 0 ? "blue" : "red"; 
   my $rec=$svg->rectangle(
      x=>$x,y=>$y,
      width=>$w,height=>$h,
      style=>{
           'fill'=>$color,
           #'stroke' => "black",
           'stroke-width' => 0.1
      }
   );
}
}

##
###
sub drawbinline
{
my ($svg,$xstart,$ystart,$yend,$line,$bp_per_pix)=@_;
my @pos;
push @pos, keys %$line;
for(my $i=0;$i<@pos;$i++){
   my $x=$pos[$i];
   my $y1=$ystart;
   my $y2=$yend;
   my $l=$svg->line(
       x1=>$x,y1=>$y1,
       x2=>$x,y2=>$y2,
       style=>{
           'stroke'=>"gray",
           'stroke-width' => 0.1,
           'stroke-dashoffset' => 1
       }
   );
   my $name=$line->{$pos[$i]};
   my $x1=$x+0.5;
   my $y=$y2+15;
   my $tag=$svg->text(
         x=> $x1, y=> $y,'transform' => "rotate(-90,$x1, $y)",
         style =>{
              'font-size' => 2,
              'text-anchor' => 'start',
         }
         #transform => 'rotate(-90)'
      )->cdata($name);
}
}




#######
sub plot_ruler()
{
        my %rulcfg = @_;
        $rulcfg{scaletype} ||= "scale";
        $rulcfg{scaletypepos} ||= "left";
        $rulcfg{scalestart} ||= "auto";
        $rulcfg{rulerstyle} ||= "3";


        my $scale_size = 6;
        my $divid = 10; #my $divid = 50
        my $unit_start;

        my $bp_len = $rulcfg{bp_end} - $rulcfg{bp_start};
        my $X_len = $rulcfg{X_end} - $rulcfg{X_start};

        my ($str,$str1,$str2,$unit);
        $str = $bp_len / $divid;
        $str = sprintf("%e",$str);
        if ($str=~/([0-9\.]+)e([0-9+-]+)/) {
                $str1 = $1;
                $str2 = $2;
        }
        $str1 = int ( $str1 + 0.5 );
        $unit = $str1 * 10 ** $str2;
        $unit = $rulcfg{bigscalesize}/10 if(defined $rulcfg{bigscalesize});

        my $g = $rulcfg{svg}->group('id'=>times().rand());

        ## draw the main axis
        $g->line('x1',$rulcfg{X_start},'y1',$rulcfg{Y},'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1);
        return if($rulcfg{bp_end} - $rulcfg{bp_start}  == 0);
        
        ##add end mark, cjinfeng
        $g->line('x1',$rulcfg{X_end},'y1',$rulcfg{Y} - $scale_size,'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
        $g->line('x1',$rulcfg{X_end},'y1',$rulcfg{Y} + $scale_size,'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '2');
        $g->line('x1',$rulcfg{X_end},'y1',$rulcfg{Y} - $scale_size,'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
        ##draw ruler mark text at the specified postion(left or right of the ruler)
        if ($rulcfg{scaletypepos} eq "left") {
                $g->text('x',$rulcfg{X_start}-$font->stringWidth($font_family,$font_size,$rulcfg{scaletype})-6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$font_family,"font-size",$font_size,"fill",'#000000');
        }
        if ($rulcfg{scaletypepos} eq "right") {
                $g->text('x',$rulcfg{X_end} + 6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$font_family,"font-size",$font_size,"fill",'#000000');
        }

        ##decide unit start
        if ($rulcfg{scalestart} eq "auto") {
                $unit_start = $rulcfg{bp_start} + ($unit - $rulcfg{bp_start} % $unit);
        }
        if ($rulcfg{scalestart} eq "force") {
                $unit_start = int($rulcfg{bp_start} / 10 + 0.5) * 10; ##....................?	
        }
        ## draw small scale lines
        for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit) {
                my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
                $g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '2');
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
                my $scale_num = int ($i / 1000000);
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$font->stringHeight($font_size),'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '1');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$scale_size+$font->stringHeight($font_size)+2,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '2');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}-$scale_size-$font->stringHeight($font_size)+6,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '3');
        }
=pod
        ## draw big scale lines and text scales
        if((int($unit_start/1000))%10<5 && (int($unit_start/1000))%10>0)
        {
                                my $left=$unit_start%1000;
                                my $big=int$unit_start/10000;
                                $unit_start=$big*10000+5000+$left;
        }
        elsif((int($unit_start/1000))%10>5)
        {
                        my $left=$unit_start%1000;
                        my $big=int$unit_start/10000;
                        $unit_start=($big+1)*10000+$left;
        }
        else{$unit_start=$unit_start;}
        for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit*10) {
                my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
                $g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1)  if ($rulcfg{rulerstyle} eq '2');
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
                my $scale_num = int ($i / 1000000);
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$font->stringHeight($font_size),'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '1');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$scale_size+$font->stringHeight($font_size)+2,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '2');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}-$scale_size-$font->stringHeight($font_size)+6,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '3');
        }
=cut
}




sub drawsvg
{
my ($svg,$name)=@_;
open OUT,">$name\.svg";
print OUT $svg->xmlify;
close OUT;
system("/rhome/cjinfeng/BigData/software/draw/svg2xxx_release/svg2xxx -t pdf -m 10000 $name.svg > step02.recombination_bin.sh.svg.draw 2> step02.recombination_bin.sh.svg.draw2");
}
 
