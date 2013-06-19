#!/usr/bin/env perl

$lagandir = $ENV{LAGAN_DIR};

# Status
#   -- extension problems

if (@ARGV < 2) {
    print ("usage:\n rechaos seqfile1 seqfile2 [-chaos \"chaos flags\"] [-recurse \"(wl1,nd1,co1),(wl2,nd2,co2),...\"] [-out \"filename\"] [-lazy] [-maskedonly] [-debug] [-translate] [-fastreject]\n");
    exit(1);
}

#$recurfl = "(12,0,25,0)x,(13,1,30,0)x,(8,1,30,0)x,(7,1,30,0)x";
$recurfl = "(12,0,25,0)x,(13,1,30,0)x,(4,0,4,3000)xt,(8,1,30,0)x,(7,1,30,0)x";
#$recurfl = "(12,0,10,200)x,(12,0,10,150)x,(3,0,10,150)xt,(8,0,10,150)x,(12,0,25,0),(13,1,30,0),(3,0,30,0)t,(8,1,30,0),(7,1,25,0)";
$minbox = 10;
$minside = 5;
$seq1 = $ARGV[0];
$seq2 = $ARGV[1];
$tofile = 0;
$masker = 1;
$lazycheck = 0;
$fastreject = 0;
$frminlevel = 0;
$frmaxlevel = 3;
@frseq1 = (150000, 50000, 30000, 15000);
@frseq2 = (150000, 50000, 30000, 15000);
#@frseq1 = (70000, 60000, 60000, 20000);
#@frseq2 = (70000, 60000, 60000, 20000);
$sentinelleft = 1.1;
$sentinelright = 1.2;
$gfc = " ";
$dounmasked = 1;
$filename = "";
$debug = 0;
$anchparams = "";
$translate = 0;

sub max {
    my ($a, $b) = @_;
    return $a if ($a > $b);
    return $b;    
}

sub min {
    my ($a, $b) = @_;
    return $a if ($a < $b);
    return $b;    
}

$i = 2;
while ($i < @ARGV) {
    if ($ARGV[$i] =~ /-\chaos/) {
	$chaosfl = $chaosfl." ".$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-ext/) {
	$chaosfl = $chaosfl." -ext ";
    }
    elsif ($ARGV[$i] =~ /-recurse/) {
	$recurfl = $ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-lazy/) {
	$lazycheck = 1;
    }
    elsif ($ARGV[$i] =~ /-nomask/) {
	$masker = 0;
    }
    elsif ($ARGV[$i] =~ /-out/) {
	$tofile = 1;
	$filename = $ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-maskedonly/) {
	$dounmasked = 0;
    }
    elsif ($ARGV[$i] =~ /-fastreject/) {
	$fastreject = 1;
    }
    elsif ($ARGV[$i] =~ /-debug/) {
	$debug = 1;
    }
    elsif ($ARGV[$i] =~ /-translate/) {
	$translate = 1;
    }
    elsif ($ARGV[$i] =~ /-gfc/) {
	$gfc = " -gfc ";
    }
    elsif ($ARGV[$i] =~ /-gap/){
	$anchparams = $anchparams." -gap ".$ARGV[++$i];
	$anchparams = $anchparams." ".$ARGV[++$i];
    }
    else { 
	die ("Unrecognized option $ARGV[$i]\n");
    }
    $i++;
}

if ($lazycheck) {
    if (-f $filename) {
	print STDERR "Output file already exists, lazy mode exit!\n";
	exit (0);
    }
}

$extracase1 = 0;
$extracase2 = 0;
if (-e "$seq1.masked") { $extra1 = $seq1; $seq1 = "$seq1.masked"; $extracase1 = 1; }
if (-e "$seq2.masked") { $extra2 = $seq2; $seq2 = "$seq2.masked"; $extracase2 = 1; }
if (! $dounmasked){ $extracase1 = 0; $extracase2 = 0; }

#open(SEQ1, "$seq1");
#open(SEQ2, "$seq2");

#$line1 = <SEQ1>;
#while ($line1 = <SEQ1>) {
#    chomp $line1;
#    $seq1len += length($line1);
#}
#
#$line2 = <SEQ2>;
#while ($line2 = <SEQ2>) {
#    chomp $line2;
#    $seq2len += length($line2);
#}

$seq1len = `$lagandir/utils/getlength $seq1`; chomp $seq1len;
$seq2len = `$lagandir/utils/getlength $seq2`; chomp $seq2len;

$b1[0] = $b2[0] = 1;
$e1[0] = $seq1len;
$e2[0] = $seq2len;

$cumanchs = 0;

$clipleft1 = 0;
$clipleft2 = 0;
$clipright1 = $seq1len + 1;
$clipright2 = $seq2len + 1;
$app_str = "";

$i = 0;
while (1) {
    $goodanchs = 0;
    $totalanchs = 0;
    
    $stillmore = ($recurfl =~ /\((\d+)\,(\d+)\,(\d+)\,(\d+)\)(\w*)(.*)/);
    if (! $stillmore) {
	if ($extracase1 || $extracase2) {
	    if ($extracase1) { $seq1 = $extra1; $extracase1 = 0; }
	    if ($extracase2) { $seq2 = $extra2; $extracase2 = 0; }
	}
	else {
	    last;
	}
    }
    else {
	$wordlen = $1;
	$degeneracy = $2;
	$cutoff = $3;
	$extcutoff = $4;
	$tail = $5;
	
	$extraparams = "";
	$extraparams = "-t ".$extraparams if ((index ($tail, "t") != -1) && ($translate));
	$extraparams = $extraparams." -rsc $extcutoff" if (index ($tail, "x") != -1);
    }

    $recurfl = $6;
    next if ((index ($tail, "t") != -1) && (!$translate));

    print STDERR "Using $seq1 $seq2 ($wordlen, $degeneracy, $cutoff, $extcutoff) $tail\n";

# PRINT OUT LIST OF REGIONS TO ALIGN

    open (PFILE, ">$$.anchs.pairs");
    for ($j = 0; $j < @b1; $j++) {
	print PFILE "-s1 $b1[$j] $e1[$j] -s2 $b2[$j] $e2[$j]\n";
    }
    close (PFILE);

#    print STDERR "PAIRS hits\n";
#    print STDERR `cat $$.anchs.pairs`;
#    print STDERR "-----------------\n";
#    print STDERR `cat $$.anchs.pairs`;
#    print STDERR "-----------------\n";
#    print STDERR "$lagandir/chaos $seq1 $seq2 -wl $wordlen -nd $degeneracy -co $cutoff $extraparams $gfc $chaosfl -pairs $$.anchs.pairs > $$.anchtemp";

# PERFORM THE ALIGNMENTS USING CHAOS
    
    $saver = "$lagandir/chaos $seq1 $seq2 $extraparams -wl $wordlen -nd $degeneracy -co $cutoff $gfc $chaosfl -pairs $$.anchs.pairs > $$.anchtemp";
    `$lagandir/chaos $seq1 $seq2 $extraparams -wl $wordlen -nd $degeneracy -co $cutoff $gfc $chaosfl -pairs $$.anchs.pairs > $$.anchtemp`;
    if ($?) { 
	print STDERR "$saver\n";
	exit(1); 
    }

# ADD IN BOUNDARIES

    $stillmore = ($recurfl =~ /\((\d+)\,(\d+)\,(\d+)\,(\d+)\)(\w*)(.*)/);
    if ($fastreject || $stillmore || $extracase1 || $extracase2){
	$temp1 = $seq1len + 1;
	$temp2 = $seq2len + 1;
	$app_str = $app_str."seq1 0 $clipleft1; seq2 0 $clipleft2; score=$sentinelleft (+)\n";
	$app_str = $app_str."seq1 $clipright1 $temp1; seq2 $clipright2 $temp2; score=$sentinelright (+)\n";
    }

# APPEND HITS FROM $app_str TO LOCAL ALIGNMENT LIST

    open (OFILE, ">>$$.anchtemp");
    print OFILE $app_str;
    close (OFILE);

#    `wc $$.anchtemp` =~ /(\d+)/x;
#    $totalanchs = $totalanchs + $1;	
#    print STDERR "CHAOS hits\n";
#    print STDERR `cat $$.anchtemp`;

# FIND MAXIMAL-SCORING CONSISTENT CHAIN

    `$lagandir/anchors $$.anchtemp $gfc $anchparams | sort -n +1 > $$.anchs.sorted`;
    if ($?) { exit(1); }

# IF WE'RE DONE, THEN QUIT!

    $stillmore = ($recurfl =~ /\((\d+)\,(\d+)\,(\d+)\,(\d+)\)(\w*)(.*)/);
    if (!$stillmore && !$extracase1 && !$extracase2) { 
	last; 
    }
    
#    `wc $$.anchs` =~ /(\d+)/x;
#    print STDERR "ANCHS hits\n";
#    print STDERR `cat $$.anchs.sorted`;
#    $goodanchs = $goodanchs + $1;

#    if ($?) { exit(1); }

# READ SORTED ANCHORS TO @anchors

    open(SFILE, "$$.anchs.sorted");
    @anchors = <SFILE>;
    close(SFILE);

    @b1new = 0;
    @b2new = 0;
    @e1new = 0;
    @e2new = 0;
    @scores = 0;

    $app_str = "";
    
    # FOR EACH UNALIGNED REGION

    $area = 0;
    $maxarea = 0;
    $k = 0;
    
    for ($m = 0; $m < @anchors; $m++){

	# SAVE OLD ANCHORS (SKIP FIRST AND LAST FAKE ANCHORS)

	if ($m >= 1 && $m < @anchors - 1){
	    $anchors[$m] =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) (.*)/;
	    $score = $5; chomp $score;
	    $app_str = $app_str."seq1 $1 $2; seq2 $3 $4; score=$score (+)\n";
	}

	if ($m == 0){ next; }

	# DETERMINE REGION BOUNDARIES

	$anchors[$m-1] =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) (.*)/;
	$gap1begin = $2 + 1;
	$gap2begin = $4 + 1;
	$prevanchorscore = $5; chomp $prevanchorscore;

	$anchors[$m] =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) (.*)/;
	$gap1end = $1 - 1;
	$gap2end = $3 - 1;
	$nextanchorscore = $5; chomp $nextanchorscore;

	# CHECK IF RECURSION NEEDED
	
	$boxarea = ($gap1end - $gap1begin + 1) * ($gap2end - $gap2begin + 1);
	$area = $area + $boxarea;
	$maxarea = $boxarea if ($boxarea > $maxarea);

	if ($boxarea >= $minbox && ($gap1end - $gap1begin + 1) > $minside &&
	    ($gap2end - $gap2begin + 1) > $minside ){

	    # FAST REJECT
	    
	    if ($fastreject && ($i >= $frminlevel) && ($i <= $frmaxlevel)){

		# SKIP MARKED ENDS OF ALIGNMENT

		if ($nextanchorscore == $sentinelleft ||
		    $prevanchorscore == $sentinelright){
		    next;
		}

		# TRIM NEW ENDS OF ALIGNMENT
		
		if ($prevanchorscore == $sentinelleft){
#		    if ($boxarea > $frseq1[$i] * $frseq2[$i]){
		    if (($gap1end - $gap1begin > $frseq1[$i]) ||
			($gap2end - $gap2begin > $frseq2[$i])){
			if (@anchors == 2){ exit(3); }
			$clipleft1 = max ($gap1begin-1, $gap1end - $frseq1[$i]);
			$clipleft2 = max ($gap2begin-1, $gap2end - $frseq2[$i]);
			$gap1begin = $clipleft1 + 1;
			$gap2begin = $clipleft2 + 1;
		    }
		}
		elsif ($nextanchorscore == $sentinelright){
#		    if ($boxarea > $frseq1[$i] * $frseq2[$i]){
		    if (($gap1end - $gap1begin > $frseq1[$i]) ||
			($gap2end - $gap2begin > $frseq2[$i])){
			if (@anchors == 2){ exit(3); }
			$clipright1 = min ($gap1end+1, $gap1begin + $frseq1[$i]);
			$clipright2 = min ($gap2end+1, $gap2begin + $frseq2[$i]);
			$gap1end = $clipright1 - 1;
			$gap2end = $clipright2 - 1;
		    }
		}
	    }

	    # ADD REGION

	    if ($gap1begin < $gap1end && $gap2begin < $gap2end){
		$b1new[$k] = $gap1begin;
		$b2new[$k] = $gap2begin;
		$e1new[$k] = $gap1end;
		$e2new[$k] = $gap2end;
		$k++;
	    }
	}
    }

    @b1 = @b1new;
    @b2 = @b2new;
    @e1 = @e1new;
    @e2 = @e2new;
    if ($debug) {
	print STDERR "Level $i Summary:\n";
	print STDERR "   Using $seq1 $seq2 ($wordlen, $degeneracy, $cutoff)\n";
	if ($totalanchs == 0) {
	    $percentage = 0;
	}
	else {
	    $percentage = $goodanchs / $totalanchs * 100.0;
	}
	print STDERR "   $goodanchs good out of $totalanchs total anchors ($percentage%)\n";
	$area = $area / 1000000;
	$maxarea = $maxarea / 1000000;
	print STDERR "   Total area left = $area (max = $maxarea)\n";
    }
    $cumanchs = $cumanchs + $goodanchs;
    $i++;
}

$res = `sort -nr +1 $$.anchs.sorted`;
if ($?) { exit(1); }

`rm $$.*`;

if($tofile) {
    open(OUTFILE, ">$filename");
    print OUTFILE "$res";
    close OUTFILE;
}
else {
    print "$res";
}

print STDERR "$cumanchs cumulative anchors\n"

