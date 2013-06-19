#!/usr/bin/env perl

$lagandir = $ENV{LAGAN_DIR};
$consrate = 45;
$consupperrate = 65;

if (@ARGV < 2) {
    print ("usage:\n lagan seqfile1 seqfile2 [-chaos \"chaos flags\"] [-order \"order flags\"] [-recurse \"(wl1,nd1,co1,rsc1),(wl2,nd2,co2,rsc2),...\"] [-bin] [-mfa] [-out \"filename\"] [-lazy] [-maskedonly] [-debug] [-usebounds] [-rc] [-translate] [-draft] [-info] [-fastreject]\n");
    exit(1);
}

$firstName = $ARGV[0];
$secondName = $ARGV[1];
$rcFlag = 0;
$arglist = "";
$contigflag = 0;
$infofile = 0;
$okformat = 0;
$binfile = 0;
$infofilename = "alignment";
$direction = "+";
$gfc = " -gfc ";
$rundraft = 0;
$draftparams = "";
$dofastreject = 0;
$doxmfa = 0;
$filename = "";
$format = "";

for ($i = 2; $i < @ARGV; $i++) {
    if ($ARGV[$i] =~ /-order/) {
	$orderfl = $ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-bin/) {
	$orderfl = $orderfl." -bin";
	$binfile = 1;
	$okformat = 1;
    }
    elsif ($ARGV[$i] =~ /-info/) {
	$infofile++;
    }
    elsif ($ARGV[$i] =~ /-mfa/) {
	$orderfl = $orderfl." -mfa";
	$okformat = 1;
    }
    elsif ($ARGV[$i] =~ /-xmfa/) {
	$orderfl = $orderfl." -xmfa";
	$doxmfa = 1;
	$okformat = 1;
    }
    elsif ($ARGV[$i] =~ /-out/) {
	$filename = $ARGV[++$i];
	$infofile++;
	$infofilename = $ARGV[$i];
    }
    elsif (($ARGV[$i] =~ /-gs/) || ($ARGV[$i] =~ /-gc/) || ($ARGV[$i] =~ /-mt/) || ($ARGV[$i] =~ /-ms/) || ($ARGV[$i] =~ /-bw/)){
	$orderfl = $orderfl." ".$ARGV[$i];
	$orderfl = $orderfl." ".$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-s1/) {
	$orderfl = $orderfl." -s1 $ARGV[++$i]";
	$orderfl = $orderfl." ".$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-maskedonly/) {
	$arglist = $arglist." -maskedonly";
    }
    elsif ($ARGV[$i] =~ /-translate/) {
	$arglist = $arglist." -translate";
	$draftparams = $draftparams." -translate";
    }
    elsif ($ARGV[$i] =~ /-fastreject/) {
    	$arglist = $arglist." -fastreject";
	$dofastreject = 1;
	$doxmfa = 1;
	$okformat = 1;
    }
    elsif ($ARGV[$i] =~ /-draftreject/) {
    	$draftparams = $draftparams." -fastreject";
    }
    elsif ($ARGV[$i] =~ /-gap/) {
	$arglist = $arglist." -gap ".$ARGV[++$i];
	$arglist = $arglist." ".$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-recurse/) {
	$arglist = $arglist." -recurse \"".$ARGV[++$i]."\"";
    }
    elsif ($ARGV[$i] =~ /-chaos/) {
	$arglist = $arglist." -chaos \"".$ARGV[++$i]."\"";
    }
    elsif ($ARGV[$i] =~ /-usebounds/) {
	$contigflag = 1;
    }
    elsif ($ARGV[$i] =~ /-rc/) {
	`$lagandir/utils/rc < $ARGV[1] > $ARGV[1].rc`;
	if ($?) { exit(1); }
	$secondName = "$ARGV[1].rc";
	if (-e "$ARGV[1].masked") { 
	    `$lagandir/utils/rc < $ARGV[1].masked > $ARGV[1].rc.masked`;
	    if ($?) { exit(1);} 
	}
	$rcFlag = 1;
	$direction = "-";
    }
    elsif ($ARGV[$i] =~ /-draft/){
	$rundraft = 1;
    }
    elsif ($ARGV[$i] =~ /-cons/){
	$draftparams = $draftparams." -cons $ARGV[$++i]";
    }
    elsif ($ARGV[$i] =~ /-draftskipfr/){
	$draftparams = $draftparams." -skipfr $ARGV[$++i]";
    }
    elsif ($ARGV[$i] =~ /-lazy/){
	$draftparams = $draftparams." -cons $ARGV[$++i]";
    }

    else {
	print "Invalid option for lagan: $ARGV[$i]";
	exit(1);
    }
}

$arglist = $arglist." -ext ";

if ($rundraft){
    `$lagandir/draft.pl $firstName $secondName $draftparams`;
    if ($?) { exit(1);} 
    $secondName = "merged_seq.fa";
}

# print STDERR "perl $lagandir/rechaos.pl $firstName $secondName $gfc $arglist > $$.anchs.final\n";
`perl $lagandir/rechaos.pl $firstName $secondName $gfc $arglist > $$.anchs.final`;

$ex_val = $? >> 8;
if ($ex_val == 3) { exit(0); }

if ($ex_val) { exit(1); }
if ($contigflag){
    @bounds = `$lagandir/utils/getbounds $$.anchs.final $firstName $secondName`;
    if ($?) { exit(1); }
    chomp $bounds[0];
    print STDERR ("Aligning with bounds: $bounds[0]\n");
    print `$lagandir/order $firstName $secondName $bounds[0] $orderfl -anc $$.anchs.final`;
    if ($?) { exit(1); }
}
else {
    if ($dofastreject){
	if (!$filename) {
	    print STDERR "-fastreject requires -out filename!\n";
	    exit(1);
	}
	open(SFILE, "$$.anchs.final");
	@anchors = <SFILE>;
	close(SFILE);

	$anchors[0] =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) (.*)/;
	$end1 = $1 - 1;
	$end2 = $3 - 1;
	$anchors[@anchors - 1] =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) (.*)/;
	$start1 = $2 + 1;
	$start2 = $4 + 1;
	$bounds = "-s1 $start1 $end1 -s2 $start2 $end2 ";

	@anchors = 0;
	$orderfl = $bounds.$orderfl." -xmfa";
    }
    if (!$okformat) {
	$format = "-bin";
    }

    `$lagandir/order $firstName $secondName $format -out $$.align $orderfl -anc $$.anchs.final`;
    if ($?) { exit(1); }

    if (!$okformat) {
	if ($filename) {
	    `$lagandir/utils/bin2bl $$.align > $filename`;
	}
	else {
	    print `$lagandir/utils/bin2bl $$.align`;
	}
    }
    else {
	if ($filename) {
	    `cat $$.align > $filename`;
	}
	else {
	    print `cat $$.align`;
	}
    }
    if ($dofastreject){
	`$lagandir/utils/scorealign $filename $consrate -ibounds -cropxmfa > $$.temp`;
	if ($?) { exit(1); }
	`mv $$.temp $filename`;
    }
}

$infofile += $okformat;
if ($infofile == 3){
    open (INFOFILE, ">$infofilename.info");
    if ($binfile){
	`$lagandir/utils/bin2mf $infofilename > $infofilename.mfa`;
	if ($?) { exit(1); }
	$infofilename = $infofilename.".mfa";
    }
    @temp = `head $secondName`;
    if ($?) { exit(1); }
    chomp $temp[0]; $temp[0] = substr $temp[0], 1;
    print INFOFILE "$temp[0]\n";

    $len = `$lagandir/utils/getlength $secondName`; chomp $len;
    if ($?) { exit(2); }
    $first = $last = $first2 = $last2 = -1;

    $score = `$lagandir/utils/scorealign $infofilename $consupperrate`; chomp $score;
    if ($?) { exit(3); }
    if ($score > 0){
	$score = `$lagandir/utils/scorealign $infofilename $consrate`; chomp $score;
	if ($?) { exit(4); }
	@temp = `$lagandir/utils/scorealign $infofilename $consrate -bounds 0`; 
	if ($?) { exit(5); }
	$temp[0] =~ /(.*) (.*)/;
	$first = $1; $last = $2;

	@temp = `$lagandir/utils/scorealign $infofilename $consrate -bounds 1`; 
	if ($?) { exit(6); }
	$temp[0] =~ /(.*) (.*)/;
	$first2 = $1; $last2 = $2;
    }

    print INFOFILE "1 $first $last 1 $len 0 0 $direction $score $first2 $last2\n";

    close (INFOFILE);
#    `$lagandir/utils/rm $infofilename` if ($binfile);
}

`rm $secondName` if ($rcflag);
`rm $$.*`;
if ($?) { exit(1); }

exit(0);


