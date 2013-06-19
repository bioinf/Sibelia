#!/usr/bin/env perl
use File::Basename;

$lazyflag = 0;
$lagandir = $ENV{LAGAN_DIR};
$recurfl = "-recurse \"(12,0,30,0)x,(13,1,30,0)x,(3,0,30,0)xt,(8,1,30,0)x,(7,1,30,0)x,(7,1,15,0)x\"";
$laganparams = "-maskedonly ";
$anchgapstart = -5;
$anchgapcont = -0.2;
$usebounds = 1;

$startingrate = 65;
$rateinc = 1;
$frlevel = "";
$pid = "mergedir";

if (@ARGV < 2) {
    if ((@ARGV == 1) && ($ARGV[0] =~ /-version/)){
	print STDERR "DRAFT version 0.1\n";
	exit (0);
    }
    else {
	print STDERR ("Usage:\n\ndraft.pl SEQFILE MFAFILE [-cons RATE] [-translate] [-version]\n");
	exit (1);
    }
}

$arglist = "";
$skipfr = 0;
for ($i = 2; $i < @ARGV; $i++) {
    if ($ARGV[$i] =~ /-recurse/){
	$recurfl = " -recurse \"".$ARGV[++$i]."\"";
    }
    elsif ($ARGV[$i] =~ /-skipfr/){
	$skipfr = 1;
	$pid = $ARGV[++$i];
	chomp $pid;
    }
    elsif ($ARGV[$i] =~ /-translate/){
	$recurfl = $recurfl." -translate";
    }
    elsif ($ARGV[$i] =~ /-cons/){
	$startingrate = $ARGV[++$i];
	chomp $startingrate;
    }
    elsif ($ARGV[$i] =~ /-lazy/){
	$lazyflag = 1;
    }
    elsif ($ARGV[$i] =~ /-fastreject/){
	$frarg = " -fastreject $frlevel";
    }
    else {
	print STDERR "Bad arg to draft: $ARGV[$i]";
    }
}

$arglist = "$arglist $recurfl -usebounds $laganparams $frarg";

# create new directory
$newdir = `pwd`;
chomp $newdir;
$newdir = "$newdir/$pid";
`mkdir $newdir` if (!(-e $newdir));

open (LOGFILE, ">$newdir/log");

print STDERR ("\n");
print STDERR ("Finding Contig Alignments\n");
print STDERR ("-------------------------\n");

print LOGFILE ("\n");
print LOGFILE ("Finding Contig Alignments\n");
print LOGFILE ("-------------------------\n");

# extract contigs;
$contigfile = basename ($ARGV[1]);
$contigdir = dirname ($ARGV[1]);

`cp $ARGV[1] $newdir`;
@contigs = `perl $lagandir/mextract.pl $newdir/$contigfile`;
if ($?) { exit(1);} 
for ($i = 0; $i < @contigs; $i++){
    chomp $contigs[$i];
    `$lagandir/utils/rc < $contigs[$i] > $contigs[$i].rc`;
    if ($?) { exit(1); }
}

# extract masked contigs
$maskedname = $ARGV[1].".masked";

if (-e $maskedname){
    $maskedcontigfile = basename ($maskedname);
    `cp $maskedname $newdir`;
    @maskedcontigs = `perl $lagandir/mextract.pl $newdir/$maskedcontigfile -masked`;
    if ($?) { exit(1);} 
    for ($i = 0; $i < @maskedcontigs; $i++){
	chomp $maskedcontigs[$i];
	`$lagandir/utils/rc < $maskedcontigs[$i] > $contigs[$i].rc.masked`;
	if ($?) { exit(1); }
    }
}

# create file storing name of contig stats
open (LFILE, ">$newdir/filenames") if (!$lazyflag);
$num = 0;

for ($i = 0; $i < @contigs; $i++){
    chomp $contigs[$i];
    $skip1 = $skip2 = 0;
    # make alignments
    if (!$lazyflag || !(-e "$contigs[$i].mfa")){
	$execute = "perl $lagandir/lagan.pl $ARGV[0] $contigs[$i] -mfa $arglist -out $contigs[$i].mfa";
	$execute = $execute." -gap $anchgapstart $anchgapcont" if ($usebounds);
	`$execute`;
	$ex_val = $? >> 8;
	if (!(-e "$contigs[$i].mfa")) { $skip1 = 1; }
	elsif ($?) { exit(1);} 

	if (!$skip1 && $usebounds){
	    # compute bounds
	    @bounds = `$lagandir/utils/getbounds anchs.final $ARGV[0] $contigs[$i]`;
	    if ($?) { exit(1);} 
	    $bounds[0] =~ /-s1 (\d+) (\d+) -s2 (\d+) (\d+)/;
	    $s1shift = $1 - 1;
	    $s2shift = $3 - 1;
	}
	`rm anchs.final`;
    }

    if (!$lazyflag || !(-e "$contigs[$i].rc.mfa")){
	$execute = "perl $lagandir/lagan.pl $ARGV[0] $contigs[$i].rc -mfa $arglist -out $contigs[$i].rc.mfa";
	$execute = $execute." -gap $anchgapstart $anchgapcont" if ($usebounds);
	`$execute`;
	$ex_val = $? >> 8;
	if (!(-e "$contigs[$i].rc.mfa")) { $skip2 = 1; }
	elsif ($?) { exit(1);} 
 	if (!$skip2 && $usebounds){
	    # compute bounds
	    @bounds = `$lagandir/utils/getbounds anchs.final $ARGV[0] $contigs[$i].rc`;
	    if ($?) { exit(1);} 
	    $bounds[0] =~ /-s1 (\d+) (\d+) -s2 (\d+) (\d+)/;
	    $s1rcshift = $1 - 1;
	    $s2rcshift = $3 - 1;
	}
	`rm anchs.final`;
    }

    if ($skip1) {
	$fscore = 0;
    }
    else {
	$fscore = `$lagandir/utils/scorealign $contigs[$i].mfa $startingrate`; chomp $fscore;
	if ($?) { exit(1);} 
    }
    if ($skip2) {
	$bscore = 0;
    }
    else {
	$bscore = `$lagandir/utils/scorealign $contigs[$i].rc.mfa $startingrate`; chomp $bscore;
	if ($?) { exit(1);} 
    }
    # pick strand

#    print LFILE "$s1shift $contigs[$i].mfa\n" if (!$lazyflag);
#    print LFILE "$s1rcshift $contigs[$i].rc.mfa\n" if (!$lazyflag);
    
#    if (0){
    if ($fscore > 0 || $bscore > 0){
	$j = $i + 1;
	if ($fscore > $bscore){
	    print STDERR ("(+) direction preferred for Contig \"$contigs[$i]\": $fscore > $bscore\n");
	    print LOGFILE ("(+) direction preferred for Contig \"$contigs[$i]\": $fscore > $bscore\n");
	    print LFILE "$j $s1shift $s2shift $contigs[$i].mfa\n" if (!$lazyflag);
	    print STDERR "$j $s1shift $s2shift $contigs[$i].mfa\n" if (!$lazyflag);
	}
	elsif ($bscore > $fscore){
	    print STDERR ("(-) direction preferred for Contig \"$contigs[$i]\": $fscore < $bscore\n");
	    print LOGFILE ("(-) direction preferred for Contig \"$contigs[$i]\": $fscore < $bscore\n");
	    print LFILE "$j $s1rcshift $s2rcshift $contigs[$i].rc.mfa\n" if (!$lazyflag);
	    print STDERR "$j $s1rcshift $s2rcshift $contigs[$i].rc.mfa\n" if (!$lazyflag);
	}
    }
#    }
    else {
	print STDERR ("Contig \"$contigs[$i]\" could not be matched: $fscore, $bscore\n");
	print LOGFILE ("Contig \"$contigs[$i]\" could not be matched: $fscore, $bscore\n");
    }
}
close (LFILE);

print STDERR ("\n");
print STDERR ("Computing Contig Ordering\n");
print STDERR ("-------------------------\n\n");

print LOGFILE ("\n");
print LOGFILE ("Computing Contig Ordering\n");
print LOGFILE ("-------------------------\n\n");

$foundorder = 0;

for ($cutoff = $startingrate; !$foundorder && ($cutoff < 100); $cutoff += $rateinc){
    `$lagandir/utils/scorecontigs /$newdir/filenames $ARGV[0] $newdir/contignames $cutoff > $newdir/ranges`;
    if ($?) { exit(1);} 
    @list = `cat $newdir/ranges`;
    $list[0] =~ /numContigs = (\d+)/;
    next if ($1 == 0);

    `$lagandir/utils/contigorder $newdir/ranges > $newdir/corder`;
    if ($?) { exit(1);} 
    @list = `cat $newdir/corder`;
    chomp $list[0];
    $foundorder = 1 if ($list[0] ne "ordering failed");
}

if ($foundorder){
    open (OFILE, ">$newdir/draft");
    print OFILE ("Draft Ordering\n");
    print OFILE ("--------------\n\n");
    
    @contignames = `cat $newdir/contignames`;
    for ($i = 0; $i < @contignames; $i++){
	$contignames[$i] =~ /(\d+) (\d+) (\d+) (.*)/;
	$num[$i] = $1; chomp $num[$i];
	$s1shifts[$i] = $2; chomp $s1shifts[$i];
	$s2shifts[$i] = $3; chomp $s2shifts[$i];
	$filenames[$i] = $4; chomp $filenames[$i];
    }

    @list = `cat $newdir/corder`;
    for ($i = 0; $i < @list; $i++){
	$list[$i] =~ /(\d+) --\> \((\d+) (\d+)\) (.*)/;
	$score = $4; chomp $score;
	print OFILE ("$filenames[$1] --> ($2 $3) score=$score, offset=($s1shifts[$1] $s2shifts[$1]), index=$num[$1]\n");
    }
    close (OFILE);
    
    print STDERR `cat $newdir/draft`;
    print LOGFILE `cat $newdir/draft`;
    close (LOGFILE);
}
else {
    print STDERR "Could not compute ordering.";
    print LOGFILE "Could not compute ordering.";
    close (LOGFILE);
    exit (0);
}

$filename1 = $ARGV[0];
$filename2 = "$newdir/$contigfile";

`$lagandir/cmerge2.pl $filename1 $filename2 $newdir/draft $filename2.merged -skipfr $pid`;
if ($?) { exit(1); }

print STDERR "EXECUTE $lagandir/cmerge2.pl $filename1 $filename2 $newdir/draft $filename2.merged -skipfr $pid\n";

`cp $filename2.merged merged_seq.fa`;
`cp $filename2.merged.masked merged_seq.fa.masked`;
`cp $newdir/minfo minfo`;
`cp $newdir/ranges ranges`;
`cp $newdir/log log`;

print STDERR ("\n");
print STDERR ("Computing Final Alignment\n");
print STDERR ("-------------------------\n\n");

# `rm -rf $newdir`;

