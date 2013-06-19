#!/usr/bin/env perl
use File::Basename;

$lagandir = $ENV{LAGAN_DIR};
$pid = $$;

# process arguments
if (@ARGV < 4 && @ARGV > 6) {
    print STDERR ("usage:\n cmerge seqfile mfafile draftfile outfile [-nocrop] [-skipfr pid]\n");
    exit(1);
}
$arglist = "";
$nocrop = 0;
for ($i = 4; $i < @ARGV; $i++) {
    if ($ARGV[$i] =~ /-nocrop/){
	$nocrop = 1;
    }
    elsif ($ARGV[$i] =~ /-skipfr/){
	$skipfr = 1;
	$pid = $ARGV[++$i];
	chomp $pid;
    }
    else {
	print STDERR "Bad arg to cmerge: $ARGV[$i]";
	exit(1);
    }
}
$arglist = "$arglist $recurfl";

if (!$skipfr) {
    exit(1);
}
$newdir = `pwd`;
chomp $newdir;
$newdir = "$newdir/$pid";

open (LOGFILE, ">>$newdir/log");
open (INFOFILE, ">$newdir/minfo");

print STDERR ("\n");
print STDERR ("Computing Contig Overlaps\n");
print STDERR ("-------------------------\n");

print LOGFILE ("\n");
print LOGFILE ("Computing Contig Overlaps\n");
print LOGFILE ("-------------------------\n");

# initialize merged file
open (OFILE, ">$ARGV[3]");
print OFILE (">merged\n");
close (OFILE);
`cp $ARGV[3] $ARGV[3].masked`;

# initialize padding file
open (OFILE, ">$newdir/padding");
print OFILE (">padding\n");
print OFILE ("NNNNNNNNNNNNNNNNNNNN.NNNNNNNNNNNNNNNNNNNN\n");
close (OFILE);
$padlength = `$lagandir/utils/getlength $newdir/padding`; chomp $padlength;

# other initialization
$totlength = `$lagandir/utils/getlength $ARGV[0]`;
chomp $totlength;
$mergedEnd = 0;

# read contig list
$numContigs = 0;
@list = `cat $ARGV[2]`;

for ($i = 3; $i < @list; $i++){
    $list[$i] =~ /(.*)\.mfa --\> \((\d+) (\d+)\) score=(\d+), offset=\((\d+) (\d+)\), index=(\d+)/;
    $filenames[$i-3] = $1;
    $seq1Begin[$i-3] = $2;
    $seq1End[$i-3] = $3;
    $score[$i-3] = $4;
    $s1shifts[$i-3] = $5;
    $s2shifts[$i-3] = $6;
    $num[$i-3] = $7;


    $temp = $seq1Begin[$i-3] - $s1shifts[$i-3];
    $seq2Begin[$i-3] = `$lagandir/utils/getcontigpos $filenames[$i-3].mfa $temp`; chomp $seq2Begin[$i-3];
    $seq2Begin[$i-3] += $s2shifts[$i-3];

    $temp = $seq1End[$i-3] - $s1shifts[$i-3];
    $seq2End[$i-3] = `$lagandir/utils/getcontigpos $filenames[$i-3].mfa $temp`; chomp $seq2End[$i-3];
    $seq2End[$i-3] += $s2shifts[$i-3];

    print STDERR "$filenames[$i-3].mfa --> $seq1Begin[$i-3] $seq1End[$i-3] $score[$i-3] $s1shifts[$i-3] $s2shifts[$i-3] $num[$i-3] $seq2Begin[$i-3] $seq2End[$i-3]\n";

    $numContigs++;
}

# extract contigs
$contigfile = basename ($ARGV[1]);
$contigdir = dirname ($ARGV[1]);
$newdir = `pwd`;
chomp $newdir;
$newdir = "$newdir/$pid";

# start out merged file with only padding
`mv $ARGV[3] $ARGV[3].new`;
`$lagandir/utils/seqmerge $ARGV[3].new $newdir/padding > $ARGV[3]`;
`mv $ARGV[3].masked $ARGV[3].masked.new`;
`$lagandir/utils/seqmerge $ARGV[3].masked.new $newdir/padding > $ARGV[3].masked`;
$contigStart[0] = 1;
$startChop[0] = 0;

`cp $filenames[0] $newdir/current`;
`cp $filenames[0].masked $newdir/current.masked`;

# merge contigs
for ($i = 1; $i < $numContigs; $i++){
    `$lagandir/rechaos.pl $newdir/current $filenames[$i] -recurse \"(12,0,40,0)x\" -maskedonly > $newdir/currentanchs`;
    # find the overlap

    `$lagandir/utils/getoverlap $newdir/currentanchs` =~ /(-?\d+) (-?\d+) (-?\d+) (-?\d+)/;
    $rangebegin1 = $1; 
    $rangeend1 = $2;
    $rangebegin2 = $3;
    $rangeend2 = $4;

    chomp $rangebegin1;
    chomp $rangeend1;
    chomp $rangebegin2;
    chomp $rangeend2;

    $thislength = `$lagandir/utils/getlength $filenames[$i-1]`; chomp $thislength;
    $nextlength = `$lagandir/utils/getlength $filenames[$i]`; chomp $nextlength;
    
    # if no overlap, flush the buffer
    if ($rangebegin1 == -1 && $rangeend1 == -1){

	print STDERR "No overlap found...\n";

	`mv $ARGV[3] $ARGV[3].new`;
	`$lagandir/utils/seqmerge $ARGV[3].new $newdir/current $newdir/padding > $ARGV[3]`;
	`cp $filenames[$i] $newdir/current`;

	`mv $ARGV[3].masked $ARGV[3].masked.new`;
	`$lagandir/utils/seqmerge $ARGV[3].masked.new $newdir/current.masked $newdir/padding > $ARGV[3].masked`;
	`cp $filenames[$i].masked $newdir/current.masked`;

	$contigEnd[$i-1] = $contigStart[$i-1] + $thislength - 1;
	$contigStart[$i] = $contigEnd[$i-1] + $padlength + 1;
	$endChop[$i-1] = 0;
	$startChop[$i] = 0;
    }
    else {
	print STDERR "Overlap detected!\n";

	# extract the overlapped region > overlap
	$j = $rangebegin1 - 1;

	if ($j > 0){
	    `$lagandir/utils/cextract $newdir/current 1 $j 0 0 > $newdir/overlap`;
	    `$lagandir/utils/cextract $newdir/current.masked 1 $j 0 0 > $newdir/overlap.masked`;
	    $overlaplength = `$lagandir/utils/getlength $newdir/overlap`; chomp $overlaplength;
	    
	    `mv $ARGV[3] $ARGV[3].new`;	
	    `$lagandir/utils/seqmerge $ARGV[3].new $newdir/overlap > $ARGV[3]`;
	    `mv $ARGV[3].masked $ARGV[3].masked.new`;	
	    `$lagandir/utils/seqmerge $ARGV[3].masked.new $newdir/overlap.masked > $ARGV[3].masked`;
	}
	    
	# extract the nonoverlapped region > current
	`$lagandir/utils/cextract $filenames[$i] $rangebegin2 $nextlength 0 0 > $newdir/current`;
	`$lagandir/utils/cextract $filenames[$i].masked $rangebegin2 $nextlength 0 0 > $newdir/current.masked`;

	$contigEnd[$i-1] = $contigStart[$i-1] + $overlaplength - 1;
	$contigStart[$i] = $contigEnd[$i-1] + 1;
	$endChop[$i-1] = $thislength - $rangeend1;
	$startChop[$i] = $rangebegin2 - 1;
    }

    if (index ($filenames[$i-1], ".rc") == -1) { $direction = "+"; } else { $direction = "-"; }
    @temp = `head $filenames[$i-1]`;
    chomp $temp[0]; $temp[0] = substr $temp[0], 1;

    print INFOFILE "$temp[0]\n";
    print INFOFILE "$num[$i-1] $seq1Begin[$i-1] $seq1End[$i-1] $contigStart[$i-1] $contigEnd[$i-1] $startChop[$i-1] $endChop[$i-1] $direction $score[$i-1] $seq2Begin[$i-1] $seq2End[$i-1]\n";

}

$thislength = `$lagandir/utils/getlength $filenames[$numContigs - 1]`; chomp $thislength;
$contigEnd[$numContigs - 1] = $contigStart[$numContigs - 1] + $thislength - 1;
$endChop[$numContigs - 1] = 0;

`mv $ARGV[3] $ARGV[3].new`;
`$lagandir/utils/seqmerge $ARGV[3].new $newdir/current $newdir/padding > $ARGV[3]`;
`mv $ARGV[3].masked $ARGV[3].masked.new`;
`$lagandir/utils/seqmerge $ARGV[3].masked.new $newdir/current.masked $newdir/padding > $ARGV[3].masked`;

if (index ($filenames[$numContigs - 1], ".rc") == -1) { $direction = "+"; } else { $direction = "-"; }
@temp = `head $filenames[$numContigs - 1]`;
chomp $temp[0]; $temp[0] = substr $temp[0], 1;
print INFOFILE "$temp[0]\n";
print INFOFILE "$num[$numContigs - 1] $seq1Begin[$numContigs - 1] $seq1End[$numContigs - 1] $contigStart[$numContigs - 1] $contigEnd[$numContigs - 1] $startChop[$numContigs - 1] $endChop[$numContigs - 1] $direction $score[$numContigs - 1] $seq2Begin[$numContigs - 1] $seq2End[$numContigs - 1]\n";


print STDERR "Merging complete!\n\n";
print LOGFILE "Merging complete!\n\n";

# 1. write getoverlap() -- given a set of chaos hits, find the beginning and end in both seqs
# 2. implement contigStart, contigStop -- positions of the contig begins/ends in the merged draft sequence
# 3. startChop, endChop -- number chopped from each end
# 4. secFrom, secTo -- pos in the chopped contig sequence
