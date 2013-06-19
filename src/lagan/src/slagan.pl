#!/usr/bin/perl -w

use strict;

my $lagandir = $ENV{LAGAN_DIR};

if (@ARGV < 2) {
	print ("Usage:\n slagan.pl seqfile1 seqfile2 [-glocal \"glocal flags\"] [-chaos \"chaos flags\"] [-order \"order flags\"] [-recurse \"(wl1,nd1,co1),(wl2,nd2,co2),...\"] [-mfa] [-out \"filename\"] [-maskedonly] [-debug] [-translate] [-fastreject]\n");
	exit(1);
}

my ($seq1, $firstName) = ($ARGV[0], $ARGV[0]);
die("$0: File not found: $seq1. Stopped") unless -f $seq1;
my ($seq2, $secondName) = ($ARGV[1], $ARGV[1]);
die("$0: File not found: $seq2. Stopped") unless -f $seq2;

my ($extra1, $extra2) =(0, 0);
if (-e "$seq1.masked") { $seq1 = "$seq1.masked"; $extra1 = 1;}
if (-e "$seq2.masked") { $seq2 = "$seq2.masked"; $extra2 = 1;}

my $max_ext = 25000;
my $ext_mul = 1;
my $arglist = "";
my $glocal_fl = " -gapopen 0,1000,2000,2000 -gapcont 0.2,0.06,0.06,0.06 -dist 0,1.0,2.5,2.5";
my $chaos_fl = " -wl 11 -nd 1 -co 10 -ext -rsc 2250 -b";
my $lagan_fl = "";
my $supermap_fl = "-glocal_out=slagan.out.glocal";
my $outfile = 0;
my $fastrej = 0;

for (my $i = 2; $i < @ARGV; $i++) {
	if ($ARGV[$i] =~ /-glocal_fl/) {
		$glocal_fl = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-chaos_fl/) {
		$chaos_fl = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-lagan_fl/) {
		$lagan_fl = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-max_ext/) {
		$max_ext = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-ext_mul/) {
		$ext_mul = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-out/) {
		$outfile = $ARGV[++$i];
		if (-e "$outfile") { system("rm $outfile") and exit(1); }
	} elsif ($ARGV[$i] =~ /-order/) {
		$arglist = $arglist." -order $ARGV[++$i]";
	} elsif (($ARGV[$i] =~ /-gs/) || ($ARGV[$i] =~ /-gc/) || ($ARGV[$i] =~ /-mt/) || ($ARGV[$i] =~ /-ms/) || ($ARGV[$i] =~ /-bw/)) {
		$arglist = $arglist." ".$ARGV[$i];
		$arglist = $arglist." ".$ARGV[++$i];
	} elsif ($ARGV[$i] =~ /-ext/) {
		$arglist = $arglist." -ext $ARGV[++$i]";
	} elsif ($ARGV[$i] =~ /-maskedonly/) {
		$arglist = $arglist." -maskedonly";
	} elsif ($ARGV[$i] =~ /-translate/) {
		$arglist = $arglist." -translate";
	} elsif ($ARGV[$i] =~ /-fastreject/) {
		$fastrej = 1;
#		$arglist = $arglist." -fastreject";
	} elsif ($ARGV[$i] =~ /-recurse/) {
		$arglist = $arglist." -recurse \"".$ARGV[++$i]."\"";
	} elsif ($ARGV[$i] =~ /-chaos/) {
		$chaos_fl = $ARGV[++$i];
	} else {
		die("$0: Invalid option for rlagan: $ARGV[$i]");
	}
}

my $seq1len = `$lagandir/utils/getlength $firstName`;
my $seq2len = `$lagandir/utils/getlength $secondName`;
chomp $seq1len;
chomp $seq2len;

`$lagandir/chaos $seq1 $seq2 $chaos_fl > chaos.$$`;
if ($?) { exit(1); }

#`$lagandir/glocal chaos.$$ $glocal_fl > out.$$`;
#@regs = `$lagandir/anal_gloc.pl < out.$$`;
#print @regs;

open(FH, "> seq1len"); print FH $firstName." ".$seq1len."\n"; close FH;
open(FH, "> seq2len"); print FH $secondName." ".$seq2len."\n"; close FH;
my $supermap_outfile = "slagan.out.smap";
my $supermap_inv = "$lagandir/supermap.pl -sizes1=seq1len -sizes2=seq2len $supermap_fl chaos.$$ -no_clust_run -f -out=$supermap_outfile 1>&2";
#print $supermap_inv."\n";
system($supermap_inv);

open(FH, "< $supermap_outfile");
my @regs = <FH>;
die("$0: Supermap generated no regions. Stopped") unless scalar @regs;
close FH;
unlink "seq1len"; unlink "seq2len"; # unlink $supermap_outfile;

#$prevend1 = $seq1len;
#$prevend2 = $seq2len;
#$nextstart1 = 1;
#$nextstart2 = 1;

for (my $k = 0; $k < @regs; $k++) {
	$regs[$k] =~ /^([^\s]+)\s([\d]+)\s([\d]+)\s\s\s([^\s]+)\s([\d]+)\s([\d]+)\s(\+|\-)\s\((DM|M1|M2),\s([\d]+)\saligns\)$/o;

	my ($startreg1, $endreg1, $startreg2, $endreg2, $strand, $type) = ($2, $3, $5, $6, $7, $8);

=head1
	$regs[$k] =~ /.* Region \[(\d+) (\d+)\]\[(\d+) (\d+)\] (.*) (.)/;
	$startreg1 = $1; $endreg1 = $2; $startreg2 = $3; $endreg2 = $4;
	$strand = $6;
	if ($k+2 < @regs) {
		$regs[$k+1] =~ /.* Region \[(\d+) (\d+)\]\[(\d+) (\d+)\] (.*) (.)/;
		$nextstart1 = $2;
	} else {
		$nextstart1 = 1;
	}
	$y1 = $prevend1-$endreg1;
	$y2 = $startreg1-$nextstart1;
	$expandback = ($max_ext < $y1)? $max_ext:$prevend1-$endreg1;
	$expandforw = ($max_ext < $y2)? $max_ext:$startreg1-$nextstart1;
	$prevend1 = $startreg1;
	$startreg1 = $startreg1 - $expandforw;
	$endreg1 = $endreg1 + $expandback;
=cut

	my $rcf = "";
	if ($strand eq "+") {
#		$endreg2 = ($endreg2 + $expandback * $ext_mul > $prevend2)? $prevend2:($endreg2 + $expandback * $ext_mul);
#		$startreg2 = ($startreg2 - $expandforw * $ext_mul < $nextstart2)? $nextstart2:($startreg2 - $expandforw * $ext_mul);
	} else {
		$rcf = "-rc";
#		$endreg2 = ($endreg2 + $expandforw * $ext_mul > $prevend2)? $prevend2:($endreg2 + $expandforw * $ext_mul);
#		$startreg2 = ($startreg2 - $expandback * $ext_mul < $nextstart2)? $nextstart2:($startreg2 - $expandback * $ext_mul);
	}

#print "$lagandir/utils/fa2xfa $firstName $startreg1 $endreg1 1 > seq1$k.$$\n";
	`$lagandir/utils/fa2xfa $firstName $startreg1 $endreg1 1 > seq1$k.$$\n`;
#print "$lagandir/utils/fa2xfa $secondName $startreg2 $endreg2 2 $rcf > seq2$k.$$\n";
	`$lagandir/utils/fa2xfa $secondName $startreg2 $endreg2 2 $rcf > seq2$k.$$\n`;
#	if ($extra1) { `$lagandir/utils/fa2xfa $seq1 $startreg1 $endreg1 1 > seq1$k.$$.masked\n`; }
#	if ($extra2) { `$lagandir/utils/fa2xfa $seq2 $startreg2 $endreg2 2 $rcf > seq2$k.$$.masked\n`; }
#print "$lagandir/lagan.pl seq1$k.$$ seq2$k.$$ $arglist $lagan_fl -mfa -out lagan.$k.$$\n";
	`$lagandir/lagan.pl seq1$k.$$ seq2$k.$$ $arglist $lagan_fl -mfa -out lagan.$k.$$\n`;

	my $suff = "";
	if ($outfile) { $suff = " >> $outfile"; }
	if (-e "lagan.$k.$$") {
		if ($fastrej) {
#print "$lagandir/utils/scorealign lagan.$k.$$ 45 -cropxmfa -ibounds $suff\n";
			print `$lagandir/utils/scorealign lagan.$k.$$ 45 -cropxmfa -ibounds $suff`;
		} else {
#print "$lagandir/utils/scorealign lagan.$k.$$ 45 -ibounds\n";
			my $sc = `$lagandir/utils/scorealign lagan.$k.$$ 45 -ibounds`;
			chomp($sc);
			if ($sc) {
				print `cat lagan.$k.$$ $suff`;
				print `echo \"=$sc $type\n\" $suff`;
			}
		}
	}
}

my ($outName1, $outName2) = ($ARGV[0], $ARGV[1]);
$outName1 =~ s/^.*\///;
$outName1 =~ s/\..*//;
$outName2 =~ s/^.*\///;
$outName2 =~ s/\..*//;

`cat chaos.$$ > ${outName1}_$outName2.chaos`;
####`cat out.$$ > ${outName1}_$outName2.mon`;
unlink(glob("*.$$"));
if ($extra1 || $extra2) { `rm *.$$.masked`; }
exit(0);


# out: .chaos .mon->.smap .xmfa
