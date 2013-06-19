#!/usr/bin/env perl


# This script requires the environment variables:
# LAGAN_DIR and VISTA_DIR

($lagandir = $ENV{LAGAN_DIR}) or die "LAGAN_DIR not set";

$paregmin = 75;
$paregmax = 100;
$pamin = 50;

$pbases = 10000;
$ptickdist = 2000;
$presolution = 25;
$pwindow = 40;
$pnumwindows = 4;


if (@ARGV < 2) {
    print ("usage:\n mviz.pl data_file param_file [plotfile]\n\n");
    exit(1);
}

$pfspec = 0;
if (@ARGV==3) {
    $pfspec = 1;
    $plotfile=@ARGV[2];
    print "Using VISTA plotfile: $plotfile\n";
}


$filename = $ARGV[1];
open(PARAMFILE, "$filename") || die "Could not open $filename.\n\n";

$i=0;
$j=0;
$k=0;
$filespec = 0;
while ($line = <PARAMFILE>) {
    chomp $line;
    if ((substr($line, 0, 1) ne "#") && ($line ne "")) {
	if (!$filespec) {
	    $seqfile = $line;
	    $filespec = 1;
	} elsif (substr($line,0,1) eq "-") {
	    if (substr($line,0,2) eq "--") {
		@vparams[$j++] = $line;
	    } else {
		@params[$i++] = $line;
	    }
	} else {
	    @targets[$k++] = $line;
	}
    }
}

$seqfile = @ARGV[0];

if ($lagandir eq "") {
    print ("Must specify environment variable LAGAN_DIR\n");
    exit(1);
}

for ($i=0; $i<@vparams; $i+=2) {
    if (@vparams[$i] eq "--regmin") { $paregmin = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--regmax") { $paregmax = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--min") { $pamin = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--bases") { $pbases = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--tickdist") { $ptickdist = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--resolution") { $presolution = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--window") { $pwindow = @vparams[$i+1]; }
    elsif (@vparams[$i] eq "--numwindows") { $pnumwindows = @vparams[$i+1]; }
}

open(FASTAFILE, "$seqfile") || die "Could not open $seqfile.\n\n";

$prefix = substr $seqfile, 0, (rindex $seqfile, ".");
if (substr($prefix, -1, 1) ne "_") {$prefix = "$prefix\_";}

$line = <FASTAFILE>;
chomp $line;

while (substr($line, 0, 1) ne ">") {
    $line = <FASTAFILE>;
    chomp $line;
}

$i=0;
%list=();

if (substr($line, 0, 1) eq ">") {
    @keys[$i] = substr($line, 1);

    $list{@keys[$i]}=$i;

    if (@targets == 0) {
	@targets[0] = @keys[$i];
	print "Setting Base Sequence: @targets[0]\n";
    }
} else {
    print ("$filename is NOT a Multi-FASTA file...\n");
    exit(1);
}

while ($line = <FASTAFILE>) {
    chomp $line;

    if (substr($line, 0, 1) eq ">") {
	$i++;
	@keys[$i] = substr($line, 1);

	$list{@keys[$i]}=$i;
    }
}

if ((@targets > 1)) { 
    
    $j=0;
    for ($i=1; $i<@targets; $i++) {
	$_ = @targets[$i];
	@bp[$j++]=/\w+/g;
	$_=$&;
	@bp[$j++]=/\w+/g;
    }
    $j=1;
    foreach $s (@bp) { 
	@targets[$j++]=$s;
    }
    if (@targets %2 != 1) {
	$c = @targets;
	print ("$c sequences: ");
	print ("Must specify single base sequence\n");
	print (" OR base sequence and pairs of sequences.\n");
	exit(1);
    }
}

$i=0;
if (@targets == 1) {
    foreach $s (@keys) {
	$s = substr $s, 0, (rindex $s, "_aligned");
	if ($s ne @targets[0]) {
	    @targets[++$i] = @targets[0];
	    @targets[++$i] = $s;	    
	}
    }
}

print "TARGETS:\n";foreach $s (@targets) { print "\"$s\"\n"; }

$prjhead = "$lagandir/utils/mproject.pl $seqfile";
$binhead = "$lagandir/utils/mf2bin.pl";
$j=0;
for($i=1; $i<@targets; $i+=2) {
    $outprefix = "$prefix@targets[$i]\_@targets[$i+1]";
    $pargs = "$targets[$i]_aligned $targets[$i+1]_aligned";
    $pstr = "$prjhead $pargs > $outprefix.prj";
    print "$pstr\n";
    if(`$pstr`) { print "\nprojection failure...\n"; exit(1); }
    $bstr = "$binhead $outprefix.prj -out $outprefix.bin";
    print "$bstr\n";
    if(`$bstr`) { print "\npacking failure...\n"; exit(1); }
    @bins[$j++] = "$outprefix.bin";
    print "\n";
}

%distinct=();
foreach $s (@targets) {
    $distinct{$s} = 0;
}

@dseqs = keys %distinct;

if (!$pfspec) {

    $plotfile = "$prefix.plotfile";
    open (PLOTFILE, ">$plotfile");

    print PLOTFILE "TITLE $prefix.fa - mlagan\n\n";
    print PLOTFILE "OUTPUT $prefix.pdf\n\n";

    print PLOTFILE "SEQUENCES ";
    foreach $s (@dseqs) {
	print PLOTFILE "$s ";
    }
    print PLOTFILE "\n\n";
    
    $i=1;
    foreach $s (@bins) {
	print PLOTFILE "ALIGN $s BINARY\n";
	print PLOTFILE " SEQUENCES @targets[$i] @targets[$i+1]\n";
	print PLOTFILE " REGIONS $paregmin $paregmax\n";
	print PLOTFILE " MIN $pamin\n";
	print PLOTFILE "END\n\n";   
	$i+=2;
    }
    
    print "touch $prefix.ann\n\n";
    `touch $prefix.ann`;
    
    print PLOTFILE "GENES $prefix.ann\n\n";    
    print PLOTFILE "LEGEND on\n\n";
    print PLOTFILE "COORDINATE @targets[0]\n\n";
    print PLOTFILE "PAPER letter\n\n";
    print PLOTFILE "BASES $pbases\n\n";
    print PLOTFILE "TICK_DIST $ptickdist\n\n";
    print PLOTFILE "RESOLUTION $presolution\n\n";
    print PLOTFILE "WINDOW $pwindow\n\n";
    print PLOTFILE "NUM_WINDOWS $pnumwindows\n\n";

}

($vistadir = $ENV{VISTA_DIR}) or die "VISTA_DIR not set";

$vistastr = "$vistadir/RunVista $plotfile";
print "$vistastr\n";
if (!`$vistastr`) { print "\nVISTA failure...\n"; exit(1); }

print "\n\nmviz.pl -- end.\n\n";


