#!/usr/bin/env perl

# This script requires the environment variables:
# LAGAN_DIR and VISTA_DIR

if (@ARGV < 1) {
    print ("usage:\n mrunfile.pl filename [-pairwise] [-vista]\n\n");
    exit(1);
}

($lagandir = $ENV{LAGAN_DIR}) or die "LAGAN_DIR not set";


$filename = $ARGV[0];
open(PARAMFILE, "$filename") || die "Could not open $filename.\n\n";

$pairwise = 0;
$dovista = 0;

for ($l=1; $l<@ARGV; $l++) {
    if ($ARGV[$l] eq "-pairwise") {
	$pairwise = 1;
    }
    elsif ($ARGV[$l] eq "-vista") {
	$dovista = 1;
    }
}

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
	    @seqs[$k++] = $line;
	}
    }
}

if ($lagandir eq "") {
    print ("Must specify environment variable LAGAN_DIR\n");
    exit(1);
}

if ($pairwise) {
    $mexecs = "mrunpairs.pl";
} else {
    $mexecs = "mrun.pl";
}

$mstr = "$lagandir/utils/$mexecs $seqfile";

foreach $s (@params) {
    $mstr = "$mstr $s"
}

foreach $s (@seqs) {
    $mstr = "$mstr $s"
}

foreach $s (@vparams) {
    $mstr = "$mstr $s"
}

print "$mstr\n";
`$mstr`;

if($dovista) {

    $prefix = substr $seqfile, 0, (rindex $filename, ".");
    $prefix = "$prefix\_";
    
    if ($pairwise) {
	$prefix="$prefix\pairwise\_";
    }
    
    $plotfile = "$prefix.plotfile";

    ($vistadir = $ENV{VISTA_DIR}) or die "VISTA_DIR not set";

    $vistastr = "$vistadir/RunVista $plotfile";
    print "$vistastr\n";
    if (!`$vistastr`) { print "\nVISTA failure...\n"; exit(1); }

}

print "\nmrunfile.pl -- end.\n\n";













