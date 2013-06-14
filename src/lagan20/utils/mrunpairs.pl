#!/usr/bin/env perl

# This script requires the environment variables:
# LAGAN_DIR and VISTA_DIR

# VISTA .plotfile defaults

($lagandir = $ENV{LAGAN_DIR}) or die "LAGAN_DIR not set";

$paregmin = 75;
$paregmax = 100;
$pamin = 50;

$pbases = 10000;
$ptickdist = 2000;
$presolution = 25;
$pwindow = 40;
$pnumwindows = 4;


if (@ARGV < 1) {
    print ("usage:\n mrunpairs.pl filename\n");
    print ("options: [base sequence name [sequence pairs]]\n");
    print ("default: [base sequence name = first sequence]\n");
    print ("other MLAGAN parameters...\n");
    print ("other VISTA parameters...\n");
    exit(1);
}

$filename = $ARGV[0];

$i = 1;
$j = 0;
$k = 0;
$l = 0;
$treespec = 0;
while ($i < @ARGV) {
    if ($ARGV[$i] eq "-tree") {
	$treepos = $j+1;
	@params[$j] = "-tree";
	@params[++$j] = "\"$ARGV[++$i]\"";
	$_ = @params[$j];
	$topen = tr/"\("/"\("/;
	$tclose = tr/"\)"/"\)"/;
	$treespec = ($topen == $tclose);
    } else {
	if (substr($ARGV[$i],0,1) eq "-") {
	    if (substr($ARGV[$i],0,2) eq "--") {
		@vparams[$l++] = $ARGV[$i++];
		@vparams[$l++] = $ARGV[$i];
	    } else {
		$j++;
		@params[$j] = $ARGV[$i];
		if ((@params[$j] eq "-gapstart") || 
		    (@params[$j] eq "-gapend") ||
		    (@params[$j] eq "-gapcont") ||
		    (@params[$j] eq "-gapperseq") ||
		    (@params[$j] eq "-match") ||
		    (@params[$j] eq "-mismatch") ||
		    (@params[$j] eq "-overlap") ||
		    (@params[$j] eq "-glwidth")) {
		    @params[++$j] = $ARGV[++$i];
		}
	    }
	} else {
	    @targets[$k++] = $ARGV[$i];
	}
    }
    $i++;
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

if (!$treespec) { 
    $j++;
    $treepos = $j+1;
    @params[$j] = "-tree";
    @params[++$j] = "\"()\"";
}

if ($lagandir eq "") {
    print ("Must specify environment variable LAGAN_DIR\n");
    exit(1);
}

$mextstr = "$lagandir/mextract.pl $filename";
print "$mextstr\n";
if(!`$mextstr`) { print "\nMulti-FASTA extraction failure...\n"; exit(1); }

if (-e "$filename.masked") {
    $mextstr = "$lagandir/mextract.pl $filename.masked -masked";
    print "$mextstr\n";
    if(!`$mextstr`) {
	print "\nMasked Multi-FASTA extraction failure...\n"; 
	exit(1);
    }
}

open(FASTAFILE, "$filename") || die "Could not open $filename.\n\n";

$line = <FASTAFILE>;
chomp $line;

while (substr($line, 0, 1) ne ">") {
    $line = <FASTAFILE>;
    chomp $line;
}

$i=0;
%list=();
$i=0;
%list=();

if (substr($line, 0, 1) eq ">") {
    $_ = substr($line, 1);
    /\w+/g;
    @keys[$i] = $&;
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
	$_ = substr($line, 1);
	/\w+/g;
	@keys[$i] = $&;
	$list{@keys[$i]}=$i;
    }
}

$fprefix = substr $filename, 0, (rindex $filename, ".");
$prefix = "$fprefix\_";
$pprefix = "$fprefix\_pairwise\_";

foreach $s (@keys) {
    @fnames[$list{$s}] = "$prefix$keys[$list{$s}].fa";
}


if ((@targets > 1)) { 
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
	if ($s ne @targets[0]) {
	    @targets[++$i] = @targets[0];
	    @targets[++$i] = $s;	    
	}
    }

}

$j=0;
for($i=1; $i<@targets; $i+=2) {

    $outprefix = "$pprefix@targets[$i]\_@targets[$i+1]";

    $mfiles = " @fnames[$list{@targets[$i]}] @fnames[$list{@targets[$i+1]}]";

    @params[$treepos]="\"(@targets[$i] @targets[$i+1])\"";

    $mparams = "";
    foreach $s (@params) {
	$mparams = "$mparams $s";
    }
    
    $mlagan = "$lagandir/mlagan$mfiles$mparams > $outprefix.out";
    print "\n$mlagan\n\n";
    if(`$mlagan`) { print "\n\n"; exit(1); }

    $binhead = "$lagandir/mpack.pl";
    $bstr = "$binhead $outprefix.out -out $outprefix.bin";
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


$plotfile = "$pprefix.plotfile";
open (PLOTFILE, ">$plotfile");

print PLOTFILE "TITLE $prefix.fa - mlagan\n\n";
print PLOTFILE "OUTPUT $pprefix.pdf\n\n";

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

#$vistadir = `echo \$VISTA_DIR`;
#chomp $vistadir;

#if ($vistadir eq "") {
#    print ("Must specify environment variable VISTA_DIR\n");
#    exit(1);
#}
    
#$vistastr = "$vistadir/RunVista $plotfile";
#print "$vistastr\n";
#if (!`$vistastr`) { print "\nVISTA failure...\n"; exit(1); }


print "\n\nmrunpairs.pl -- end.\n\n";





