#!/usr/bin/env perl

# defaults 

$linelen = 50;
$interval = 10;
$labellen = 5;
$uselabels = 1;
$useintervals = 1;
$usecounts = 1;
$usebase = 0;
$liststart = 1;
$listend = 0;
$usestart = 0;
$useend = 0;

# constants

$minlinelen = 10;
$mininterval = 10;
$minlabellen = 3;


# usage notes

if (@ARGV < 1) {
    print ("usage:\n mpretty.pl filename\n");
    print ("options:\n");
    print (" -linelen value\n");
    print ("  (min: $minlinelen, default: $linelen)\n");
    print (" -interval value\n");
    print ("  (min: $mininterval, default: $interval, none: 0)\n");
    print (" -labellen value\n");
    print ("  (min: $labellen, default: $labellen, none: 0)\n");
    print (" -base sequence_name\n");
    print ("  (if used, must specify a sequence on which to base counting\n");
    print (" -start value\n");
    print ("  (if used, must specify a start coordinate (>=1)\n");
    print (" -end value\n");
    print ("  (if used, must specify an end coordinate (>=start)\n");
    print (" -nocounts\n");
    exit(1);
}


# parse parameters

for ($i=1; $i<@ARGV; $i++) {
    if ($ARGV[$i] eq "-nocounts") {
	$usecounts = 0;
    }
    if ($ARGV[$i] eq "-linelen") {
	$linelen = $ARGV[++$i];
	if ($linelen < $minlinelen) {
	    $linelen = $minlinelen;
	}
    }
    if ($ARGV[$i] eq "-interval") {
	$interval = $ARGV[++$i];
	if ($interval <= 0) {
	    $useintervals = 0;
	}
	if ($interval < $mininterval) {
	    $interval = $mininterval;
	}
    }
    if ($ARGV[$i] eq "-labellen") {
	$labellen = $ARGV[++$i];
	if ($labellen <= 0) {
	    $uselabels = 0;
	}
	if ($labellen < $minlabellen) {
	    $labellen = $minlabellen;
	}
    }
    if ($ARGV[$i] eq "-base") {
	$baseseq = $ARGV[++$i];
	$usebase = 1;
    }
    if ($ARGV[$i] eq "-start") {
	$usestart = 1;
	$liststart = $ARGV[++$i];
    }
    if ($ARGV[$i] eq "-end") {
	$useend = 1;
	$listend = $ARGV[++$i];
    }
}

# preprocessing for labels

if ($uselabels) {
    $labtail = "";
    for ($i=0; $i<$labellen; $i++) {
	$labtail="$labtail ";
    }
}

if (($usestart && ($liststart<1)) || ($useend && ($listend<$liststart))) {
    die "Invalid range specified: [$liststart, $listend].\n\n"; 
}

# read in Multi-FASTA file

$filename = $ARGV[0];
open(FASTAFILE, "$filename") || die "Could not open $filename.\n\n";
$line = <FASTAFILE>;
chomp $line;

while (substr($line, 0, 1) ne ">") {
    $line = <FASTAFILE>;
    chomp $line;
}

$i=0;
%list=();
@seqs=(());

if (substr($line, 0, 1) eq ">") {
    $_ = substr($line, 1);
    /\w+/g;
    @keys[$i] = $&;
    @count[$i]=0;
    @label[$i] = substr("@keys[$i]$labtail", 0, $labellen);
    $list{@keys[$i]}=$i;
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
	@count[$i]=0;
	@label[$i] = substr("@keys[$i]$labtail", 0, $labellen);
	$list{@keys[$i]}=$i;
	push @seqs, ();
    } else {
	push @{$seqs[$i]}, "$line";
    }
}

$i=0;
$maxlen = 0;
for $row (@seqs) {
    @strs[$i++] = join "", @$row;
    $templen = length @strs[$i-1];
    if ($templen > $maxlen) {
	$maxlen = $templen;
    }
}

$foundseq=0;
if ($usebase) {
    foreach $s (@keys) {
	$foundseq = ($s eq $baseseq) || $foundseq;
    }
if (!$foundseq) { die "Could not find Base Sequence: <$baseseq>\n\n"; }
}	

# preprocessing for counts

if ($usecounts) {
    foreach $s (@keys) {
	$_ = @strs[$list{$s}];
	$ls = tr/ATCGNatcgn/ATGCNatcgn/;
	@tot[$list{$s}] = $ls;
    }
}

# length of sequence display
$l=$maxlen; 
if ((!$listend) || ($listend>$maxlen)) {
    $listend = $maxlen;
}

if ($maxlen < $liststart) { die "Starting out of bounds...\b\b"; }


if ($usebase) {

# find base sequence position

    $i=0;
    $j=0;
    while ($j<$liststart) {
	if (substr(@strs[$list{$baseseq}], $i, 1) ne "-") {
	    $j++;
	}
	$i++;
    }
    $liststart = $i;
    while ($j<$listend) {
	if (substr(@strs[$list{$baseseq}], $i, 1) ne "-") {
	    $j++;
	}
	$i++;
    }
    $listend = $i;
}

# pretty print

if ($usecounts) {
    foreach $s (@keys) {
	$_ = substr(@strs[$list{$s}], 0, $liststart-1);
	$lc = tr/ATCGN/ATGCN/;
	@count[$list{$s}]+=$lc;
    }
}

for ($i=$liststart-1; $i<$listend; $i+=$linelen) {
    if ($listend-$i<$linelen) { $linelen = $listend-$i;}
    foreach $s (@keys) {
	if ($uselabels) {
	    print "@label[$list{$s}] : ";
	}
	$p = substr(@strs[$list{$s}], $i, $linelen);
	print "$p";
	    
	if ($usecounts) {
	    $_ = $p;
	    $lc = tr/ATCGN/ATGCN/;
	    @count[$list{$s}]+=$lc;
	    print " @ @count[$list{$s}]/@tot[$list{$s}]";
	}
	    
	print "\n";
    }
	
    if ($useintervals) {
	if ($uselabels) {
	    print "$labtail = ";
	}
	for ($j=$i+1; $j<=$i+$linelen && $j<=$l; $j+=$interval) {
	    $ct = "$j";
	    print $ct;
	    for ($k=0; $k<($interval-(length $ct)); $k++) {
		print " ";
	    }
	}
	print "\n";
    }
    print "\n";
}














