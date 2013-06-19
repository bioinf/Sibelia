#!/usr/bin/env perl

# defaults 
# constants

# usage notes

if (@ARGV < 1) {
    print ("usage:\n mf2bin.pl inputfile [-out outputfile] \n");
    exit(1);
}

# parse parameters

$tofile = 0;
for ($i=1; $i<@ARGV; $i++) {
    if ($ARGV[$i] eq "-out") {
	$tofile = 1;
	$outfilename = $ARGV[++$i];
    }
}

if ($tofile) {
    open(OUTFILE, ">$outfilename");
}

# read in Multi-FASTA file

$infilename = $ARGV[0];
open(FASTAFILE, "$infilename") || die "Could not open $infilename.\n\n";
$line = <FASTAFILE>;
chomp $line;

$i=0;
%list=();
@seqs=(());

if (substr($line, 0, 1) eq ">") {
    $_ = substr($line, 1);
    /\w+/g;
    @keys[$i] = $&;
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
	$list{@keys[$i]}=$i;
	push @seqs, ();
    } else {
	push @{$seqs[$i]}, "$line";
    }
}

$i=0;
for $row (@seqs) {
    @strs[$i++] = join "", @$row;
}

if (@keys != 2) {
    print ("mpack needs two FASTA sequences\n");
    exit(1);
}


# pack bin
# format from Alex Poliakov's glass2bin.pl script

%base_code = ('-' => 0, 'A' => 1, 'C' => 2, 'T' => 3, 'G' => 4, 'N' => 5,
	      'a' => 1, 'c' => 2, 't' => 3, 'g' => 4, 'n' => 5);
$l = length @strs[0]; # $l--;
$s1 = reverse(@strs[0]);
$s2 = reverse(@strs[1]);


for ($i=0; $i<$l; $i++) {
    if ($tofile) {
	print OUTFILE pack("H2", 
			   $base_code{chop($s1)} . $base_code{chop($s2)});
    } else {
	print pack("H2", 
		   $base_code{chop($s1)} . $base_code{chop($s2)});
    }
}


