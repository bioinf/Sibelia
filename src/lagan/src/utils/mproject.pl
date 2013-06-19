#!/usr/bin/env perl

if (@ARGV < 2) {
    print ("usage:\n mproject.pl filename seqname1 [seqname2 ... ]\n");
    exit(1);
}

$filename = $ARGV[0];

$i = 1;
while ($i < @ARGV) {
    @targets[$i-1] = $ARGV[$i];
    $i++;
}

open(FASTAFILE, "$filename") || die "Could not open $filename.\n\n";

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
	$_ = substr($line,1);
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

$seqlen = length $strs[0];
# $seqlen--;

for ($i=0; $i<$seqlen; $i++) {
    @isgap[$i] = 1;
    foreach $s (@targets) {
	if (substr(@strs[$list{$s}], $i, 1) ne "-") {
	    @isgap[$i] = 0;
	    break;
	}
    }
}

foreach $s (@targets) {
    print ">@keys[$list{$s}]\n";
    $j=0;
    for ($i=0; $i<$seqlen; $i++) {
	if(!@isgap[$i]) {
	    print substr(@strs[$list{$s}], $i, 1);
	    $j++;
	    if (($j % 60) == 0) {
		print "\n";
	    }
	}
    }
    print "\n";
} 










