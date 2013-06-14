#!/usr/bin/env perl

$savname1 = "";
$savname2 = "";
$skip = 0;
$endblock = 0;
$score = 0;
$strand = "";
$initstrnd;
$s1s = 999999999;
$s2s = 999999999;
$first = 1;
$plus_sc = 0;
$minus_sc = 0;


while ($line = <STDIN>) {

    if ($line =~ /^>/) {
	if (!$first) { 
	    if ($strand eq "+") {
		print STDOUT "  Region [$s1s $reg1s][$s2s $reg2s] $score $strand\n";
	    }
	    else {
		print STDOUT "  Region [$s1s $reg1s][$reg2s $s2s] $score $strand\n";
	    }

	    if ($strand ne $initstrnd) {
		print STDOUT "INV\n"
	    }
	    if ($strand eq "+") { $plus_sc += $score; }
	    else  { $minus_sc += $score; }	  
	    if ($plus_sc > $minus_sc) {
		print STDOUT "Main score (+) $plus_sc; Inverted $minus_sc\n";
	    }
	    else {
		print STDOUT "Main score (-) $minus_sc; Inverted $plus_sc\n";
	    }
	    $plus_sc = 0;
	    $minus_sc = 0;
	    $score = 0;
	    $s1s = 999999999;
	    $s2s = 999999999;
	    $strand = "";
	}
	$first = 1;
	$name1 = $line;
	chomp $name1;
	$line = <STDIN>;
	if ($line !~ /^>/) {
	    print STDERR "Expecting a name, but got $line";
	    exit (1);
	}
	$name2 = $line;
	chomp $name2;
	$inblock = 1;
	$skip = 0;
	if (($name1 eq $savname1) && ($name2 eq $savname2)) {
	    $skip = 1;
	}
	else { 	print STDOUT "$name1 $name2\n"; }

	$savname1 = $name1;
	$savname2 = $name2;
    }
    elsif (!$skip) {
	$endblock = 0;
	$line =~ /\((\d+) (\d+)\)=\((\d+) (\d+)\) ([0-9\.]*) (.) (.*)/;
	if ($1 == 0 || $3 == 0) {
	    next;
	}
#	print STDOUT "strand $strand $s2s $4\n";
	if (($strand eq "+") && ($6 eq "+") && ($s2s + 20 < $4) ) {
	    $endblock += 2;
	}
	if (($strand eq "-") && ($6 eq "-") && ($s2s > $4 + 20) ) {
	    $endblock += 2;
	}
	if ($strand eq "") { $strand = $6; }
	if ($6 ne $strand) {
	    $endblock += 1;
	}

	if (!$endblock) {
	    $s2s = $3;  
	    $s1s = $1;
	    $s1e = $2;
	    $s2e = $4;
	    $score += $5;
	    if ($first) {
		print STDOUT "    "; 
		print STDOUT "    "; 
		$initstrnd = $strand;
		$reg1s = $2;
		$reg2s = $4;
		$first = 0;
	    }
	}
	else {
	    if ($strand eq "+") {
		print STDOUT "  Region [$s1s $reg1s][$s2s $reg2s] $score $strand\n";
	    }
	    else {
		print STDOUT "  Region [$s1s $reg1s][$reg2s $s2s] $score $strand\n";
	    }

	    if ($strand eq "+") { $plus_sc += $score; }
	    else  { $minus_sc += $score; }	  

	    if ($endblock %2) { print STDOUT "INV "; }
	    else {print STDOUT "    "; }
	    if ($endblock > 1) { print STDOUT "TRL "; }
	    else {print STDOUT "    "; }
	    $s2s = $3;  
	    $s1s = $1;
	    $s1e = $2;
	    $s2e = $4; 
	    $reg1s = $s1e; 
	    $reg2s = $s2e; 
	    $score = $5;
	    $strand = $6;
	    #	    print STDOUT "strand $strand\n";
	}
    }
}
if (!$first){
    if ($strand eq "+") {
	print STDOUT "  Region [$s1s $reg1s][$s2s $reg2s] $score $strand\n";
    }
    else {
	print STDOUT "  Region [$s1s $reg1s][$reg2s $s2s] $score $strand\n";
    }
    if ($strand eq "+") { $plus_sc += $score; }
    else  { $minus_sc += $score; }
}

if ($plus_sc > $minus_sc) {
    print STDOUT "Main score (+) $plus_sc; Inverted $minus_sc\n";
}
else {
    print STDOUT "Main score (-) $minus_sc; Inverted $plus_sc\n";
}
