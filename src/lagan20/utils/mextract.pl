#!/usr/bin/env perl

if (@ARGV < 1) {
    print ("usage:\n mextract.pl filename [-masked]\n");
    exit(1);
}

$masked=0;
$filename = $ARGV[0];
if(@ARGV==2) {
    if ($ARGV[1] eq "-masked") {
	$masked = 1;
    }
}

open(FASTAFILE, "$filename") || die "Could not open $filename.\n\n";
$prefix = substr $filename, 0, (rindex $filename, ".");
if ($masked || index ($filename, ".masked") != -1) {
    $prefix = substr $filename, 0, (rindex $prefix, ".");
}

$line = <FASTAFILE>;
chomp $line;

while (substr($line, 0, 1) ne ">") {
    $line = <FASTAFILE>;
    chomp $line;
}

$suffix = "fa";
if ($masked) {
    $suffix = "$suffix.masked";
}

if (substr($line, 0, 1) eq ">") {
    $name = substr($line, 1);
    if (index ($name, " ") != -1){
	$name = substr($name, 0, index ($name, " "));
    }
    if (substr ($name, length ($name) - 1) eq ","){
	$name = substr($name, 0, length ($name) - 1);
    }
#    $name = substr($line, 1);
#    $_ = substr($line, 1);
#    /\w+/g;
#    $name = $&;

#    substr($line, 1)." " =~ /(.+)[,]\s+/g;
#    $name = $1;

    $fname = "$prefix\_$name.$suffix";
    print("$fname\n");
    open(OUTFILE, ">$fname");
    print OUTFILE ">$name\n";
} else {
    print ("$filename is NOT a Multi-FASTA file...\n");
    exit(1);
}

while ($line = <FASTAFILE>) {
    chomp $line;
    if (substr($line, 0, 1) eq ">") {
	close OUTFILE;

#	substr($line, 1)." " =~ /(.+)[,]\s/g;
#	$name = $1;

	$name = substr($line, 1);
	if (index ($name, " ") != -1){
	    $name = substr($name, 0, index ($name, " "));
	}
	if (substr ($name, length ($name) - 1) eq ","){
	    $name = substr($name, 0, length ($name) - 1);
	}
#	$_ = substr($line, 1);
#	/\w+/g;
#	$name = $&;

	$fname = "$prefix\_$name.$suffix";
	print("$fname\n");
	open(OUTFILE, ">$fname");
	print OUTFILE ">$name\n";
    } else {
	print OUTFILE "$line";
    }
}

close OUTFILE;
