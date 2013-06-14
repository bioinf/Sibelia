#!/usr/bin/perl

use strict;

$0 = rindex($0, "/") > -1 ? substr($0, rindex($0, "/")+1) : $0;

die("$0: LAGAN_DIR not defined. Stopped") unless defined $ENV{"LAGAN_DIR"};
my $LAGAN_DIR = $ENV{LAGAN_DIR};

my ($outfile, $base);

foreach my $arg (@ARGV) {
	if ($arg =~ /-out\s+([^\s]+)/) {
		$outfile = $1;
		$arg =~ s/-out\s+([^\s]+)//;
	} elsif ($arg =~ /-base[\s\=]+([^\s]+)/) {
		$base = $1;
		$arg =~ s/-base[\s\=]+([^\s]+)//;
		die("$0: Invalid base parameter (expected 1 or 2). Stopped") unless $base eq "1" or $base eq "2";
	}
}

if (@ARGV < 2) {
	print ("Usage:\n$0 seqfile1 seqfile2 [-glocal \"glocal flags\"] [-chaos \"chaos flags\"] [-order \"order flags\"] [-recurse \"(wl1,nd1,co1),(wl2,nd2,co2),...\"] [-mfa] [-out \"filename\"] [-maskedonly] [-debug] [-translate] [-fastreject]\n");
	exit(1);
}

my $args = join(" ", @ARGV);
system($LAGAN_DIR."/slagan.pl $args > slagan.pl.out");
die("$0: slagan.pl returned error $?. Stopped") if $?;

system($LAGAN_DIR."/xmfa2mfa.pl ".($base eq "2" ? "2" : "1")." < slagan.pl.out ".($outfile ? "> $outfile" : ""));
die("$0: xmfa2mfa.pl returned error $?. Stopped") if $?;

unlink "slagan.pl.out";
