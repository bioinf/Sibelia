#!/usr/bin/perl

use strict;

$0 = rindex($0, "/") > -1 ? substr($0, rindex($0, "/")+1) : $0;

my (@lines, @filt_lines);
my ($line, $line_in, $type);

my $mode = ($ARGV[0] eq "1" ? "M1" : ($ARGV[0] eq "2" ? "M2" : die("$0: Invalid base genome argument (expected 1 or 2)")));

die("$0: LAGAN_DIR not defined. Stopped") unless defined $ENV{"LAGAN_DIR"};

while (<STDIN>) {
	$line_in = $_;
	if ($line_in =~ /^\=.*(DM|M1|M2)$/) {
		$type = $1; $line .= $line_in;
		$lines[$#lines+1] = $line if $type eq "DM" or $type eq $mode;
		undef $line; undef $type;
	} else {
		$line .= $line_in;
	}
}

foreach my $line (@lines) {
	if ($mode eq "M2") {
		$line =~ /(\>[^\s\n]+\s([\+\-])[^\n]+)\n(.+)\n(\>[^\s\n]+\s([\+\-])[^\n]+)\n(.+)\n(\=.+?)\n/s;
#		$line =~ /(\>[^\s\n]+\s([\+\-])[^\n]+)\n([^\n]+)\n(\>[^\s\n]+\s([\+\-])[^\n]+)\n([^\n]+)\n(\=.+?\n)/s;
		
		my ($head1, $strand1, $seq1, $head2, $strand2, $seq2, $foot) = ($1, $2, $3, $4, $5, $6, $7);
		
		die if $strand1 ne $strand2;
		if ($strand1 eq "-") {
			$seq1 =~ s/\n//g;
			$seq2 =~ s/\n//g;
			$seq1 = reverse($seq1);
			$seq2 = reverse($seq2);
			$seq1 =~ s/(.{80})/$1\n/g;
			$seq2 =~ s/(.{80})/$1\n/g;
		}
		$line = $head2."\n".$seq2."\n".$head1."\n".$seq1."\n".$foot."\n";
	}
	push @filt_lines, $line;
}

open(OUT, "> tmp.xmfa");
foreach my $line (@filt_lines) { print OUT $line; }
close OUT;

system($ENV{"LAGAN_DIR"}."/utils/Glue tmp.xmfa > glue.out 2> glue.err");

open(IN, "< glue.out");
my @glue_out = <IN>;
close IN;

open(IN, "< glue.err");
my @glue_err = <IN>;
close IN;

unlink("tmp.xmfa");
unlink("glue.out");
unlink("glue.err");

print STDOUT @glue_out;
print STDERR @glue_err;
