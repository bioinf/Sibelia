#!/usr/bin/perl

# Supermap: Piecewise monotonic alignment map generator for Shuffle-LAGAN
# Author: Andrey Kislyuk (kislyuk@ocf.berkeley.edu)

package Supermap;
require 5.005;
my ($VERSION) = ('$Id: supermap.pl,v 1.50 2005/06/15 22:40:04 kislyuk Exp $' =~ /,v\s+(\d+\S+)/o);

# Default constant values
my $overlap_factor = 0.8; # Aligns will be discarded if another align overlaps them by this factor or more in both seqs and has the same orientation
my $max_asym = 10; # Chains will be formed only if the resulting region's lengths differ by at most this factor
my $min_seq_score; # All aligns for sequences with this total score will be discarded. See getMinSeqScore
my $max_expand_len = 30000; # Aligns will be expanded or contracted on both sides on both strands by this amount up to the total length below
my $expand_factor = 4; # When one of an align's sequences is constrained in its expansion by a neighbor/start/end, the other one will be expanded by this times more than the first one
my $max_chainlen = 1500000; # Aligns will not be joined if the total length on either strand exceeds this. Set 0 to disable (no chain length limit)
my $max_job_size = 50000; # Maximum job size, in blat hits, for chunking when running glocal in parallel
my $erode_align = 15; # Amount by which to erode the coords of each align loaded (to avoid overlap problems when chaining)
my ($c1, $c2, $c3, $c4) = (100, 50, 400, 25); # BLAT->CHAOS score conversion parameters
#my $max_dist_y = 10000; # Join x-monotonic into same single-chain only if at most that apart in y-species.
my $default_lagan_dir = "/home/genome/glocal";
my $glocal_name = (0 ? "SLAGAN" : "glocal");

use Getopt::Long;
use File::Path;
use File::Copy;
use Cwd;
use IPC::Open2;
use IO::Handle;
#use Carp;
use strict;
use warnings;
no warnings "uninitialized";

sub main();
sub init();
sub getSeqSizes($$$);
sub prepareHits();
sub runSLAGAN();
sub reprintInputHits($$$);
sub processResults();
sub removeSLAGANOutput();
sub seqBelowMinScore($);
sub alignHashID($);
sub printChainToTemp($$$$);
sub chainBase1Hits($$);
sub chainBase2Hits($);
sub load2MHashes($);
sub loadBase2Hashes($);
sub postProcessRegions();
sub workerRun($$$$);
sub dequeueClustJobs($);
sub get_all_seqs($$);
sub isBLAT($);
sub useIf($$);
sub writeSizes($$);
sub getMinSeqScore($);
sub checkAlignCoords($);
sub expandSeq1($$);
sub expandSeq2($$);
sub finalExpand($$);
sub expSeq1Reg($$$$$);
sub expSeq2Reg($$$$$);
sub finalExpReg($$$$$);

# array index constants
use constant START1 =>  0; use constant END1   =>  1;
use constant START2 =>  2; use constant END2   =>  3;
use constant SEQ1   =>  4; use constant SEQ2   =>  5;
use constant ORIENT =>  6; use constant ORIGIN =>  7;
use constant SCORE  =>  8; use constant TOTSC  =>  9;
use constant HASHID => 10; use constant FLIPPED=> 11;
use constant CHALO1 => 12; use constant CHAHI1 => 13;
use constant CHALO2 => 14; use constant CHAHI2 => 15;
use constant CHALO1E=> 16; use constant CHAHI1E=> 17;
use constant CHALO2E=> 18; use constant CHAHI2E=> 19;
#use constant PREV1  =>  8; use constant NEXT1  =>  9;
#use constant PREV2  => 10; use constant NEXT2  => 11;
#use constant OSTART1=> 12; use constant OEND1  => 13;
#use constant OSTART2=> 14; use constant OEND2  => 15;

$SIG{'INT'} = $SIG{'QUIT'} = $SIG{'HUP'} = $SIG{'TRAP'} = $SIG{'ABRT'} = $SIG{'STOP'} = $SIG{'TERM'} = \&dequeueClustJobs;

my ($debug, $quiet, $outfile, $proflip, $skip, $no_pid, $input_glob, $input_dir,
	$server, $db, $gen1, $gen2, $gen1sizefile, $gen2sizefile, $write_sizes1, $write_sizes2,
	$score_file, $cfg, $cfg_file, $sizes1, $sizes2, $dbh, $tmp_dir, $tmp_prefix, $nodelete,
	$clust_run_pid, $print_chains, $no_aligntotals, $no_clust_run, $num_jobs, $input_is_blat,
	$force_overwrite, $print_csv, $using_GP, $slagan_params, $tmp_existed, $print_stats, $lagan_dir, $glocal_out_logfile);
my (@input_files);
my (%offsets1, %offsets2, %aligns1, %aligns2, %flipped_aligns);

my $supermapexec = $0; my $mycwd = getcwd(); $supermapexec =~ s/^\./$mycwd/ unless $supermapexec =~ /^\.\./; $supermapexec = $mycwd."/".$supermapexec if $supermapexec =~ /^\.\./;
die("$0: Problem resolving my name, \'$supermapexec\' is not a file") unless -f $supermapexec or $ARGV[0] eq "worker";
$0 = rindex($0, "/") > -1 ? substr($0, rindex($0, "/")+1) : $0;

$lagan_dir = $ENV{"LAGAN_DIR"} if defined $ENV{"LAGAN_DIR"};
$lagan_dir = $ENV{"LAGAN_DIR"} = $default_lagan_dir unless defined $ENV{"LAGAN_DIR"};
$lagan_dir =~ s/^\.\./$mycwd\/\.\./;
$lagan_dir =~ s/^\./$mycwd\//;
$ENV{"LAGAN_DIR"} = $lagan_dir;
print STDERR "$0: Warning: LAGAN_DIR=$lagan_dir is not a valid directory\n" unless -d $lagan_dir;
push @INC, $lagan_dir;

my $SLAGAN = $lagan_dir."/".$glocal_name;
my $error_file = "./$0.$$.error.log";
my $default_score_file = $lagan_dir."/test.score";
my $default_outfile = "$0.out";
my $worker_tmp_dir = "/tmp/$0.$$.worker/"; # The directory where workers store their intermediate files (two workers should not use the same directory)

my $usage = "
-infile=file \t Name of input file containing all hits for the two genomes
-outfile=file \t Output filename (default: $default_outfile)
-gen1=id \t First genome ID (must exist in the GPDB)
-gen2=id \t Second genome ID (must exist in the GPDB)
-sizes1=file \t File with sequence sizes for first genome
-sizes2=file \t File with sequence sizes for second genome
-bacteria \t Rearrange circular DNA to find a better alignment map
-server=hostname GPDB server (default: lemur)
-db=dbname \t GPDB name (default: GP)
-config=file \t GPDB config file (default: ~/.gprc)
-score=file \t Score file for SLAGAN (default: $default_score_file)
-glocal_out=file \t Save intermediate GLOCAL alignment hits to this file
-no_clust_run \t Run CPU/memory intensive jobs locally, not on the GP cluster
-tmp_dir=dir \t Working directory (default: /tmp/$0.pid)
-f \t\t Overwrite output file without prompting if it exists
-v \t\t Verbose mode
-q \t\t Quiet mode
-k \t\t Keep all temporary files
-expand_length=N Maximum length by which to expand alignments (default: $max_expand_len)
-max_length=N \t Maximum length for any alignment chain in either strand
\t\t (default: $max_chainlen)
-min_seq_score=N Sequences with total align score below this threshold will be
\t\t discarded (default: U penalty in SLAGAN score file)
-max_job_size=N  Threshold, in hits, for splitting workload into separate jobs
\t\t for clust_run (default: $max_job_size)
-c1, c2, c3, c4=N: Score factors for BLAT->CHAOS conversion
\t\t (default: $c1, $c2, $c3, $c4)

Options may be abbreviated.
Input file format is BLAT or CHAOS. Sequence names should not contain spaces.
Alignments with negative scores are discarded.
Sequence size file format, one sequence per line: seq_name seq_size
";

exit(main());

# ___ Subroutines _______________

sub main() {
	if ($ARGV[0] eq "worker") { workerRun($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]); exit(0); } # Running SLAGAN in distributed mode
	init();

	print("$0: Retrieving sequence info...\n") unless $quiet;
	$sizes1 = getSeqSizes($dbh, $gen1, $gen1sizefile);
	(writeSizes($sizes1, $write_sizes1), exit(0)) if defined $write_sizes1;
	$sizes2 = getSeqSizes($dbh, $gen2, $gen2sizefile);
	(writeSizes($sizes2, $write_sizes2), exit(0)) if defined $write_sizes2;

	die("$0: No sequence size data found. Stopped") if (keys(%$sizes1) < 1 or keys(%$sizes2) < 1);
	die("$0: Flip mode is only applicable for two single-sequence organisms. Stopped") if ($proflip and not (keys(%$sizes1) == 1 and keys(%$sizes2) == 1));

	# Sort and separate the alignments, run SLAGAN on them
	prepareHits();
	runSLAGAN();

	# Chain SLAGAN alignments into supermonotonic chain and save the intermediate results
	my ($dc, $sc1, $sc2) = processResults();

	# Load the results back and expand regions, then print them
	postProcessRegions();

	print "$0: Output written to $outfile\n" unless $quiet;
	print "$0: Intermediate files kept in $tmp_dir\n" if $nodelete and not $quiet;
	rmdir $tmp_dir unless $tmp_existed or $nodelete;

	return 0;
}


# Startup tasks
sub init() {
	system('export LC_ALL="C"'); # Things may misbehave if locale is set to UTF-8

	# Berkeley Genome Pipeline functionality is used if corresponding Perl modules are found in @INC
	foreach my $dir (@INC) {
		$using_GP = 1 if -f $dir."/GPDBI.pm" and -f $dir."/GPutils.pm";
	}

	useIf $using_GP, "GPDBI";
	useIf $using_GP, "GPutils";
	useIf 1, "Utils";
#	useIf 1, "Desoverlap";

	die("$0: GetOptions failed to retrieve options. Check the input options. Usage:".$usage) unless
	GetOptions(
		"server=s"		=> \$server,
		"gen1=s"		=> \$gen1,
		"gen2=s"		=> \$gen2,
		"sizes1=s"		=> \$gen1sizefile,
		"sizes2=s"		=> \$gen2sizefile,
		"blatfile=s"	=> \$input_glob,
		"infile=s"		=> \$input_glob,
		"outfile=s"		=> \$outfile,
		"glocal_out=s"	=> \$glocal_out_logfile,
		"bacteria"		=> \$proflip,
		"server=s"		=> \$server,
		"db=s"			=> \$db,
		"config=s"		=> \$cfg_file,
		"tmp_dir=s"		=> \$tmp_dir,
		"skip"			=> \$skip,
		"no_pid"		=> \$no_pid,
		"no_clust_run"	=> \$no_clust_run,
		"print_chains"	=> \$print_chains,
		"print_stats"	=> \$print_stats,
		"no_aligntotals"=> \$no_aligntotals,
		"print_csv"		=> \$print_csv,
		"max_job_size"	=> \$max_job_size,
		"max_length=i"	=> \$max_chainlen,
		"expand_length=i"=>\$max_expand_len,
		"min_seq_score=i"=>\$min_seq_score,
		"max_asym=i"	=> \$max_asym,
		"overlap_factor"=> \$overlap_factor,
		"score=s"		=> \$score_file,
		"c1=i"			=> \$c1,
		"c2=i"			=> \$c2,
		"c3=i"			=> \$c3,
		"c4=i"			=> \$c4,
		"slagan_params"	=> \$slagan_params,
		"write_sizes1=s"=> \$write_sizes1,
		"write_sizes2=s"=> \$write_sizes2,
		"keep"			=> \$nodelete,
		"f"				=> \$force_overwrite,
		"v"				=> \$debug,
		"q"				=> \$quiet
	);

	undef $quiet if $debug;
	my @uinfo = getpwuid($>);
	print("$0: Version ".$VERSION." started ".localtime()." by ".$uinfo[0]."\n") unless $quiet;
	$tmp_prefix = $0.($no_pid ? "" : ".".$$);

	unless ($no_clust_run) {
		$no_clust_run = `which clust_run 2> /dev/null`; $no_clust_run = not $no_clust_run;
		print("$0: clust_run not found - cluster operation disabled\n") if $no_clust_run and not $quiet;
	}

	if ($tmp_dir) {
		$tmp_existed = 1 if -d $tmp_dir;
		mkdir $tmp_dir unless -d $tmp_dir;
		$tmp_dir .= "/" unless /\/^Z/;
	} else {
		$tmp_dir = "/tmp/".$tmp_prefix;
		mkdir $tmp_dir;
		$tmp_dir .= "/";
	}
	die("$0: No write permissions in working directory $tmp_dir. Stopped") unless -w $tmp_dir;
	die("$0: Genome IDs or size files not specified. Usage:".$usage) unless ($gen1 or $gen1sizefile) and ($gen2 or $gen2sizefile);
	die("$0: '-gen' options are invalid because GPDB is not available. Use '-sizes'. Stopped") if (($gen1 or $gen2) and not $using_GP);
	die("$0: Sequence size file $gen1sizefile not found. Stopped") unless -f $gen1sizefile or $gen1;
	die("$0: Sequence size file $gen2sizefile not found. Stopped") unless -f $gen2sizefile or $gen2;
	die("$0: Maximum job size too small, must exceed 10000 hits. Stopped") if $max_job_size < 10000;
	die("$0: Overlap factor must be between 0 and 1. Stopped") if $overlap_factor < 0 or $overlap_factor > 1;
	print("$0: SLAGAN score file not specified, using default $default_score_file\n") unless $score_file or $quiet;
	print("$0: Output file not specified, using default $default_outfile\n") unless $outfile or $quiet;

	# Check input file or glob
	if (defined $input_glob) {
		if ($input_glob =~ /\//) { ($input_dir, $input_glob) = ($input_glob =~ /\A(.*\/)([^\/]+)\Z/); }
		$input_glob .= "\$" unless $input_glob =~ /\$$/;
		$input_glob = "^".$input_glob unless $input_glob =~ /^\^/;
		@input_files = Utils::safe_glob($input_glob, $input_dir);
	} elsif (@ARGV > 0) {
		foreach my $file (@ARGV) {
			if ($file =~ /\//) { ($input_dir, $file) = ($file =~ /\A(.*\/)([^\/]+)\Z/); }
			push @input_files, $file;
		}
	} else { # TODO: split stdin for >2GB input
		open(FH, "> $tmp_dir$tmp_prefix.in");
		print FH while <STDIN>;
		close FH;
		push @input_files, "$tmp_prefix.in";
		$input_dir = $tmp_dir;
	}
	unless ($input_dir =~ /\A\//) { $input_dir = $mycwd."/".$input_dir; }
	die("$0: No input files matching \"$input_dir$input_glob\" found. Stopped") unless @input_files > 0;
	print "$0: ".@input_files." input file(s)\n" if $debug;

	# Check output file
	$outfile = $default_outfile unless $outfile;
	if (-f $outfile and not $force_overwrite and -t STDERR) {
		print STDERR "$0: $outfile exists. Overwrite? (y/N, '-f' to force) ";
		my $overwrite = <STDIN>; chomp $overwrite;
		(print("Move \"$outfile\" or use option '-f'.\n"), exit(1)) unless ($overwrite eq "Y" or $overwrite eq "y" or $overwrite eq "yes");
	}
	open(FH, "> ".$outfile) or die("$0: Cannot open $outfile for writing: $!");
	close FH;

	# Check SLAGAN score file
	$score_file = $default_score_file unless $score_file;
	unless ($score_file =~ /\A\//) { $score_file = $mycwd."/".$score_file; }
	$max_expand_len += $erode_align;
	die("$0: max_length cannot be less than 0. Stopped") if $max_chainlen < 0;
	$max_chainlen = 1000000000 if $max_chainlen == 0;
	$max_chainlen -= 2*$max_expand_len;
	# SLAGAN output for a given sequence will be discarded if the total score for the sequence is below this threshold. Default value is the SLAGAN unrelated gap penalty.
	$min_seq_score = getMinSeqScore($score_file) unless defined $min_seq_score;

	# Connect to GPDB
	if ($using_GP) {
		$GPutils::Error = "";
		$cfg = read_gp_config(Get_Abs_Path($cfg_file)) or die($GPutils::Error);
		$server ||= $cfg->Get_Val("DB", "server");
		$db ||= $cfg->Get_Val("DB", "main_db");
		$dbh = GPDBI->connect($server, 0, $db, undef, undef, "gp_cgi", undef, {PrintError => 0, RaiseError => 1});
	}
}


# Load sequence names and sizes either from GPDB or from file
sub getSeqSizes($$$) {
	my ($dbh, $dataset, $gen_size_file) = @_;
	if ($dataset) {
		return get_all_seqs($dbh, $dataset);
	} else {
		my %sizes;
		open(FH, "< ".$gen_size_file) or die("$0: Could not open file $gen_size_file for reading: ".$!);
		while (<FH>) {
			chomp;
			my ($seq, $size) = split;
			die("$0: Invalid format in file $gen_size_file") unless $seq and $size;
			$sizes{$seq} = $size;
		}
		close FH;
		return \%sizes;
	}
}


# Convert BLAT to CHAOS if necessary
# Flip hits on circular sequence if necessary
sub prepareHits() {
	my ($cur_align);
	local (*FH, *OUT1);

	print "$0: Preparing files...\n" unless $quiet;
	$input_is_blat = 1 if isBLAT($input_dir.$input_files[0]);

	if ($input_is_blat) {
		foreach my $file (@input_files) {
			system('awk \'{$13=($13+$15)?$13:1; print $1,$2,$3";",$5,$6,$7"; '.
				'score = "' . $c1 . '*$8-' . $c2 . '*$9-' . $c3 . '*($12+$14)-' . $c4 .
				'*log($13+$15),"("$4")"}\''.
				"< $input_dir$file > $tmp_dir$file.chaos");
		}
	} else {
		foreach my $file (@input_files) {
			system('ln -s "'.$input_dir.$file.'" "'.$tmp_dir.$file.'.chaos"');
		}
	}

	if ($proflip) {
		open(FH, "< ".$tmp_dir.$input_files[0].".chaos") or die("$0: Could not open file ".$tmp_dir.$input_files[0].".chaos for reading: ".$!);
		open(OUT1, "> ".$tmp_dir.$input_files[0].".flipped.chaos") or die("$0: Could not open file ".$tmp_dir.$input_files[0].".flipped.chaos for writing: ".$!);

		my (@seq1s, @seq1e, @seq2s, @seq2e, @scores, @orientations, @seqn1, @seqn2);
		my ($seq1center, $seq2center, $seq1median, $seq2median);
		my $i = 0;
		while (<FH>) {
			/\A[\s]*.*\s([\d]+)\s([\d]+)\;\s.*\s([\d]+)\s([\d]+)\;\sscore\s\=\s([e\d\.\+\-]+)\s\(([\+\-]+)\)/;
#			($seqn1[$i], $seq1s[$i], $seq1e[$i], $seqn2[$i], $seq2s[$i], $seq2e[$i], $scores[$i], $orientations[$i]) = ($1, $2, $3, $4, $5, $6, $7, $8);
			($seq1s[$i], $seq1e[$i], $seq2s[$i], $seq2e[$i], $scores[$i], $orientations[$i]) = ($1, $2, $3, $4, $5, $6);
			if ($seq1s[$i] > $seq1e[$i]) { my $j = $seq1s[$i]; $seq1s[$i] = $seq1e[$i]; $seq1e[$i] = $j; }
			if ($seq2s[$i] > $seq2e[$i]) { my $j = $seq2s[$i]; $seq2s[$i] = $seq2e[$i]; $seq2e[$i] = $j; }
			$i++;
		}

		# For each interval pair,
		# if the seq1 interval median is greater than seq1 median, and the corresponding interval median in seq2 is less than seq2 median,
		# OR if the seq1 interval median is less than seq1 median, and the corresponding interval median in seq2 is greater than seq2 median,
		# set start of interval in seq1 to 2CoM1 - previous end of interval
		# set end of interval in seq1 to 2CoM1 - previous start of interval
		# flip the orientation (+/-)
		$seq1center = $$sizes1{(keys(%$sizes1))[0]} / 2;
		$seq2center = $$sizes2{(keys(%$sizes2))[0]} / 2;
		my $flip_counter = 0;
		foreach $i (0..@seq1s-1) {
			$seq1median = ($seq1s[$i] + $seq1e[$i]) / 2;
			$seq2median = ($seq2s[$i] + $seq2e[$i]) / 2;
			if (($seq1median > $seq1center and $seq2median < $seq2center)
				or ($seq1median < $seq1center and $seq2median > $seq2center)) {
				my $j = $seq2s[$i];
				$seq2s[$i] = (2 * $seq2center) - $seq2e[$i];
				$seq2e[$i] = (2 * $seq2center) - $j;
				if ($orientations[$i] eq "+") { $orientations[$i] = "-"; } else { $orientations[$i] = "+"; }
				$cur_align = [];
				$$cur_align[START1] = $seq1s[$i]; $$cur_align[START2] = $seq2s[$i];
				$$cur_align[END1] = $seq1e[$i]; $$cur_align[END2] = $seq2e[$i];
				$$cur_align[SCORE] = $scores[$i]; $$cur_align[ORIENT] = $orientations[$i];
$$cur_align[SEQ1] = (keys(%$sizes1))[0]; $$cur_align[SEQ2] = (keys(%$sizes2))[0];
$$cur_align[START1] += $erode_align; $$cur_align[END1] -= $erode_align;
$$cur_align[START2] += $erode_align; $$cur_align[END2] -= $erode_align;
				$flipped_aligns{alignHashID($cur_align)} = $cur_align;
				$flip_counter++;
			}
			print OUT1 "seq1 ".$seq1s[$i]." ".$seq1e[$i]."; seq2 ".$seq2s[$i]." ".$seq2e[$i]."; score = ".$scores[$i]." (".$orientations[$i].")\n";
		}
		close FH; close OUT1;
		print "$0: Single-sequence flip mode: ".($flip_counter+0)." hits flipped\n" if $debug;
	}
}


# Load all hits into a hash table, then write the hits for each sequence into a file
# Run SLAGAN on each of these files, via worker instances either on the cluster or sequentially
sub runSLAGAN() {
	my ($clust_run_invoke, $num_jobs, $sort_pid1, $sort_pid2, $sort_pid3, $one_seq_mode,
		$cur_align, $next_align, $curlen1, $curlen2, $nextlen1, $nextlen2, $overlap1, $overlap2, $dump_count);
	local (*RH1, *WH1, *RH2, *WH2, *RH3, *WH3, *IN, *DUPES);
#	my $filter = Desoverlap->new($overlap_factor, $debug);

	print "$0: Sorting input hits...\n" if $debug;
	open(DUPES, "> supermap.duplicates") if $debug;

	$one_seq_mode = 1 if (keys(%$sizes1) == 1 and keys(%$sizes2) == 1);

	$sort_pid1 = open2(\*RH1, \*WH1, "sort --key=1,1 --key=2,2n"); # pre-scan
	$sort_pid2 = open2(\*RH2, \*WH2, "sort --key=1,1 --key=2,2n"); # gen1base
	$sort_pid3 = open2(\*RH3, \*WH3, "sort --key=4,4 --key=5,5n"); # gen2base

	# Sort input on seq1
	foreach my $file (@input_files) {
		open(IN, "< $tmp_dir$file".($proflip?".flipped":"").".chaos");
		print WH1 while <IN>;
		close IN;
	}
	close WH1;

	# Scan input, check if start2, end2 are ascending for sorting, erode alignments
	while (<RH1>) {
		/\A[\s]*(.*)\s([\d]+)\s([\d]+)\;\s(.*)\s([\d]+)\s([\d]+)\;\sscore\s\=\s([e\d\.\+\-]+)\s\(([\+\-]+)\)/o;

		$next_align=[];
		($$next_align[SEQ1], $$next_align[START1], $$next_align[END1], $$next_align[SEQ2], $$next_align[START2], $$next_align[END2], $$next_align[SCORE], $$next_align[ORIENT])
		= ($1, $2, $3, $4, $5, $6, $7, $8);
		next if $$next_align[SCORE] <= 0;
		if ($one_seq_mode) { $$next_align[SEQ1] = (keys(%$sizes1))[0]; $$next_align[SEQ2] = (keys(%$sizes2))[0]; }
		checkAlignCoords($next_align);

		unless ($$next_align[END1]-$$next_align[START1] <= $erode_align*2 or $$next_align[END2]-$$next_align[START2] <= $erode_align*2) {
				$$next_align[START1] += $erode_align; $$next_align[END1] -= $erode_align;
				$$next_align[START2] += $erode_align; $$next_align[END2] -= $erode_align;
		}

=head1
		# Overlap scan
		if ($$next_align[START1] <= $$cur_align[END1] and $$next_align[END1] >= $$cur_align[START1] # overlap in seq1
		and $$next_align[START2] <= $$cur_align[END2] and $$next_align[END2] >= $$cur_align[START2] # overlap in seq2
		and $$cur_align[SEQ1] eq $$next_align[SEQ1] and $$cur_align[SEQ2] eq $$next_align[SEQ2]
		and $$cur_align[ORIENT] eq $$next_align[ORIENT]) {
			($curlen1, $curlen2, $nextlen1, $nextlen2)
				= ($$cur_align[END1] - $$cur_align[START1] + 1, $$cur_align[END2] - $$cur_align[START2] + 1,
				   $$next_align[END1] - $$next_align[START1] + 1, $$next_align[END2] - $$next_align[START2] + 1);

			if ($$next_align[START1] <= $$cur_align[START1] and $$next_align[END1] >= $$cur_align[END1]) {
				$overlap1 = $$cur_align[END1] - $$cur_align[START1] + 1; # next covers cur
			} elsif ($$next_align[START1] <= $$cur_align[START1]) {
				$overlap1 = $$next_align[END1] - $$cur_align[START1] + 1; # next is to the left
			} elsif ($$next_align[END1] >= $$cur_align[END1]) {
				$overlap1 = $$cur_align[END1] - $$next_align[START1] + 1; # next is to the right
			} else {
				$overlap1 = $$next_align[END1] - $$next_align[START1] + 1; # cur covers next
			}
			if ($$next_align[START2] <= $$cur_align[START2] and $$next_align[END2] >= $$cur_align[END2]) {
				$overlap2 = $$cur_align[END2] - $$cur_align[START2] + 1;
			} elsif ($$next_align[START2] <= $$cur_align[START2]) {
				$overlap2 = $$next_align[END2] - $$cur_align[START2] + 1;
			} elsif ($$next_align[END2] >= $$cur_align[END2]) {
				$overlap2 = $$cur_align[END2] - $$next_align[START2] + 1;
			} else {
				$overlap2 = $$next_align[END2] - $$next_align[START2] + 1;
			}
			die("$0: Bad internal state") if $overlap1 < 0 or $overlap2 < 0;

			if (($overlap1 / $curlen1 > $overlap_factor) and ($overlap2 / $curlen2 > $overlap_factor)
			and $$cur_align[SCORE] <= $$next_align[SCORE]) {
				$dump_count++;
				print DUPES "Cur: (".$$cur_align[START1]."-".$$cur_align[END1].")(".$$cur_align[START2]."-".$$cur_align[END2].") ".$$cur_align[SCORE]." over with (".$$next_align[START1]."-".$$next_align[END1].")(".$$next_align[START2]."-".$$next_align[END2].") ".$$next_align[SCORE]."\n" if $debug;
				$cur_align = $next_align; next; # discard current align
			} elsif (($overlap1 / $nextlen1 > $overlap_factor) and ($overlap2 / $nextlen2 > $overlap_factor)
			and $$cur_align[SCORE] >= $$next_align[SCORE]) {
				$dump_count++;
				print DUPES "Nxt: (".$$next_align[START1]."-".$$next_align[END1].")(".$$next_align[START2]."-".$$next_align[END2].") ".$$next_align[SCORE]." over with (".$$cur_align[START1]."-".$$cur_align[END1].")(".$$cur_align[START2]."-".$$cur_align[END2].") ".$$cur_align[SCORE]."\n" if $debug;
				next; # discard next align
			}
		}
=cut
		foreach my $cur_align ($next_align){ # (@{$filter->put($next_align)}) {
			print WH2 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n";
			print WH3 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n";
		}
		
#		print WH2 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n" if @$cur_align;
#		print WH3 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n" if @$cur_align;
#		$cur_align = $next_align;
	}
#	$filter->printAll();
	# Flush alignments remaining in filter buffer
#	foreach my $cur_align (@{$filter->getBuffer()}) {
#		print WH2 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n" if $cur_align != 0;
#		print WH3 $$cur_align[SEQ1]." ".$$cur_align[START1]." ".$$cur_align[END1]."; ".$$cur_align[SEQ2]." ".$$cur_align[START2]." ".$$cur_align[END2]."; "."score = ".$$cur_align[SCORE]." (".$$cur_align[ORIENT].")\n" if $cur_align != 0;
#	}

	close RH1; waitpid $sort_pid1, 0;

	close WH2;
	$num_jobs = reprintInputHits(1, 1, \*RH2);
	close RH2; waitpid $sort_pid2, 0;

	close WH3;
	$num_jobs = reprintInputHits(2, $num_jobs, \*RH3);
	close RH3; waitpid $sort_pid3, 0;

	close DUPES if defined fileno DUPES;
#	print STDERR "$0: Warning: ".$filter->{dump_count}." near duplicate alignments discarded (overlap factor $overlap_factor)\n" if $filter->{dump_count} and not $quiet;

	open(FH, "> ".$tmp_dir."CLUSTER_JOB_PARAMS") or die;
	foreach my $i (1..$num_jobs-1) {
		print FH "worker JOB".$i.".tar ".$score_file." ".$SLAGAN.($debug ? " -v" : "");
		print FH " << JOB$i.tar > CLUSTER_JOB_MESSAGES.$i >> CLUSTER_JOB_ERRMSG.$i" unless $no_clust_run;
		print FH "\n";
	}
	close FH;

	if ($no_clust_run) {
		open(FH, "< ".$tmp_dir."CLUSTER_JOB_PARAMS") or die;
		print "$0: Running ".($num_jobs-1)." SLAGAN jobs locally...\n" unless $quiet;
		while (<FH>) {
			chomp;
			print("Job $.: \"$0 $_\"\n") if $debug;
			system("cd $tmp_dir; $supermapexec ".$_);
		}
		close FH;
	} else {
		$clust_run_invoke = "clust_run -program=".$supermapexec." -parameters=".$tmp_dir."CLUSTER_JOB_PARAMS -init_dir=$tmp_dir -wait";
		print "$0: Running ".($num_jobs-1)." distributed SLAGAN jobs with clust_run...\n" unless $quiet;
		print "$0: \"$clust_run_invoke\"\n" if $debug;

		if ($clust_run_pid = fork()) { # I am the parent
			waitpid($clust_run_pid, 0);
		} elsif (not defined $clust_run_pid) {
			die("$0: Could not fork");
		} else { # I am the child
			die("$0: Could not exec \"$clust_run_invoke\"") unless exec($clust_run_invoke);
		}
		undef $clust_run_pid;
	}

	foreach my $i (1..$num_jobs-1) {
		system("cd $tmp_dir; tar -xf ".$tmp_dir."JOB".$i.".results.tar");
		unlink $tmp_dir."JOB".$i.".tar" unless $nodelete;
		unlink $tmp_dir."JOB".$i.".results.tar" unless $nodelete;
		unlink $tmp_dir."CLUSTER_JOB_MESSAGES.$i" unless $nodelete;
		unlink $tmp_dir."CLUSTER_JOB_ERRMSG.$i" unless $nodelete;
	}

	unlink "$tmp_dir$input_glob.chaos" unless $nodelete;
	unlink $tmp_dir."CLUSTER_JOB_PARAMS" unless $nodelete;
	
	foreach my $file (@input_files) {
		unlink $tmp_dir.$file.".chaos" unless $nodelete;
	}
}


sub reprintInputHit($$$) {
	my ($base_gen, $align, $FH) = @_;
	if ($base_gen == 1 and $$align[ORIENT] eq "+") {
		print $FH $$align[SEQ1]." ".$$align[START1]." ".$$align[END1]."; ".$$align[SEQ2]." ".$$align[START2]." ".$$align[END2]."; "."score = ".$$align[SCORE]." (".$$align[ORIENT].")\n";
	} elsif ($base_gen == 1 and $$align[ORIENT] eq "-") {
		print $FH $$align[SEQ1]." ".$$align[START1]." ".$$align[END1]."; ".$$align[SEQ2]." ".$$align[END2]." ".$$align[START2]."; "."score = ".$$align[SCORE]." (".$$align[ORIENT].")\n";
	} elsif ($base_gen == 2 and $$align[ORIENT] eq "+") {
		print $FH $$align[SEQ2]." ".$$align[START2]." ".$$align[END2]."; ".$$align[SEQ1]." ".$$align[START1]." ".$$align[END1]."; "."score = ".$$align[SCORE]." (".$$align[ORIENT].")\n";
	} elsif ($base_gen == 2 and $$align[ORIENT] eq "-") {
		print $FH $$align[SEQ2]." ".$$align[START2]." ".$$align[END2]."; ".$$align[SEQ1]." ".$$align[END1]." ".$$align[START1]."; "."score = ".$$align[SCORE]." (".$$align[ORIENT].")\n";
	} else {
		die("$0: Bad internal state from hit ".$$align[SEQ1]." ".$$align[START1]." ".$$align[END1]."; ".$$align[SEQ2]." ".$$align[START2]." ".$$align[END2]."; "."score = ".$$align[SCORE]." (".$$align[ORIENT].")");
	}
}


sub writeJobFile($$) {
	my ($job_id, $seq_list) = @_;
	local *LIST;

	open(LIST, "| cd $tmp_dir; xargs tar --append --file=".$tmp_dir."JOB".$job_id.".tar");
	foreach my $file (sort alnum keys(%$seq_list)) { $file =~ /\/([^\/]+)$/; print LIST $1." "; }
	close LIST;

	foreach my $file (sort alnum keys(%$seq_list)) { unlink $file unless $nodelete; }
}


# Separate input into files based on sequence name and reverse order in gen2base hits
sub reprintInputHits($$$) {
	my ($base_gen, $job_id, $RH) = @_;
	my ($one_seq_mode, $line_count, $prev_seq, $cur_seq, $cur_align);
	my (%cur_seq_list, %pruned_sizes);
	local (*OUT, *LIST);

	$one_seq_mode = 1 if (keys(%$sizes1) == 1 and keys(%$sizes2) == 1);

	print "$0: Reprinting hits (base genome $base_gen)..." if $debug;

	$line_count = 0;
	while (<$RH>) {
		/\A[\s]*(.*)\s([\d]+)\s([\d]+)\;\s(.*)\s([\d]+)\s([\d]+)\;\sscore\s\=\s([e\d\.\+\-]+)\s\(([\+\-]+)\)/o;

		$cur_align=[];
		($$cur_align[SEQ1], $$cur_align[START1], $$cur_align[END1], $$cur_align[SEQ2], $$cur_align[START2], $$cur_align[END2], $$cur_align[SCORE], $$cur_align[ORIENT])
			= ($1, $2, $3, $4, $5, $6, $7, $8);

		$cur_seq = ($base_gen == 1 ? $$cur_align[SEQ1] : $$cur_align[SEQ2]);

		if ($cur_seq ne $prev_seq) {
			$pruned_sizes{$cur_seq} = ($base_gen == 1 ? $$sizes1{$cur_seq} : $$sizes2{$cur_seq});
			print " ".$cur_seq if $debug;
			close OUT if defined fileno OUT;
			open(OUT, "> ".$tmp_dir.$input_files[0].".gen".$base_gen."base.".$cur_seq.".chaos") or die("$0: Could not open file ".$tmp_dir.$input_files[0].".gen".$base_gen."base.".$cur_seq.".chaos for writing: ".$!);
			if ($line_count > $max_job_size) {
				writeJobFile($job_id, \%cur_seq_list);
				undef %cur_seq_list; $line_count = 0; $job_id++;
			}
			$cur_seq_list{$tmp_dir.$input_files[0].".gen".$base_gen."base.".$cur_seq.".chaos"} = 1;
		}
		reprintInputHit($base_gen, $cur_align, \*OUT) if @$cur_align;

		$prev_seq = $cur_seq;
#		$cur_align = $next_align;
		$line_count++;
	}

#	reprintInputHit($base_gen, $next_align, \*OUT) if @$next_align;
	writeJobFile($job_id, \%cur_seq_list);
	$job_id++;

	close OUT;
	print "\n" if $debug;
	$sizes1 = \%pruned_sizes if $base_gen == 1;
	$sizes2 = \%pruned_sizes if $base_gen == 2;
	return $job_id;
}


sub seqBelowMinScore($) {
	my ($line) = @_;
	$line =~ /\A[\s]*\([\d]+\s[\d]+\)\=\([\d]+\s[\d]+\)\s([\d\.\-]+)\s[\+\-]+\s\[([\d\.\-]+)\][\s]*s1\:.*[\s]*s2\:.*\n\Z/;
	die("$0: Unable to extract score values from SLAGAN output:\n$line") if not defined $2;
	return ($2 < $min_seq_score);
}

sub processResults() {
	my ($cur_seq, $input_prefix, $dropped_seqs, $sort_pid, $sort_pid2);
	local (*RH, *WH, *IN, *OUT, *hashesDM_RH, *hashesDM_WH);
	print "$0: Loading SLAGAN output...\n" unless $quiet;
	open(GLOCAL_OUT_LOG, "> ".$glocal_out_logfile) if $glocal_out_logfile;

	# Sort gen2base aligns on seq1, then seq2, then start2, then print them to separate files, one file per gen1 seq
	# These files will be loaded on demand when scanning gen1base aligns (chainBase1Hits())
	$sort_pid = open2(\*RH, \*WH, "sort --key=9,9 --key=7,7 --key=1.2,1n"); # input is base 2, key is 9 because a space is expected between s2: and seq2name
	$input_prefix = $tmp_dir.$input_files[0].".gen2base";
	foreach my $seq (sort alnum keys(%$sizes2)) {
		open(IN, "< $input_prefix.$seq.chaos.glocal-out") or (delete($$sizes2{$seq}), next);
		my $line = <IN>;
		die("$0: Empty SLAGAN output file $input_prefix.$seq.chaos.glocal-out, check corresponding job logs. Stopped") unless $line;
		if (seqBelowMinScore($line)) { print "$0: Discarding file $input_prefix.$seq.chaos.glocal-out - score too low ($1<$min_seq_score)\n" if $debug; next; }
		seek IN, 0, 0; # back to start
		print WH while <IN>;
		close IN;
	}
	close WH or die("$0: Error executing sort");
	while (<RH>) {
		/\ss2\:[\s]*([^\s]+)[\s]*\n\Z/;
		if ($1 ne $cur_seq or not defined $cur_seq) {
			next unless $1;
			close OUT if defined fileno OUT;
			$cur_seq = $1;
			open(OUT, "> $input_prefix.sorted-gen1.$cur_seq.chaos.glocal-out") or die("$0: Could not open file $input_prefix.sorted-gen1.$cur_seq.chaos.glocal-out for writing: ".$!);
		}
		print OUT $_;
	}
	close RH; close OUT if defined fileno OUT;
	waitpid $sort_pid, 0;

	# Sort gen1base aligns on seq1, then start1
	$sort_pid = open2(\*RH, \*WH, "sort --key=7,7 --key=1.2,1n"); # input is base 1
	$input_prefix = $tmp_dir.$input_files[0].".gen1base";
	foreach my $seq (sort alnum keys(%$sizes1)) {
		open(IN, "< $input_prefix.$seq.chaos.glocal-out") or (delete($$sizes1{$seq}), next);
		my $line = <IN>;
		if (seqBelowMinScore($line)) { $dropped_seqs++; print "$0: Discarding file $input_prefix.$seq.chaos.glocal-out - score too low ($1<$min_seq_score)\n" if $debug; next; }
		seek IN, 0, 0; # back to start
		print WH while <IN>;
		if ($glocal_out_logfile) { seek IN, 0, 0; print GLOCAL_OUT_LOG while <IN>; }
		close IN;
		unlink "$input_prefix.$seq.chaos.glocal-out" unless $nodelete;
	}
	unlink $input_prefix.".chaos" unless $nodelete;
	close WH or die("$0: Error executing sort");

	# Feed the gen1base aligns to the 2M/1M1 chain scanner (chainBase1Hits())
	# The hashesDM handle is used to write 2M aligns' hashes to be sorted in seq2 order
	print "$0: Generating supermonotonic map...\n" unless $quiet;
	$sort_pid2 = open2(\*hashesDM_RH, \*hashesDM_WH, "sort --key=2,2");
	chainBase1Hits(*RH, *hashesDM_WH);
	close RH;
	waitpid $sort_pid, 0;
	close hashesDM_WH or die("$0: Error executing sort");

	# Print sorted 2M aligns' hashes, one file per gen2 seq
	undef $cur_seq;
	while(<hashesDM_RH>) {
		my $line = $_;
		$line =~ /\A[^\s]+\s([^\s]+)\s[^\s]+\n\Z/;
		if ($1 ne $cur_seq or not defined $cur_seq) {
			close OUT if defined fileno OUT;
			$cur_seq = $1;
			open(OUT, "> $tmp_dir".$input_files[0].".hashesDM.gen2.$cur_seq") or die("$0: Could not open file $tmp_dir".$input_files[0].".hashesDM.gen2.$cur_seq for writing: ".$!);
		}
		print OUT $line;
	}
	close hashesDM_RH;
	waitpid $sort_pid2, 0;

	# Sort gen2base aligns on seq2, then start2
	$sort_pid = open2(\*RH, \*WH, "sort --key=7,7 --key=1.2,1n"); # input is base 2
	$input_prefix = $tmp_dir.$input_files[0].".gen2base";
	foreach my $seq (sort alnum keys(%$sizes2)) {
		open(IN, "< $input_prefix.$seq.chaos.glocal-out") or next;
		my $line = <IN>;
		if (seqBelowMinScore($line)) { $dropped_seqs++; print "$0: Discarding file $input_prefix.$seq.chaos.glocal-out - score too low ($1<$min_seq_score)\n" if $debug; next; }
		seek IN, 0, 0; # back to start
		print WH while <IN>;
		close IN;
		unlink "$input_prefix.$seq.chaos.glocal-out" unless $nodelete;
	}
	unlink $input_prefix.".chaos" unless $nodelete;
	close WH or die("$0: Error executing sort");

	# Feed the gen2base aligns to the 1M2 chain scanner (chainBase2Hits())
	chainBase2Hits(*RH);
	close RH;
	waitpid $sort_pid, 0;

	close GLOCAL_OUT_LOG if defined fileno GLOCAL_OUT_LOG;

	removeSLAGANOutput();
	print STDERR "$0: Warning: Alignments for $dropped_seqs sequences discarded due to total score below cutoff ($min_seq_score)\n" if $dropped_seqs and not $quiet;
}


sub removeSLAGANOutput() {
	my $input_prefix = $tmp_dir.$input_files[0].".gen1base";
	foreach my $seq (sort alnum keys(%$sizes1)) { unlink "$input_prefix.$seq.chaos.glocal-out" unless $nodelete; }
	unlink $input_prefix.".chaos" unless $nodelete;

	$input_prefix = $tmp_dir.$input_files[0].".gen2base";
	foreach my $seq (sort alnum keys(%$sizes2)) { unlink "$input_prefix.$seq.chaos.glocal-out" unless $nodelete; }
	unlink $input_prefix.".chaos" unless $nodelete;

	rmdir $tmp_dir;
}


sub alignHashID($) {
	my ($align) = @_;
#	return 23*$$align[START1] + 41*$$align[START2] + 61*$$align[END1] + 83*$$align[END2];
	return $$align[SEQ1].":".$$align[START1]."-".$$align[END1]."=".$$align[SEQ2].":".$$align[START2]."-".$$align[END2];
}


# The chain writer lags the chainer by two chains because the full contents of neighboring chains must be known.
sub printChainToTemp($$$$) {
	my ($FH, $prev_chain, $cur_chain, $next_chain) = @_;
	return unless defined $cur_chain;

	my $type = ${$$cur_chain[0]}[ORIGIN];
	my ($first_align, $last_align) = ($$cur_chain[0], $$cur_chain[@$cur_chain-1]);
	print $FH ${$$cur_chain[0]}[ORIGIN]." ".@$cur_chain." ".
		$$first_align[START1]." ".$$first_align[END1]." ".$$first_align[START2]." ".$$first_align[END2]." ".
		$$first_align[SEQ1]." ".$$first_align[SEQ2]." ".$$first_align[ORIENT]." ".$$first_align[SCORE]." ".
		$$last_align[START1]." ".$$last_align[END1]." ".$$last_align[START2]." ".$$last_align[END2]." ".
		$$last_align[SEQ1]." ".$$last_align[SEQ2]." ".$$last_align[ORIENT]." ".$$last_align[SCORE];
	if ($print_chains) {
		foreach my $align (@$cur_chain) {
			print $FH " ".$$align[START1]." ".$$align[END1]." ".$$align[START2]." ".$$align[END2];
		}
	}
	print $FH "\n";
}


sub chainBase1Hits($$) {
	my ($FH, $hashesDM) = @_;
	local *OUT;
	my ($cur_align, $prev_align, $cur_chain, $prev_chain, $pre_prev_chain, $chain_start_2M, $chain_start_1M1,
		$cur_seq, $align_peers, $flip_counter);
	my @bad_aligns; my %base2peers;

	while (<$FH>) {
		/\A[\s]*\(([\d]+)\s([\d]+)\)\=\(([\d]+)\s([\d]+)\)\s([\d\.\-]+)\s([\+\-]+)\s\[([\d\.\-]+)\][\s]*s1\:(.*)[\s]*s2\:(.*)\n\Z/;

		next if ($1==$2); # skip null alignments
		(push(@bad_aligns, $_), next) unless $1 and $2 and $3 and $4 and $5 and $6;

		$cur_align = [];
		($$cur_align[START1], $$cur_align[END1], $$cur_align[START2], $$cur_align[END2], $$cur_align[SCORE], $$cur_align[ORIENT], $$cur_align[TOTSC], $$cur_align[SEQ1], $$cur_align[SEQ2])
			= ($1, $2, $3, $4, $5, $6, $7, $8, $9);
		$$cur_align[SEQ1] =~ s/^\s+//; $$cur_align[SEQ1] =~ s/\s+$//;
		$$cur_align[SEQ2] =~ s/^\s+//; $$cur_align[SEQ2] =~ s/\s+$//;
#warn("Seen: ".$_) if $$cur_align[SEQ1] eq "AC002301.1";
		checkAlignCoords($cur_align);
				
		if ($proflip and defined $flipped_aligns{alignHashID($cur_align)}) {
			my $seq2center = $$sizes2{(keys(%$sizes2))[0]} / 2;
			my $j = $$cur_align[START2];
			$$cur_align[START2] = (2 * $seq2center) - $$cur_align[END2];
			$$cur_align[END2] = (2 * $seq2center) - $j;
			if ($$cur_align[ORIENT] eq "+") { $$cur_align[ORIENT] = "-"; } else { $$cur_align[ORIENT] = "+"; }
			$$cur_align[FLIPPED]=1;
			$flip_counter++;
		}

		$$cur_align[HASHID] = alignHashID($cur_align);

		if ($$cur_align[SEQ1] ne $cur_seq) {
#warn("Handling seq trans") if $prev_align and $$prev_align[SEQ1] eq "AC002301.1";
printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain);# unless defined $cur_seq;
printChainToTemp(*OUT, $prev_chain, $cur_chain, undef);# unless defined $cur_seq;

			undef $chain_start_2M; undef $chain_start_1M1; undef $prev_align;
			undef $pre_prev_chain; undef $prev_chain; undef $cur_chain;
			$cur_seq = $$cur_align[SEQ1];
			%base2peers = %{loadBase2Hashes($tmp_dir.$input_files[0].".gen2base.sorted-gen1.$cur_seq.chaos.glocal-out")};
			close OUT if defined fileno OUT;
			open(OUT, "> ".$tmp_dir.$input_files[0].".2MM1.$cur_seq");
		}

		$align_peers = $base2peers{$$cur_align[HASHID]};
		$$cur_align[ORIGIN] = defined($align_peers) ? 2 : 1;

		if ($chain_start_2M and defined $align_peers and defined $prev_align # continue open 2M chain
			and (($$cur_align[ORIENT] eq "+" and $$cur_align[START2] > $$prev_align[END2]
						and $$prev_align[HASHID] eq $$align_peers[0])
					or ($$cur_align[ORIENT] eq "-" and $$cur_align[END2] < $$prev_align[START2]
						and $$prev_align[HASHID] eq $$align_peers[1])
				or ($$cur_align[FLIPPED] and ($$cur_align[ORIENT] eq "+" and $$cur_align[START2] < $$prev_align[END2]
						and $$prev_align[HASHID] eq $$align_peers[0])
					or ($$cur_align[ORIENT] eq "-" and $$cur_align[END2] > $$prev_align[START2]
						and $$prev_align[HASHID] eq $$align_peers[1])))
			and $$cur_align[ORIENT] eq $$prev_align[ORIENT]
			and $$cur_align[FLIPPED] eq $$prev_align[FLIPPED]
			and $$cur_align[SEQ2] eq $$prev_align[SEQ2]
			and ($$cur_align[START1] > $$prev_align[END1] or ($$cur_align[FLIPPED] and $$cur_align[START1] > $$prev_align[END1]))
			and abs($$cur_align[END1] - $$chain_start_2M[START1]) < $max_chainlen
			and abs($$cur_align[END2] - $$chain_start_2M[START2]) < $max_chainlen
#and abs($$cur_align[END1] - $$chain_start_2M[START1])/abs($$cur_align[END2] - $$chain_start_2M[START2]) < $max_asym
#and abs($$cur_align[END2] - $$chain_start_2M[START2])/abs($$cur_align[END1] - $$chain_start_2M[START1]) < $max_asym
			) {
				push(@$cur_chain, $cur_align);
				print $hashesDM $$cur_align[SEQ1]."\t".$$cur_align[SEQ2]."\t".$$cur_align[HASHID]."\n";
		} elsif (defined $align_peers) { # start new 2M chain
			printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain);
			$chain_start_2M = $cur_align; undef $chain_start_1M1;
			$pre_prev_chain = $prev_chain; $prev_chain = $cur_chain;
			$cur_chain = [$cur_align];
			print $hashesDM $$cur_align[SEQ1]."\t".$$cur_align[SEQ2]."\t".$$cur_align[HASHID]."\n";
		} elsif ($chain_start_1M1 and defined $prev_align # continue open 1M1 chain
			and ((($$cur_align[ORIENT] eq "+" and $$cur_align[START2] > $$prev_align[END2])
					or ($$cur_align[ORIENT] eq "-" and $$cur_align[END2] < $$prev_align[START2]))
				or ($$cur_align[FLIPPED] and (($$cur_align[ORIENT] eq "+" and $$cur_align[START2] < $$prev_align[END2])
					or ($$cur_align[ORIENT] eq "-" and $$cur_align[END2] > $$prev_align[START2]))))
			and $$cur_align[ORIENT] eq $$prev_align[ORIENT]
			and $$cur_align[FLIPPED] eq $$prev_align[FLIPPED]
			and $$cur_align[SEQ2] eq $$prev_align[SEQ2]
			and ($$cur_align[START1] > $$prev_align[END1] or ($$cur_align[FLIPPED] and $$cur_align[START1] > $$prev_align[END1]))
			and abs($$cur_align[END1] - $$chain_start_1M1[START1]) < $max_chainlen
			and abs($$cur_align[END2] - $$chain_start_1M1[START2]) < $max_chainlen
#and abs($$cur_align[END1] - $$chain_start_1M1[START1])/abs($$cur_align[END2] - $$chain_start_1M1[START2]) < $max_asym
#and abs($$cur_align[END2] - $$chain_start_1M1[START2])/abs($$cur_align[END1] - $$chain_start_1M1[START1]) < $max_asym
			) {
				push(@$cur_chain, $cur_align);
		} else { # start new 1M1 chain
			printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain);
			$chain_start_1M1 = $cur_align; undef $chain_start_2M;
			$pre_prev_chain = $prev_chain; $prev_chain = $cur_chain;
			$cur_chain = [$cur_align];
		}
		$prev_align = $cur_align;
	}
	printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain);
	printChainToTemp(*OUT, $prev_chain, $cur_chain, undef);
	print "$0: Single-sequence flip mode: ".($flip_counter+0)." gen1base hits backflipped\n" if $debug and $proflip;
	warn "$0: Warning: ".@bad_aligns." bad SLAGAN alignments discarded" if @bad_aligns > 0;
}


# Input is base 2, i.e. (start2 end2)=(start1 end1)...
sub chainBase2Hits($) {
	my ($FH) = @_;
	local *OUT;
	my ($cur_align, $prev_align, $cur_chain, $prev_chain, $pre_prev_chain, $chain_start_2M, $chain_start_1M2,
		$cur_seq, $align_is_2M, $flip_counter);
	my @bad_aligns; my %aligns2M;

	while(<$FH>) {
		/\A[\s]*\(([\d]+)\s([\d]+)\)\=\(([\d]+)\s([\d]+)\)\s([\d\.\-]+)\s([\+\-]+)\s\[([\d\.\-]+)\][\s]*s1\:(.*)[\s]*s2\:(.*)\n\Z/;

		next if ($1==$2); # skip null alignments
		(push(@bad_aligns, $_), next) unless $1 and $2 and $3 and $4 and $5 and $6;

		$cur_align = [];
		($$cur_align[START2], $$cur_align[END2], $$cur_align[START1], $$cur_align[END1], $$cur_align[SCORE], $$cur_align[ORIENT], $$cur_align[TOTSC], $$cur_align[SEQ2], $$cur_align[SEQ1])
			= ($1, $2, $3, $4, $5, $6, $7, $8, $9);
		$$cur_align[SEQ1] =~ s/^\s+//; $$cur_align[SEQ1] =~ s/\s+$//;
		$$cur_align[SEQ2] =~ s/^\s+//; $$cur_align[SEQ2] =~ s/\s+$//;
		checkAlignCoords($cur_align);

		if ($proflip and defined $flipped_aligns{alignHashID($cur_align)}) {
			my $seq2center = $$sizes2{(keys(%$sizes2))[0]} / 2;
			my $j = $$cur_align[START2];
			$$cur_align[START2] = (2 * $seq2center) - $$cur_align[END2];
			$$cur_align[END2] = (2 * $seq2center) - $j;
			if ($$cur_align[ORIENT] eq "+") { $$cur_align[ORIENT] = "-"; } else { $$cur_align[ORIENT] = "+"; }
			$$cur_align[FLIPPED] = 1;
			$flip_counter++;
		}

		$$cur_align[HASHID] = alignHashID($cur_align);

		if ($$cur_align[SEQ2] ne $cur_seq) {
			printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain) if $$prev_chain[0][ORIGIN] == 3;# and not defined $cur_seq;
			printChainToTemp(*OUT, $prev_chain, $cur_chain, undef) if $$cur_chain[0][ORIGIN] == 3;# and not defined $cur_seq;
			undef $chain_start_1M2; undef $prev_align;
			undef $pre_prev_chain; undef $prev_chain; undef $cur_chain;
			$cur_seq = $$cur_align[SEQ2];
			%aligns2M = %{load2MHashes($tmp_dir.$input_files[0].".hashesDM.gen2.$cur_seq")};
			close OUT if defined fileno OUT;
			open(OUT, "> ".$tmp_dir.$input_files[0].".M2.$cur_seq");
		}
		$$cur_align[ORIGIN] = defined($aligns2M{$$cur_align[HASHID]}) ? 2 : 3;

		if (defined $aligns2M{$$cur_align[HASHID]}) { # align is 2M
			my $prev_ch_last_al = $prev_chain ? $$prev_chain[scalar(@$prev_chain)-1] : [];
			printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain) if $$prev_chain[0][ORIGIN] == 3;
			undef $chain_start_1M2; # close 1M2 chain
			$chain_start_2M = $cur_align;
			$pre_prev_chain = $prev_chain; $prev_chain = $cur_chain;
			$cur_chain = [$cur_align];
		} elsif ($chain_start_1M2 # continue open 1M2 chain
			and ((($$cur_align[ORIENT] eq "+" and $$cur_align[START1] > $$prev_align[END1])
					or ($$cur_align[ORIENT] eq "-" and $$cur_align[END1] < $$prev_align[START1]))
				or ($$cur_align[FLIPPED] and (($$cur_align[ORIENT] eq "+" and $$cur_align[START1] < $$prev_align[END1])
											or ($$cur_align[ORIENT] eq "-" and $$cur_align[END1] > $$prev_align[START1]))))
			and $$cur_align[ORIENT] eq $$prev_align[ORIENT]
			and $$cur_align[SEQ1] eq $$prev_align[SEQ1]
			and $$cur_align[FLIPPED] == $$prev_align[FLIPPED]
			and ($$cur_align[START2] > $$prev_align[END2] or ($$cur_align[FLIPPED] and $$cur_align[START2] < $$prev_align[END2]))
			and abs($$cur_align[END1] - $$chain_start_1M2[START1]) < $max_chainlen
			and abs($$cur_align[END2] - $$chain_start_1M2[START2]) < $max_chainlen
#and abs($$cur_align[END1] - $$chain_start_1M2[START1])/abs($$cur_align[END2] - $$chain_start_1M2[START2]) < $max_asym
#and abs($$cur_align[END2] - $$chain_start_1M2[START2])/abs($$cur_align[END1] - $$chain_start_1M2[START1]) < $max_asym
			) {
				push(@$cur_chain, $cur_align);
		} else { # start new 1M2 chain
			my $prev_ch_last_al = $prev_chain ? $$prev_chain[scalar(@$prev_chain)-1] : [];
			printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain) if $$prev_chain[0][ORIGIN] == 3;
			$chain_start_1M2 = $cur_align;
			$pre_prev_chain = $prev_chain; $prev_chain = $cur_chain;
			$cur_chain = [$cur_align];
		}
		$prev_align = $cur_align;
	}
	my $prev_ch_last_al = $prev_chain ? $$prev_chain[scalar(@$prev_chain)-1] : [];
	printChainToTemp(*OUT, $pre_prev_chain, $prev_chain, $cur_chain) if $$prev_chain[0][ORIGIN] == 3;
	printChainToTemp(*OUT, $prev_chain, $cur_chain, undef) if $$cur_chain[0][ORIGIN] == 3;
	print "$0: Single-sequence flip mode: ".($flip_counter+0)." gen2base hits backflipped\n" if $debug and $proflip;
	warn "$0: Warning: ".@bad_aligns." bad SLAGAN alignments discarded" if @bad_aligns > 0;
}


# Input: file with lines of the form "seq1 seq2 hash" (seq2 should be the same per file)
# Output: hash(key->align hash ID, value->1). Input file is deleted.
sub load2MHashes($) {
	my ($file) = @_;
	my %hashes;
	local *FH;
	open(FH, "< $file") or return {};
	while (<FH>) {
		/\A[^\s]+\t[^\s]+\t([^\s]+)\n\Z/;
		warn("Hash collision in \"$_\" vs. \"".$hashes{$1}."\"") if defined $hashes{$1};
		$hashes{$1} = 1;
	}
	close FH;
	unlink $file unless $nodelete;
	return \%hashes;
}


# Input: file with gen2base alignments which should have the same seq1 ordered by start2 or not exist
# Output: hash(key->align hash ID, value->[prev align hash ID, next align hash ID]). Input file is deleted.
# Input is base 2, i.e. (start2 end2)=(start1 end1)...
sub loadBase2Hashes($) {
	my ($file) = @_;
	my ($prev_align, $cur_align, $next_align);
	my %hashes;
	local *FH;
	open(FH, "< $file") or return {};
	while (<FH>) { # Scan 1 line ahead because the next align must also be seen
		/\A[\s]*\(([\d]+)\s([\d]+)\)\=\(([\d]+)\s([\d]+)\)\s.*s1\:(.*)[\s]*s2\:(.*)/;

		$next_align = [];
		# Hits are gen2base
		($$next_align[START2], $$next_align[END2], $$next_align[START1], $$next_align[END1], $$next_align[SEQ2], $$next_align[SEQ1]) = ($1, $2, $3, $4, $5, $6);
		checkAlignCoords($next_align);
		$$next_align[SEQ1] =~ s/^\s+//; $$next_align[SEQ1] =~ s/\s+$//;
		$$next_align[SEQ2] =~ s/^\s+//; $$next_align[SEQ2] =~ s/\s+$//;
		$$next_align[HASHID] = alignHashID($next_align);
		warn("LB2H: Hash collision in \"$_\"") if defined $cur_align and defined $hashes{$$cur_align[HASHID]};
		$hashes{$$cur_align[HASHID]} =
			[$prev_align ? $$prev_align[HASHID] : 1,
			 $next_align ? $$next_align[HASHID] : 1] if $cur_align;
		$prev_align = $cur_align; $cur_align = $next_align;
	}
	$hashes{$$cur_align[HASHID]} = [$prev_align ? $$prev_align[HASHID] : 1, undef] if $cur_align;
	close FH;
	unlink $file unless $nodelete;
	return \%hashes;
}


# Load chained regions and expand them according to the expansion rules, then print them out and display some chain statistics
sub postProcessRegions() {
	local (*IN, *OUT, *RH1, *WH1, *RH2, *WH2, *RH3, *WH3);
	my ($first_align, $last_align, $type, $num_aligns, $sort_pid1, $sort_pid2, $sort_pid3);
	my (@line, @min_lengths, @max_lengths, @means, @pos_counts, @neg_counts);

	$sort_pid1 = open2(\*RH1, \*WH1, "sort --key=7,7 --key=3,3n"); # sort on seq1, start1
	$sort_pid2 = open2(\*RH2, \*WH2, "sort --key=8,8 --key=5,5n"); # sort on seq2, start2
	$sort_pid3 = open2(\*RH3, \*WH3, "sort --key=7,7 --key=3,3n"); # sort on seq1, start1
#	open(WH1, "> ".$outfile) or die("$0: Could not open output file $outfile for writing: ".$!);

	open(OUT, "> ".$outfile) or die("$0: Could not open output file $outfile for writing: ".$!);
#	open(OUT, "| sort --key=1,1 --key=2,2n > ".$outfile) or die("$0: Could not open output file $outfile for writing: ".$!);
	foreach my $seq (sort alnum keys %$sizes1) {
		open(IN, "< ".$tmp_dir.$input_files[0].".2MM1.$seq") or next;
		print WH1 while <IN>;
		close IN;
		unlink $tmp_dir.$input_files[0].".2MM1.$seq" unless $nodelete;
	}

	foreach my $seq (sort alnum keys %$sizes2) {
		open(IN, "< ".$tmp_dir.$input_files[0].".M2.$seq") or next;
		print WH1 while <IN>;
		close IN;
		unlink $tmp_dir.$input_files[0].".M2.$seq" unless $nodelete;
	}

	close WH1;
	expandSeq1(\*RH1, \*WH2);
	close RH1; waitpid $sort_pid1, 0;
	close WH2;
	expandSeq2(\*RH2, \*WH3);
	close RH2; waitpid $sort_pid2, 0;
	close WH3;
	finalExpand(\*RH3, \*OUT);
	close RH3; waitpid $sort_pid3, 0;
	close OUT;
}


# Input: chains ordered by seq1, start1
# Output: chains expanded on seq1
sub expandSeq1($$) {
	my ($RH, $WH) = @_;
	my ($first_align, $last_align, $type, $num_aligns,
		$cur_seq, $preexpand1, $postexpand1,
		$prev_chain, $cur_chain, $next_chain);
	my (@line);

	while (<$RH>) {
		chomp; @line = split;

		# skip M2 regions
		if ($line[0] == 3) {
			$,= " "; print $WH @line[0..17]; print $WH " 0 0 0 0 "; print $WH @line[18..$#line]; print $WH "\n"; undef $,; next;
		}

		$prev_chain = $cur_chain;
		$cur_chain = $next_chain;

		$first_align = []; $last_align = [];
		($type, $num_aligns, $$first_align[START1], $$first_align[END1], $$first_align[START2], $$first_align[END2],
		$$first_align[SEQ1], $$first_align[SEQ2],$$first_align[ORIENT], $$first_align[SCORE],
		$$last_align[START1], $$last_align[END1], $$last_align[START2], $$last_align[END2],
		$$last_align[SEQ1], $$last_align[SEQ2], $$last_align[ORIENT], $$last_align[SCORE]) = @line;

		$$first_align[CHALO1] = ($$first_align[START1] < $$last_align[START1] ? $$first_align[START1] : $$last_align[START1]);
		$$first_align[CHAHI1] = ($$first_align[END1] > $$last_align[END1] ? $$first_align[END1] : $$last_align[END1]);

		my @saved_line = @line;
		$next_chain = [$first_align, $last_align, $type, $num_aligns, \@saved_line];
		next unless defined $cur_chain;

		expSeq1Reg($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq);
# TODO
#		if ($cur_seq ne $$first_align[SEQ1]) {
#			undef $cur_chain;
#			$cur_seq = $$first_align[SEQ1];
#		}
	}
	expSeq1Reg($WH, $cur_chain, $next_chain, undef, $cur_seq);
}


sub expSeq1Reg($$$$$) {
	my ($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq) = @_;
	my ($preexpand1, $postexpand1);

	$preexpand1 = $$cur_chain[0][CHALO1] - (defined $prev_chain ? $$prev_chain[0][CHAHI1] : 0);
	$preexpand1 = $max_expand_len if $preexpand1 > $max_expand_len;
#$preexpand1 = 0 if $preexpand1 < 0;
	$preexpand1 = $max_expand_len if $preexpand1 < 0; # !!!
	$postexpand1 = $$next_chain[0][CHALO1] - $$cur_chain[0][CHAHI1];
	$postexpand1 = $max_expand_len if $postexpand1 > $max_expand_len;
#$postexpand1 = 0 if $postexpand1 < 0;
	$postexpand1 = $max_expand_len if $postexpand1 < 0;
#$postexpand1 = 0 if defined $prev_chain and $$prev_chain[0][CHAHI1] > $$cur_chain[0][CHAHI1]; # don't expand if covered by another align
	$$cur_chain[0][CHALO1E] = $$cur_chain[0][CHALO1] - $preexpand1;
	$$cur_chain[0][CHALO1E] = 1 if $$cur_chain[0][CHALO1E] < 1;
	$$cur_chain[0][CHAHI1E] = $$cur_chain[0][CHAHI1] + $postexpand1;
	$$cur_chain[0][CHAHI1E] = $$sizes1{$$cur_chain[0][SEQ1]} if $$cur_chain[0][CHAHI1E] > $$sizes1{$$cur_chain[0][SEQ1]};

	$cur_seq = $$cur_chain[0][SEQ1] if not defined $cur_seq;
	if ($cur_seq ne $$cur_chain[0][SEQ1]) { # Correct upper expansion
		$$cur_chain[0][CHAHI1E] = $$cur_chain[0][CHAHI1] + $max_expand_len;
		$$cur_chain[0][CHAHI1E] = $$sizes1{$$cur_chain[0][SEQ1]} if $$cur_chain[0][CHAHI1E] > $$sizes1{$$cur_chain[0][SEQ1]};
	}

	print $WH $$cur_chain[2]." ".$$cur_chain[3]." ".
		$$cur_chain[0][START1]." ".$$cur_chain[0][END1]." ".$$cur_chain[0][START2]." ".$$cur_chain[0][END2]." ".
		$$cur_chain[0][SEQ1]." ".$$cur_chain[0][SEQ2]." ".$$cur_chain[0][ORIENT]." ".$$cur_chain[0][SCORE]." ".
		$$cur_chain[1][START1]." ".$$cur_chain[1][END1]." ".$$cur_chain[1][START2]." ".$$cur_chain[1][END2]." ".
		$$cur_chain[1][SEQ1]." ".$$cur_chain[1][SEQ2]." ".$$cur_chain[1][ORIENT]." ".$$cur_chain[1][SCORE]." ".
		$$cur_chain[0][CHALO1]." ".$$cur_chain[0][CHAHI1]." ".$$cur_chain[0][CHALO1E]." ".$$cur_chain[0][CHAHI1E];

	if ($print_chains) {
		my $i = 18;
		while (1) {
			print $WH " ".${$$cur_chain[4]}[$i]." ".${$$cur_chain[4]}[$i+1]." ".${$$cur_chain[4]}[$i+2]." ".${$$cur_chain[4]}[$i+3];
			last if @{$$cur_chain[4]} <= $i+4;
			$i+=4;
		}
	}
	print $WH "\n";
}


# Input: chains ordered by seq2, start2
# Output: chains expanded on seq1 and seq2 (final output)
sub expandSeq2($$) {
	my ($RH, $WH) = @_;
	my ($first_align, $last_align, $type, $num_aligns,
		$cur_seq, $preexpand1, $postexpand1, $preexpand2, $postexpand2,
		$prev_chain, $cur_chain, $next_chain);
	my (@line);

	while (<$RH>) {
		chomp; @line = split;

		# skip M1 regions
		if ($line[0] == 1) {
			$,= " "; print $WH @line[0..21]; print $WH " 0 0 0 0 "; print $WH @line[22..$#line]; print $WH "\n"; undef $,; next;
		}

		$prev_chain = $cur_chain;
		$cur_chain = $next_chain;

		$first_align = []; $last_align = [];
		($type, $num_aligns, $$first_align[START1], $$first_align[END1], $$first_align[START2], $$first_align[END2],
		$$first_align[SEQ1], $$first_align[SEQ2],$$first_align[ORIENT], $$first_align[SCORE],
		$$last_align[START1], $$last_align[END1], $$last_align[START2], $$last_align[END2],
		$$last_align[SEQ1], $$last_align[SEQ2], $$last_align[ORIENT], $$last_align[SCORE],
		$$first_align[CHALO1], $$first_align[CHAHI1], $$first_align[CHALO1E], $$first_align[CHAHI1E]) = @line;

		$$first_align[CHALO2] = ($$first_align[START2] < $$last_align[START2] ? $$first_align[START2] : $$last_align[START2]);
		$$first_align[CHAHI2] = ($$first_align[END2] > $$last_align[END2] ? $$first_align[END2] : $$last_align[END2]);

		my @saved_line = @line;
		$next_chain = [$first_align, $last_align, $type, $num_aligns, \@saved_line];

		next unless defined $cur_chain;
		expSeq2Reg($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq);
#		if ($cur_seq ne $$first_align[SEQ2]) {
#			undef $cur_chain;
#			$cur_seq = $$first_align[SEQ2];
#		}
	}
	expSeq2Reg($WH, $cur_chain, $next_chain, undef, $cur_seq);
}


sub expSeq2Reg($$$$$) {
	my ($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq) = @_;
	my ($preexpand1, $postexpand1, $preexpand2, $postexpand2);

	$preexpand1 = $$cur_chain[0][CHALO1] - $$cur_chain[0][CHALO1E];
	$postexpand1 = $$cur_chain[0][CHAHI1E] - $$cur_chain[0][CHAHI1];

	$preexpand2 = $$cur_chain[0][CHALO2] - (defined $prev_chain ? $$prev_chain[0][CHAHI2] : 0);
	$preexpand2 = $preexpand1 * $expand_factor if $preexpand2 > $preexpand1 * $expand_factor and $$cur_chain[2] != 3;
	$preexpand2 = $max_expand_len if $preexpand2 > $max_expand_len;
#$preexpand2 = 0 if $preexpand2 < 0;
	$preexpand2 = $max_expand_len if $preexpand2 < 0;
	$preexpand1 = $preexpand2 * $expand_factor if $preexpand1 > $preexpand2 * $expand_factor and $$cur_chain[2] != 3;
	$preexpand1 = $max_expand_len if $preexpand1 > $max_expand_len;

	$postexpand2 = $$next_chain[0][CHALO2] - $$cur_chain[0][CHAHI2];
	$postexpand2 = $postexpand1 * $expand_factor if $postexpand2 > $postexpand1 * $expand_factor and $$cur_chain[2] != 3;
	$postexpand2 = $max_expand_len if $postexpand2 > $max_expand_len;
#$postexpand2 = 0 if $postexpand2 < 0;
	$postexpand2 = $max_expand_len if $postexpand2 < 0;
	$postexpand1 = $postexpand2 * $expand_factor if $postexpand1 > $postexpand2 * $expand_factor and $$cur_chain[2] != 3;
	$postexpand1 = $max_expand_len if $postexpand1 > $max_expand_len;

	$$cur_chain[0][CHALO1E] = $$cur_chain[0][CHALO1] - $preexpand1;
	$$cur_chain[0][CHALO1E] = 1 if $$cur_chain[0][CHALO1E] < 1;
	$$cur_chain[0][CHAHI1E] = $$cur_chain[0][CHAHI1] + $postexpand1;
	$$cur_chain[0][CHAHI1E] = $$sizes1{$$cur_chain[0][SEQ1]} if $$cur_chain[0][CHAHI1E] > $$sizes1{$$cur_chain[0][SEQ1]};

	$$cur_chain[0][CHALO2E] = $$cur_chain[0][CHALO2] - $preexpand2;
	$$cur_chain[0][CHALO2E] = 1 if $$cur_chain[0][CHALO2E] < 1;
	$$cur_chain[0][CHAHI2E] = $$cur_chain[0][CHAHI2] + $postexpand2;
	$$cur_chain[0][CHAHI2E] = $$sizes2{$$cur_chain[0][SEQ2]} if $$cur_chain[0][CHAHI2E] > $$sizes2{$$cur_chain[0][SEQ2]};
	if ($cur_seq ne $$cur_chain[0][SEQ2]) { # Correct upper expansion
		$postexpand2 = $postexpand1 * $expand_factor;
		$postexpand2 = $max_expand_len if $postexpand2 > $max_expand_len;
		$postexpand2 = 0 if $postexpand2 < 0;
		$$cur_chain[0][CHAHI2E] = $$cur_chain[0][CHAHI2] + $postexpand2;
		$$cur_chain[0][CHAHI2E] = $$sizes2{$$cur_chain[0][SEQ2]} if $$cur_chain[0][CHAHI2E] > $$sizes2{$$cur_chain[0][SEQ2]};
	}

	print $WH $$cur_chain[2]." ".$$cur_chain[3]." ".
		$$cur_chain[0][START1]." ".$$cur_chain[0][END1]." ".$$cur_chain[0][START2]." ".$$cur_chain[0][END2]." ".
		$$cur_chain[0][SEQ1]." ".$$cur_chain[0][SEQ2]." ".$$cur_chain[0][ORIENT]." ".$$cur_chain[0][SCORE]." ".
		$$cur_chain[1][START1]." ".$$cur_chain[1][END1]." ".$$cur_chain[1][START2]." ".$$cur_chain[1][END2]." ".
		$$cur_chain[1][SEQ1]." ".$$cur_chain[1][SEQ2]." ".$$cur_chain[1][ORIENT]." ".$$cur_chain[1][SCORE]." ".
		$$cur_chain[0][CHALO1]." ".$$cur_chain[0][CHAHI1]." ".$$cur_chain[0][CHALO1E]." ".$$cur_chain[0][CHAHI1E]." ".
		$$cur_chain[0][CHALO2]." ".$$cur_chain[0][CHAHI2]." ".$$cur_chain[0][CHALO2E]." ".$$cur_chain[0][CHAHI2E];
	if ($print_chains) {
		my $i = 22;
		while (1) {
			print $WH " ".${$$cur_chain[4]}[$i]." ".${$$cur_chain[4]}[$i+1]." ".${$$cur_chain[4]}[$i+2]." ".${$$cur_chain[4]}[$i+3];
			last if @{$$cur_chain[4]} <= $i+4;
			$i+=4;
		}
	}
	print $WH "\n";
}


sub finalExpReg($$$$$) {
	my ($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq) = @_;
	my ($preexpand1, $postexpand1, $preexpand2, $postexpand2);
	if ($$cur_chain[2] == 1) { # M1: expand in seq1 on seq2 expands * factor only
		$preexpand1 = $$cur_chain[0][CHALO1] - $$cur_chain[0][CHALO1E];
		$preexpand2 = $preexpand1 * $expand_factor;
		$preexpand2 = $max_expand_len if $preexpand2 > $max_expand_len;
		$postexpand1 = $$cur_chain[0][CHAHI1E] - $$cur_chain[0][CHAHI1];
		$postexpand2 = $postexpand1 * $expand_factor;
		$postexpand2 = $max_expand_len if $postexpand2 > $max_expand_len;
		$$cur_chain[0][CHALO2E] = $$cur_chain[0][CHALO2] - $preexpand2;
		$$cur_chain[0][CHALO2E] = 1 if $$cur_chain[0][CHALO2E] < 1;
		$$cur_chain[0][CHAHI2E] = $$cur_chain[0][CHAHI2] + $postexpand2;
		$$cur_chain[0][CHAHI2E] = $$sizes2{$$cur_chain[0][SEQ2]} if $$cur_chain[0][CHAHI2E] > $$sizes2{$$cur_chain[0][SEQ2]};
	} elsif ($$cur_chain[2] == 3) { # M2: expand in seq2 on seq1 expands * factor only
		$preexpand2 = $$cur_chain[0][CHALO2] - $$cur_chain[0][CHALO2E];
		$preexpand1 = $preexpand2 * $expand_factor;
		$preexpand1 = $max_expand_len if $preexpand1 > $max_expand_len;
		$postexpand2 = $$cur_chain[0][CHAHI2E] - $$cur_chain[0][CHAHI2];
		$postexpand1 = $postexpand2 * $expand_factor;
		$postexpand1 = $max_expand_len if $postexpand1 > $max_expand_len;
		$$cur_chain[0][CHALO1E] = $$cur_chain[0][CHALO1] - $preexpand1;
		$$cur_chain[0][CHALO1E] = 1 if $$cur_chain[0][CHALO1E] < 1;
		$$cur_chain[0][CHAHI1E] = $$cur_chain[0][CHAHI1] + $postexpand1;
		$$cur_chain[0][CHAHI1E] = $$sizes1{$$cur_chain[0][SEQ1]} if $$cur_chain[0][CHAHI1E] > $$sizes1{$$cur_chain[0][SEQ1]};
	}

	print $WH $$cur_chain[0][SEQ1]." ".$$cur_chain[0][CHALO1E]." ".$$cur_chain[0][CHAHI1E]."   ".
		$$cur_chain[0][SEQ2]." ".$$cur_chain[0][CHALO2E]." ".$$cur_chain[0][CHAHI2E]." ".$$cur_chain[0][ORIENT];
	print $WH " (".($$cur_chain[2]==1?"M1, ":$$cur_chain[2]==2?"DM, ":"M2, ").$$cur_chain[3]." aligns)" unless $no_aligntotals;
	if ($print_chains) {
		my $i = 26;
		while (1) {
			print $WH " [".${$$cur_chain[4]}[$i]."-".${$$cur_chain[4]}[$i+1]."=".${$$cur_chain[4]}[$i+2]."-".${$$cur_chain[4]}[$i+3]."]";
			last if @{$$cur_chain[4]} <= $i+4;
			$i+=4;
		}
	}
	print $WH "\n";
}


sub finalExpand($$) {
	my ($RH, $WH) = @_;
	my ($first_align, $last_align, $type, $num_aligns,
		$cur_seq, $preexpand1, $postexpand1, $preexpand2, $postexpand2,
		$prev_chain, $cur_chain, $next_chain);
	my %stats;
	my (@line);

	while (<$RH>) {
		chomp; @line = split;

		$prev_chain = $cur_chain;
		$cur_chain = $next_chain;

		$first_align = []; $last_align = [];
		($type, $num_aligns, $$first_align[START1], $$first_align[END1], $$first_align[START2], $$first_align[END2],
		$$first_align[SEQ1], $$first_align[SEQ2],$$first_align[ORIENT], $$first_align[SCORE],
		$$last_align[START1], $$last_align[END1], $$last_align[START2], $$last_align[END2],
		$$last_align[SEQ1], $$last_align[SEQ2], $$last_align[ORIENT], $$last_align[SCORE],
		$$first_align[CHALO1], $$first_align[CHAHI1], $$first_align[CHALO1E], $$first_align[CHAHI1E],
		$$first_align[CHALO2], $$first_align[CHAHI2], $$first_align[CHALO2E], $$first_align[CHAHI2E]) = @line;

		if ($type == 1) {
			$$first_align[CHALO2] = ($$first_align[START2] < $$last_align[START2] ? $$first_align[START2] : $$last_align[START2]);
			$$first_align[CHAHI2] = ($$first_align[END2] > $$last_align[END2] ? $$first_align[END2] : $$last_align[END2]);
		} elsif ($type == 3) {
			$$first_align[CHALO1] = ($$first_align[START1] < $$last_align[START1] ? $$first_align[START1] : $$last_align[START1]);
			$$first_align[CHAHI1] = ($$first_align[END1] > $$last_align[END1] ? $$first_align[END1] : $$last_align[END1]);
		}

		my @saved_line = @line;
		$next_chain = [$first_align, $last_align, $type, $num_aligns, \@saved_line];

		next unless defined $cur_chain;

		finalExpReg($WH, $prev_chain, $cur_chain, $next_chain, $cur_seq);

		if ($debug or $print_stats) {
			if ($type == 1) {
				$$cur_chain[0][ORIENT] eq "+" ? $stats{"M1+"}++ : $stats{"M1-"}++;
				$stats{"M1min"} = $num_aligns if $stats{"M1min"} > $num_aligns or not defined $stats{"M1min"};
				$stats{"M1max"} = $num_aligns if $stats{"M1max"} < $num_aligns or not defined $stats{"M1max"};
				$stats{"M1mean"} += $num_aligns;
			} elsif ($type == 2) {
				$$cur_chain[0][ORIENT] eq "+" ? $stats{"DM+"}++ : $stats{"DM-"}++;
				$stats{"DMmin"} = $num_aligns if $stats{"DMmin"} > $num_aligns or not defined $stats{"DMmin"};
				$stats{"DMmax"} = $num_aligns if $stats{"DMmax"} < $num_aligns or not defined $stats{"DMmax"};
				$stats{"DMmean"} += $num_aligns;
			} else {
				$$cur_chain[0][ORIENT] eq "+" ? $stats{"M2+"}++ : $stats{"M2-"}++;
				$stats{"M2min"} = $num_aligns if $stats{"M2min"} > $num_aligns or not defined $stats{"M2min"};
				$stats{"M2max"} = $num_aligns if $stats{"M2max"} < $num_aligns or not defined $stats{"M2max"};
				$stats{"M2mean"} += $num_aligns;
			}
		}
		if ($cur_seq ne $$first_align[SEQ2]) {
			undef $cur_chain;
			$cur_seq = $$first_align[SEQ2];
		}
	}
	finalExpReg($WH, $cur_chain, $next_chain, undef, $cur_seq);

	if ($debug or $print_stats) {
		foreach my $i ("DM", "M1", "M2") {
			$stats{$i."mean"} /= ($stats{$i."+"} + $stats{$i."-"}) unless ($stats{$i."+"} + $stats{$i."-"} == 0);
			print $i.": ".($stats{$i."+"} + $stats{$i."-"})." chains (".$stats{$i."+"}."+, ".$stats{$i."-"}."-); ".
				"length min ".$stats{$i."min"}.", avg ".$stats{$i."mean"}.", max ".$stats{$i."max"}."\n";
		}
	}
}


# Called only in a "$0 worker" invocation
sub workerRun($$$$) {
	my ($tar_file, $score_file, $SLAGAN, $debug) = @_;
	my ($tmp_dir, $io_dir) = ($worker_tmp_dir, getcwd);
	local *FH;

	mkdir($tmp_dir) or die("$0 (worker): Could not create directory $tmp_dir: ".$!);

	copy($score_file, $tmp_dir);
	$score_file =~ /.*\/([^\/]+)$/;
	$score_file = $tmp_dir.$1;

	print("$0 (worker): Version ".$VERSION." started ".localtime()."\n") if $debug;
	print("$0 (worker): Jobfile=$tar_file, scorefile=$score_file, tmpdir=$tmp_dir, iodir=$io_dir, SLAGAN=$SLAGAN\n") if $debug;

	move($io_dir."/".$tar_file, $tmp_dir);
	my @files = `cd $tmp_dir; tar -xvf $tar_file` or warn("$0 (worker): Error extracting $tar_file");
	foreach my $file (@files) {
		chomp $file;
#print "$SLAGAN $tmp_dir$file $score_file > $tmp_dir$file.glocal-out 2> $tmp_dir$file.glocal-err\n";
		system("$SLAGAN $tmp_dir$file $score_file ".
			"> $tmp_dir$file.glocal-out ".
			"2> $tmp_dir$file.glocal-err");
	}

	$tar_file =~ /(.*)\.tar$/; $tar_file = $1;
	open(FH, "| cd $tmp_dir; xargs tar --append --file=$io_dir/$tar_file.results.tar");
	foreach my $file (glob("$tmp_dir/*glocal-out")) { $file =~ /\/([^\/]+)$/; print FH $1." "; }
	close FH;

	rmtree $tmp_dir;
	opendir(DIR, "."); if (my @x = grep(/core\./,readdir(DIR))) { warn("$0 (worker): WARNING: $SLAGAN crashed ".@x." times"); } closedir DIR;
	unlink(glob("core.*")) unless $nodelete;
}


# Interrupt handler
sub dequeueClustJobs($) {
	print "\n$0: Received SIG".$_[0].". Cleaning up... ";
	if ($clust_run_pid) {
		# send SIGQUIT to clust_run so it can dequeue cluster jobs
		kill "QUIT", $clust_run_pid;
	}
	unless ($debug or $nodelete) {
		print "Removing job files...";
		foreach my $i (1..$num_jobs-1) {
			unlink $tmp_dir."JOB".$i.".tar";
			unlink $tmp_dir."JOB".$i.".results.tar";
			unlink $tmp_dir."CLUSTER_JOB_MESSAGES.$i";
			unlink $tmp_dir."CLUSTER_JOB_ERRMSG.$i";
		}

		unlink "$tmp_dir$input_glob.chaos";
		unlink $tmp_dir."CLUSTER_JOB_PARAMS";
		rmtree($tmp_dir) if $ARGV[0] eq "worker";
	}
	print "\n";
	exit(1);
}


# Retrieve sequence length data from GPDB
sub get_all_seqs($$) {
	my ($dbh, $genome) = @_;
	my ($dset, $annot_db, $family, $check_chroms, %sizes, $chroms, @real_chroms,
	$ctgs);

	($dset, $annot_db, $family) = ($genome =~ /^\d+$/o) ?
	($genome + 0, ($dbh->get_data_set($genome))[4,14]) :
	($dbh->get_family_dset($genome))[0,4,14];
	print "$0: Genome $genome, dataset $dset, annotation db \"$annot_db\", family \"$family\"\n" if $debug;
	$annot_db and $check_chroms = 1;
	if ($check_chroms) {
		$chroms = $dbh->get_chroms(($dbh->get_data_set($dset))[2]);
		foreach my $chrom (@$chroms) {
			$$chrom[1] == 1 or next;
			my $name = "chr$$chrom[2]";
			my ($chr_id, $chr_type, $ctg_id, $size) =
			$dbh->find_seq($name, $dset, $annot_db);
			$chr_id and $sizes{$name} = $size;
		}
	}
	$ctgs = $dbh->selectcol("SELECT name FROM dset$dset\_contigs " .
							"WHERE name is not null and name != ? group by name", undef, "");
	foreach my $ctg (@$ctgs) {
		$sizes{$ctg} = $dbh->get_contig_size($dset, $ctg);
	}
	return \%sizes;
}


sub alnum {
	my ($i);
	my ($len1, $len2) = (length($a), length($b));
	for ($i = 0; ($i < $len1) && ($i < $len2); ++$i) {
		my $c1 = substr($a, $i, 1);
		my $c2 = substr($b, $i, 1);
		($c1 =~ /^\d/o) || ($c2 =~ /^\d/o) || ($c1 ne $c2) and last;
	}
	my $a_r = ($i < $len1) ? substr($a, $i) : "";
	my $b_r = ($i < $len2) ? substr($b, $i) : "";
	my ($a_n, $a_s) = ($a_r =~ /^(\d+)(.*)$/o);
	my ($b_n, $b_s) = ($b_r =~ /^(\d+)(.*)$/o);
	return (defined($a_n) && defined($b_n)) ?
	(($a_n <=> $b_n) || ($a_s cmp $b_s)) : ($a cmp $b);
}


sub isBLAT($) {
	my ($file) = @_;
	local *FH;
	open(FH, "< ".$file) or die("$0: Cannot open input file $file:  ".$!);
	my $line = <FH>;
	close FH;
	if ($line =~ /\A.+\s[\d]+\s[\d]+\;\s.+\s[\d]+\s[\d]+\;\sscore/) {
		return 0;
	} elsif ($line =~ /\A[^\s]+\s[\d]+\s[\d]+\s[^\s]+\s/) {
		return 1;
	} else {
		die("$0: Unknown input format in $file. Stopped");
	}
}


sub getMinSeqScore($) {
	my ($file) = @_;
	my $score; local *FH;
	open(FH, "< ".$file) or die("$0: Could not open SLAGAN scorefile $file: $!");
	while (<FH>) {
		# sample line: {+U+;+U-;-U+;-U-}{70000 0 0 0}
		/\{\+U\+\;.+\}.*\{(\d+)\s.+\}/;
		$score = $1 if $1;
	}
	close FH;
	die("$0: Could not determine min_seq_score from SLAGAN scorefile $file. Stopped") unless $score;
	print "$0: min_seq_score: $score\n" if $debug;
	return $score;
}


sub writeSizes($$) {
	my ($sizes, $outfile) = @_; local *FH;
	open(FH, "> ".$outfile) or die("$0: Could not open file $outfile for writing: ".$!);
	foreach my $key (sort alnum keys %$sizes1) {
		print FH $key."\t".$$sizes1{$key}."\n";
	}
	close FH;
}


# Borrowed from if.pm to enable standalone conditional module loading on earlier versions of Perl
sub useIf($$) {
	my $method = 'import';
	return unless shift; # CONDITION

	my $package = $_[0];
	(my $file = $package.".pm") =~ s!::!/!g;
	require $file;
	my $method_entry_point = $package->can($method);
	goto &$method_entry_point if $method_entry_point;
}


sub checkAlignCoords($) {
	my $cur_align = $_[0];
	if ($$cur_align[START1] > $$cur_align[END1]) { my $i = $$cur_align[START1]; $$cur_align[START1] = $$cur_align[END1]; $$cur_align[END1] = $i; }
	if ($$cur_align[START2] > $$cur_align[END2]) { my $i = $$cur_align[START2]; $$cur_align[START2] = $$cur_align[END2]; $$cur_align[END2] = $i; }

#	if ($$cur_align[OSTART1] > $$cur_align[OEND1]) { my $i = $$cur_align[OSTART1]; $$cur_align[OSTART1] = $$cur_align[OEND1]; $$cur_align[OEND1] = $i; }
#	if ($$cur_align[OSTART2] > $$cur_align[OEND2]) { my $i = $$cur_align[OSTART2]; $$cur_align[OSTART2] = $$cur_align[OEND2]; $$cur_align[OEND2] = $i; }
}


=head1 NAME

Supermap: Piecewise monotonic alignment map generator for shuffle-lagan

=head1 SYNOPSIS

supermap.pl (gen2=id | sizes2=filename) (gen1=id | sizes1=filename)
[-infile=<file>] [-outfile=<file>] [-bacteria] [-score=filename] [-f]
[file1 file2 ...]

=head1 EXAMPLES

supermap.pl -sizes1=human.sizes -sizes2=mouse.sizes hm.chr*.chaos

=head1 DESCRIPTION

Supermap is a whole-genome alignment map generator. It is an extension to the
Shuffle-LAGAN suite (Brudno et al., 2003). Supermap removes the asymmetry between
the query genomes by running multiple SLAGAN passes and combining them into a full
two-genome alignment.

To run Supermap without the Berkeley Genome Pipeline functionality, you will need
sequence length files for each of the genomes. Each file should contain one sequence
length entry per line, of the form "sequence_name sequence_length".

In the CHAOS output format (this program's input), negative orientation always means second pair of coords is inverted.
In this program's output, negative orientation does not invert coordinates (coordinate pairs are always ascending).

Run supermap.pl with no arguments to see a further description.

The terms "hit" and "anchor" usually refer to local alignments produced by CHAOS or another program.
The term "chain" refers to an extended union of a number of these local alignments.

=head1 DEPENDENCIES

Supermap depends on Utils.pm, SLAGAN, and a number of Unix utilities.

To use the Berkeley Genome Pipeline and cluster functionality, Supermap needs
GPutils.pm, GPDBI.pm, and clust_run.

=head1 LIMITATIONS

Supermap is designed to allow the manipulation of large datasets in a reasonable memory footprint.
To do this, it allows multiple files on input and keeps most of its intermediate data in small temporary files.
However, one current limitation is that the alignments for any sequence in either genome must fit into the largest
addressable file size (typically 2GB), and the output alignments must also fit in that size (the remainder will be truncated).

=head1 BUGS

=head1 TODO

TODO: bacteria description, examples, other input formats
TODO: installer routine
TODO: discuss input glob parameters
TODO: local multithreading
TODO: ignore escaped slashes when splitting dir/file (copy Alex)
TODO: check for ++ etc in SLAGAN out
TODO: .supermaprc file for score files, etc
TODO: hazelton.lbl.gov/bugzilla for supermap

=head1 AUTHOR

Andrey Kislyuk L<mailto:kislyuk@ocf.berkeley.edu>.

=cut
