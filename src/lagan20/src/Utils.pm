#!/usr/bin/env perl

package Utils;
require 5.000;

use strict;
use Exporter;
use Cwd;
use IO::File;
use POSIX qw(setsid);
use Sys::Syslog qw(:DEFAULT setlogsock);

sub Trim( @ );
sub Lock_File( $ ; $ $ $ );
sub Unlock_File( $ );
sub Write_Log( $ $ ; $ $ );
sub Parse_Filename( $ );
sub Get_Abs_Path( $ );
sub Expand_Path( $ );
sub Get_Random_Key( ; $ );
sub Hex2Ascii( $ );
sub Ascii2Hex( $ );
sub Get_Config_Record( $ $ );
sub Round( $ );
sub Set_Log( $ $ );
sub Log( $ $ );
sub Min( $ $ );
sub Max( $ $ );
sub Reg_Diff( $ $ ; $ $ $ $ $ );
sub Reg_Rem_Overlap( $ ; $ $ $ );
sub Reg_Sort( $ ; $ $ $ );
sub Reg_Intersect( $ $ ; $ $ $ $ $ );
sub Reg_Merge( $ ; $ $ $ );

use vars qw(@ISA @EXPORT $VERSION $JOB $Error $Syslog $Facility $Msg_Prefix);

@ISA = qw(Exporter);
@EXPORT = qw(Trim Lock_File Unlock_File Write_Log Parse_Filename 
  Get_Abs_Path Expand_Path Hex2Ascii Ascii2Hex Get_Config_Record 
  Get_Random_Key Round Set_Log Log Min Max Reg_Diff Reg_Rem_Overlap
  Reg_Sort Reg_Intersect Reg_Merge redirect_err2log openlogs safe_glob
  daemon wr_log wr_err start_watcher confirm $JOB);

my $Id = '$Id: Utils.pm,v 1.21 2005/01/07 23:08:59 poliakov Exp $';
($VERSION) = ($Id =~ /,v\s+(\d+\S+)/o);
$JOB = '^(\S+)\@(\S+?)_(\d{4})(?:_(.+)|)$';

$Error = 0;
$Syslog = 0;
$Facility = "user";
$Msg_Prefix = undef;

my $E_FORK = "cannot fork";
my @LOG_FILE = ();
my %Locks = ();

sub Trim( @ ) {
  for (my $i = 0; $i <= $#_; ++$i) {
    $_[$i] =~ s/^\s+//;
    $_[$i] =~ s/\s+$//
  }
}

sub Lock_File( $ ; $ $ $ ) {
  my ($file, $retry, $timeout, $max_mtime) = @_;
  my ($lock_fh, $start_time, $mtime);

  if (!$file || ($file =~ /\/$/o)) {
    $Error = "Invalid filename";
    return 0;
  }
  $file = Get_Abs_Path("$file.lock");
  if (exists($Locks{$file})) { $Error = "Already locked"; return 1; }
  if (!-w (Parse_Filename($file))[0]) {
    $Error = "Permission denied";
    return 0;
  }
  if (!defined($retry)) { $retry = 1; }
  if (!defined($timeout)) { $timeout = 1200; }
  if (!defined($max_mtime)) {
    $max_mtime = ($timeout > 0) ? int($timeout / 2) : 0;
  }
  $start_time = time();
  LOCK: {
    if (!($lock_fh = IO::File->new($file, O_RDWR|O_CREAT|O_EXCL))) {
      if (!$retry || (($timeout > 0) && ((time() - $start_time) > $timeout))) {
        $Error = "Locked by someone else";
	return 0;
      }
      if ($max_mtime > 0) {
        $mtime = (stat($file))[9];
        if ($mtime && ((time() - $mtime) > $max_mtime)) { unlink($file); }
      }
      redo LOCK;
    }
  }
  $lock_fh->close();
  $Locks{$file} = 1;
  return 1;
}

sub Unlock_File( $ ) {
  my ($file) = @_;

  if (!$file) { $Error = "Invalid filename"; return 0; }
  $file = Get_Abs_Path("$file.lock");
  if (!exists($Locks{$file})) { $Error = "Not locked"; return 0; }
  if (!unlink($file)) { $Error = "Cannot unlock"; return 0; }
  delete($Locks{$file});
  return 1;
}

{
  my $Uname;
  foreach my $dir ('/bin', '/sbin', '/usr/bin', '/usr/sbin') {
    -x "$dir/uname" and $Uname = "$dir/uname", last;
  }
  my $Host = $Uname ? `$Uname -n` : 'localhost';
  chomp($Host);
  ($Host) = ($Host =~ /^([^\.]+)(\..*)?$/);

sub Write_Log( $ $ ; $ $ ) {
  no strict "refs";
  my ($log_file, $msg, $name, $pid) = @_;
  my $error = 0;
  my $date;
  local *LOG;

  if (!defined($log_file) || !defined($msg)) { return 0; }
  if (*{$log_file}{IO}) {
    *LOG = *{$log_file}{IO};
  } elsif ($log_file eq '/dev/null') {
    return 1;
  } else {
    if (!Lock_File($log_file)) { return 0; }
    if (!open(LOG, ">> $log_file")) { $error = 1; }
  }
  if (!$error) {
    chomp($msg);
    $date = localtime(time());
    if (!$name) { $name = $0; }
    if (!$pid) { $pid = $$; }
    if (!print LOG "$date $Host $name\[$pid\]:  $msg\n") { $error = 1; }
    if (!*{$log_file}{IO}) { close(LOG); }
  }
  if ($error && $!) { $Error = "$!"; }
  if (!*{$log_file}{IO}) { Unlock_File($log_file); }
  return !$error;
}}

sub Parse_Filename( $ ) {
  my ($name) = @_;
  my ($last_slash_pos, $dir, $file);
  
  if (!defined($name)) { return (); }
  $last_slash_pos = rindex($name, "/");
  if ($last_slash_pos >= 0) {
    $dir = substr($name, 0, $last_slash_pos + 1);
    $file = substr($name, $last_slash_pos + 1);
  } else {
    $dir = "";
    $file = $name;
  }
  return ($dir, $file);
}

sub Expand_Path( $ ) {
  my ($path) = @_;
  my $home_dir;
  
  $path && ($path =~ /^~/o) or return $path;
  $path =~ /^~([^\/]*)(.*)$/o;
  $home_dir = $1 ? (getpwnam($1))[7] :
    ($ENV{"HOME"} || $ENV{"LOGDIR"} || (getpwuid($>))[7]);
  defined($home_dir) and $path = "$home_dir$2";
  return $path;
}

sub Get_Abs_Path( $ ) {
  my ($path) = @_;

  defined($path) or return $path;
  $path = Expand_Path($path);
  $path =~ /^\//o or $path = getcwd() . "/$path";
  $path =~ s(/{2,})(/)g;
  
# get rid of "/./"

  while ($path =~ /^(.*?)\/\.(?:|\/(.*))$/o) {  
    $path = "$1/" . ($2 ? $2 : "");
  }
  
# get rid of "/../"

  while ($path =~ /^(((?:.*?\/)*?)[^\/]+){0,1}?\/\.\.(?:|\/(.*))$/o) {
    $path = ($1 ? $2 : "/") . ($3 ? $3 : "");
  }
  return $path;
}

{
  my @Chars = ("A" .. "Z", "a" .. "z", 0 .. 9);
  srand();

sub Get_Random_Key( ; $ ) {
  my ($len) = @_;
  
  if (!defined($len) || ($len !~ /^\d+$/o) || ($len < 2) || ($len > 1024)) {
    $len = 8;
  }
  return join("", @Chars[map {rand @Chars } (1 .. 8)]);
}}

sub Hex2Ascii( $ ) {
  my ($str) = @_;
  
  if ($str) { $str =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg; }
  return $str;
}

{
  my $a2h = {
    "\t" => "%29",
    "+" => "%2B",
    "," => "%2C",
    "." => "%2E",
    ";" => "%3B",
    "/" => "%2F",
    "?" => "%3F",
    ":" => "%3A",
    "@" => "%40",
    "=" => "%3D",
    "&" => "%26",
    " " => "%20",
    "<" => "%3C",
    ">" => "%3E",
    "\"" => "%22",
    "%" => "%25",
    "#" => "%23",
    "[" => "%5B",
    "]" => "%5D",
    "{" => "%7B",
    "}" => "%7D",
    "|" => "%7C",
    "\\" => "%5C",
    "^" => "%5E",
    "~" => "%7E",
    "`" => "%60"};	

sub Ascii2Hex( $ ) {
  my ($str) = @_;
  my $new_str = "";

  if (!$str) { return $str; }
  foreach my $char (split(//, $str)) {
    if (exists($a2h->{$char})) { $char = $a2h->{$char}; }
    $new_str .= $char;
  }
  return $new_str;
}}

sub Get_Config_Record( $ $ ) {
  my ($conf_file, $rec) = @_;
  my ($db, $field, $value);
  my @result = ();

  if (!($db = Registry->New($conf_file, "r", 1))) {
    $Error = "$Registry::Error", return ();
  }
  if (!$db->Record_Exists($rec)) {
    $Error = qq("$rec" record not found);
    return ();
  }
  foreach my $field (qw(dir users log)) {
    if (!($value = Expand_Path($db->Get_Val($rec, $field)))) {
      if ($field eq "log") {
        $value = "";
      } else {
        $Error = qq("$field" field of "$rec" record is missing), return ();
      }
    } elsif ($value !~ /^\//o) {
      $Error = qq("$field" field of "$rec" record should be absolute path);
      return ();
    }
    push(@result, $value);
  }
  foreach my $field (qw(max_down grace_period)) {
    if (!($value = $db->Get_Val($rec, $field)) ||
        ($value !~ /^\d+$/o)) {
      $value = 0;
    }
    push(@result, $value);
  }
  return @result;
}

sub Round( $ ) {
  my ($num) = @_;
  
  return int($num + 0.5);
}

sub Log( $ $ ) {
  my ($log_num, $msg) = @_;

  (defined($log_num) && ($log_num >= 0) && $LOG_FILE[$log_num]) and
    Write_Log($LOG_FILE[$log_num], $msg);
}

sub Set_Log( $ $ ) {
  my ($log_num, $file) = @_;
  
  (defined($log_num) && ($log_num >= 0) && $file) and
    $LOG_FILE[$log_num] = $file;
}

sub Min( $ $ ) {
  my ($i, $j) = @_;
  
  return ($i < $j) ? $i : $j;
}

sub Max( $ $ ) {
  my ($i, $j) = @_;
  
  return ($i > $j) ? $i : $j;
}

sub Reg_Diff( $ $ ; $ $ $ $ $ ) {
  my ($regs1, $regs2, $strict, $s1, $e1, $s2, $e2) = @_;
  my (@new_regs, $start, $end, $new_reg);
  
  $regs1 && $regs2 or return $regs1;
  $s1 ||= 0;
  defined($e1) or $e1 = 1;
  $s2 ||= 0;
  defined($e2) or $e2 = 1;
  for (my $i = 0; $i < @$regs1; ++$i) {
    $start = $$regs1[$i][$s1];
    $end = $$regs1[$i][$e1];
    for (my $j = 0; $j < @$regs2; ++$j) {
      $$regs2[$j][$s2] > $end and last;
      $$regs2[$j][$e2] < $start and next;
      if (($$regs2[$j][$s2] <= $start) && ($$regs2[$j][$e2] >= $end)) {
        undef($start), last;
      }
      if (($$regs2[$j][$s2] > $start) && ($$regs2[$j][$e2] >= $end)) {
        $end = $$regs2[$j][$s2] - 1, last;
      }
      if (($$regs2[$j][$s2] <= $start) && ($$regs2[$j][$e2] < $end)) {
        $start = $$regs2[$j][$e2] + 1, next;
      }
      ($start < ($$regs2[$j][$s2] - 1)) || !$strict and
        $new_reg = [@{$$regs1[$i]}],
        $$new_reg[$s1] = $start,
        $$new_reg[$e1] = $$regs2[$j][$s2] - 1,
        push(@new_regs, $new_reg);
      $start = $$regs2[$j][$e2] + 1;
    }
    !defined($start) || ($start > $end) and next;
    ($start < $end) || !$strict and
      $new_reg = [@{$$regs1[$i]}],
      $$new_reg[$s1] = $start,
      $$new_reg[$e1] = $end,
      push(@new_regs, $new_reg);
  }
  return \@new_regs;
}

sub Reg_Rem_Overlap( $ ; $ $ $ ) {
  my ($regs, $strict, $s, $e) = @_;
  my (@new_regs);
  
  $regs or return $regs;
  $s ||= 0;
  defined($e) or $e = 1;
  for (my $i = 0; $i < @$regs; ++$i) { push(@new_regs, [@{$$regs[$i]}]); }
  for (my $i = 0; $i < @new_regs; ++$i) {
    if (($i < $#new_regs) && ($new_regs[$i + 1][$s] <= $new_regs[$i][$e])) {
      $new_regs[$i + 1][$e] <= $new_regs[$i][$e] and
        splice(@new_regs, $i + 1, 1),
        --$i, next;
      $new_regs[$i + 1][$s] = $new_regs[$i][$e] + 1;
    }
    ($new_regs[$i][$s] < $new_regs[$i][$e]) || !$strict and next;
    splice(@new_regs, $i, 1);
    --$i;
  }
  return \@new_regs;
}

sub Reg_Sort( $ ; $ $ $ ) {
  my ($regs, $rev, $s, $e) = @_;
  my (@new_regs);
  
  $regs or return $regs;
  $s ||= 0;
  defined($e) or $e = 1;
  if ($rev) {
    @new_regs = sort { ($$b[$s] <=> $$a[$s]) || ($$b[$e] <=> $$a[$e]) } @$regs;
  } else {
    @new_regs = sort { ($$a[$s] <=> $$b[$s]) || ($$a[$e] <=> $$b[$e]) } @$regs;
  }
  return \@new_regs;
}

sub Reg_Intersect( $ $ ; $ $ $ $ $ ) {
  my ($regs1, $regs2, $strict, $s1, $e1, $s2, $e2) = @_;
  
  $regs1 && $regs2 or return undef;
  $s1 ||= 0;
  defined($e1) or $e1 = 1;
  $s2 ||= 0;
  defined($e2) or $e2 = 1;
  return Reg_Diff($regs1, Reg_Diff($regs1, $regs2, $strict, $s1, $e1,
    $s2, $e2), $strict, $s1, $e1, $s1, $e1);
}

sub Reg_Merge( $ ; $ $ $ ) {
  my ($regs, $strict, $s, $e) = @_;
  my (@new_regs);
  
  $regs or return $regs;
  $s ||= 0;
  defined($e) or $e = 1;
  for (my $i = 0; $i < @$regs; ++$i) { push(@new_regs, [@{$$regs[$i]}]); }
  for (my $i = 0; $i < @new_regs; ++$i) {
    ($i < $#new_regs) &&
        ($new_regs[$i + 1][$s] == ($new_regs[$i][$e] + 1)) and
      $new_regs[$i][$e] = $new_regs[$i + 1][$e],
      splice(@new_regs, $i + 1, 1),
      --$i, next;
  }
  for (my $i = 0; $i < @new_regs; ++$i) {
    ($new_regs[$i][$s] < $new_regs[$i][$e]) || !$strict and next;
    splice(@new_regs, $i, 1);
    --$i;
  }
  return \@new_regs;
}

sub safe_glob {
  my ($regexp, $dir) = @_;
  my (@files);
  local (*DIR);
  
  $dir ||= ".";
  $regexp ||= ".*";
  opendir(DIR, $dir) or return;
  @files = grep { /$regexp/ } readdir(DIR);
  closedir(DIR);
  return wantarray() ? @files : scalar(@files);
}

sub redirect_err2log {
  my ($facility) = @_;
  
  $Facility = $facility;
  stderr2log();
}

sub stderr2log {
  my ($oldfh);
  
  open(STDERR, "> /dev/null");
  open(STDERR, "| logger -p $Facility.err -t '$0\[$$\]'");
  $oldfh = select(STDERR); $| = 1; select($oldfh);  
}

sub openlogs {
  my ($facility) = @_;
  
  $facility and $Facility = $facility;
  stderr2log();
  setlogsock("unix");
  openlog($0, "pid", $Facility);
  $Syslog = 1;
}

sub daemon {
  my ($facility) = @_;
  my ($pid);

  if ($pid = fork()) {
    exit(0);
  } elsif (!defined($pid)) {
    wr_err("$E_FORK: $!");
    die;
  } else {
    setsid();
    close(STDIN);
    close(STDOUT);
    open(STDOUT, "> /dev/null");
    openlogs($facility);
  }
}

sub start_watcher {
  my ($watcher, $facility, @params) = @_;
  my ($pid, $parent);

  $parent = $$;
  if ($pid = fork()) {
    return;
  } elsif (!defined($pid)) {
    wr_err("$E_FORK: $!");
    die;
  } else {
    setsid();
    close(STDIN);
    close(STDOUT);
    open(STDOUT, "> /dev/null");
    $0 .= "_watcher";
    openlogs($facility);
    &$watcher($parent, @params);
  }
}

sub wr_log {
  my $msg = shift;
  
  chomp($msg);
  $msg = ( $Msg_Prefix ? &$Msg_Prefix : "") . $msg;
  if ($Syslog) {
    syslog("info", "%s", $msg);
  } else {
    print "$msg\n";
  }
}

sub wr_err {
  my $msg = shift;
  
  chomp($msg);
  print STDERR (( $Msg_Prefix ? &$Msg_Prefix : ""), "$msg\n");
  return 1;
}

sub confirm {
  my ($msg) = @_;
  my ($ans);
  
  print $msg;
  $ans = <STDIN>;
  chomp($ans);
  return ($ans =~ /^(y|yes)$/io) ? 1 : 0;
}

END {
  foreach my $lock (keys(%Locks)) { unlink($lock); }
}

1;
