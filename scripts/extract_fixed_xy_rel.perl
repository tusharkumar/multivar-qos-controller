#!/usr/bin/perl -w

use strict;

my $usage =
"extract_fixed_xy_rel.pl <log-file-suffix> [index=0]
   <log-file-suffix> identifies a set of log files by excluding the 'log_fixed<i>_' prefix.
   That is, pattern 'log_fixed*_<log-file-suffix>' is used to name the log-files to scan
   and identify the fixed settings '<i>'.

   If the log-files have multiple outputs, a specific output can be chosen via index (=0 default).

   Two output files will be produced showing a plot between the identified fixed values and first output dim.
     - 'fixed_xy_script_<log-file-suffix>'
     - 'fixed_xy_graphdata_<log-file-suffix>'
";

die $usage if(@ARGV < 1 || @ARGV > 2);

my $suffix = $ARGV[0];
print "suffix=${suffix}\n";

my $index = 0;
if(@ARGV == 2) {
  $index = $ARGV[1] + 0;
  print "Using specified index=${index}\n";
}

my @fixedlogfiles = <log_fixed*_${suffix}>;

die "No log files match suffix=${suffix}\n" if(@fixedlogfiles == 0);

my $scriptfilename = "fixed_xy_script_${suffix}";
open SFP, "> $scriptfilename" or die "extract_fixed_xy_rel.perl: ERROR: Could not open script-file for writing: $!\n";

my $graphdatafilename = "fixed_xy_graphdata_${suffix}";
open GFP, "> $graphdatafilename" or die "extract_fixed_xy_rel.perl: ERROR: Could not open graphdata-file for writing: $!\n";

foreach my $logfilename (@fixedlogfiles) {
  my $fixedval;
  if($logfilename =~ /log_fixed(.+)_${suffix}/) {
    $fixedval = $1;
  }
  else {
    die "Log file $logfilename does not have suffix=${suffix}";
  }
  print "fixedval=$fixedval\n";


  open LFP, $logfilename or die "Could not open log-file: $logfilename: $!";
  my @yseq;
  while(<LFP>) {
    if(/READ frame-normalized-observations=\[(.*)\]/) {
      my @read_data = split(/,/, $1);
      die "Invalid index=${index}" if($index < 0 || $index >= @read_data);
      my $y = $read_data[$index];
      push @yseq, $y;
    }
  }

  my @sorted_yseq = sort {$a <=> $b} @yseq;
  my $len = @sorted_yseq;
  die "Data sequence too short: len=${len}" if($len < 20);

  my $sum = 0;
  foreach my $y (@sorted_yseq) {
    $sum += $y;
  }
  my $mean = $sum / $len;

  print GFP "$fixedval $sorted_yseq[$len * 0.10], $sorted_yseq[$len * 0.25], $mean, $sorted_yseq[$len * 0.75], $sorted_yseq[$len * 0.90]\n";
}

my $plot_command;

$plot_command = "plot '$graphdatafilename' using 1:3:2:6:5 title 'spread' with candlesticks, '$graphdatafilename' using 1:4:4:4:4 title 'mean' with candlesticks;";
print SFP "$plot_command\n\n";

