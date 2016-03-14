#!/usr/bin/perl -w

use strict;

my $usage = "
./graph_scores.perl <plotname> <path1>/score_<suffix1> benchmark1 expconfig1 [datalabel1 ...] <path2>/score_<suffix2> benchmark2 expconfig2 [datalabel2 ...] ...

Generates a plot of scores from benchmark1, benchmark2, ... . The score file for each should either have a single score,
   or one or more scores corresponding to datalabel1..., etc.

If benchmark1 == benchmark2, etc., the plot will continue to group the corresponding expconfigs and datalabels.

Additionally while benchmark1 == benchmark2 and expconfig1 == expconfig2, datalabels are accumulated
   for (benchmark1, expconfig1).

If a benchmark has multiple expconfigs, i.e., benchmark1 == benchmark2, but expconfig1 != expconfig2,
     the datalabels accumulated for expconfig1 and expconfig2 must be identical and occur in the same order.

score_<plotname>_script.gnuplot and score_<plotname>_graphdata.gnuplot are produced.

Example:
  ./graph_scores.perl summaryplot \
     m/dolbycanyon_encode/score_sr_46.1 mpeg2enc '1X' canyon \
     m/dolbycity640x480_encode/score_sr_46.1 mpeg2enc '1X' city640x480 \
     m/dolbycanyon_encode/score_sr_46.1 mpeg2enc '2X' canyon \
     m/dolbycity640x480_encode/score_sr_46.1 mpeg2enc '2X' city640x480
";

die $usage if (@ARGV < 3);

print "# args -> @ARGV\n";

my $plotname = shift @ARGV;

my @bm = (); #array of distinct benchmark names

my @ec = (); #entries correspond to @bm. Each entry $ec[$i] is an array of experiment configurations (expconfigs) for $bm[$i]
my @ha = (); #entries correspond to @bm. Each entry $ha[$i] is an array with entries corresponding to expconfigs in $ec[$i]
             #$ha[$i][$j] is an hashmap for benchmark $bm[$i] and expconfig $ec[$i]$[j]:
             #     {"dl", "sc", "as", "ar"} -> {datalabel,
             #                                  score for datalabel,
             #                                  average metric value for corresponding score
             #                                  average closest reference metric value for corresponding score}
#my @sc = (); #nested entries correspond to @bm and @dl. Each leaf entry is a score for corresponding datalabel.
#my @as = (); #average metric value for corresponding entry of @sc
#my @ar = (); #average closest reference metric value for corresponding entry of @sc

# states
# - find_path (initial state)
# - find_benchname
# - find_expconfig
# - read_label (final state)

my $state = "find_path";

my $active_fh;
my $active_file = "";
my $active_bm = "";
my $active_ec = "";
my $prev_bm = "";
my $prev_ec = "";
foreach my $a (@ARGV) {
  if($state eq "find_path") {
    if($a =~ /(.*score_.*)/) {

      $active_file = $1;
      open $active_fh, "<", $active_file || die "Could not open file $active_file for reading: $!";

      $active_bm = "";

      $state = "find_benchname";
    }
    else {
      die "Expected a path to benchmark score file";
    }
  }
  elsif($state eq "find_benchname") {
    $active_bm = $a;

    if($active_bm ne $prev_bm) {
      push @bm, $active_bm;
      push @ec, [];
      push @ha, [];
    }

    $state = "find_expconfig";
  }
  elsif($state eq "find_expconfig") {
    $active_ec = $a;

    if($active_ec ne $prev_ec) {
      push @{$ec[$#ec]}, $active_ec;
      push @{$ha[$#ha]}, [];
    }

    $state = "read_label";
  }
  elsif($state eq "read_label") {
    if($a =~ /(.*score_.*)/) { # actually the next benchmark path
      my $last_index = $#{$ha[$#ha]};
      if(@{$ha[$#ha][$last_index]} == 0) { # no labels provided, read one score for entire benchmark
        my ($score, $avg_metric_s, $avg_metric_r) = read_next_score();
        save_score("", $score, $avg_metric_s, $avg_metric_r);
      }

      $prev_bm = $active_bm;
      $prev_ec = $active_ec;

      close $active_fh || die "Could not close file $active_file after reading: $!";

      $active_file = $1;
      open $active_fh, "<", $active_file || die "Could not open file $active_file for reading: $!";

      $active_bm = "";
      $active_ec = "";
      $state = "find_benchname";
    }
    else { # a data label
      my ($score, $avg_metric_s, $avg_metric_r) = read_next_score();
      save_score($a, $score, $avg_metric_s, $avg_metric_r);
    }
  }
  else {
    die "Invalid state: $state\n";
  }

  print "arg=$a file=$active_file name=$active_bm\n";
}
die "Incomplete score set" if($state ne "read_label");
my $last_index = $#{$ha[$#ha]};
if(@{$ha[$#ha][$last_index]} == 0) { # no labels provided, read one score for entire benchmark
  my ($score, $avg_metric_s, $avg_metric_r) = read_next_score();
  save_score("", $score, $avg_metric_s, $avg_metric_r);
}

sub read_next_score {
  my $fp = "(-?\\d+?(.\\d+)?(e[+-]\\d+)?)";
  my $pat = "score .+: ${fp}, ${fp}, ${fp},";
  while(1) {
    die "EOF when score for label was expected in $active_file" if(eof($active_fh));

    my $line = <$active_fh>;
    if($line =~ /$pat/) {
      my $score = $1;
      my $avg_metric_s = $4;
      my $avg_metric_r = $7;
      return ($score, $avg_metric_s, $avg_metric_r);
    }
  }
}

sub save_score {
  my ($label, $score, $avg_metric_s, $avg_metric_r) = @_;

  my %dl_entry = ();
  $dl_entry{"dl"} = $label;
  $dl_entry{"sc"} = $score;
  $dl_entry{"as"} = $avg_metric_s;
  $dl_entry{"ar"} = $avg_metric_r;

  my $last_index = $#{$ha[$#ha]};
  push @{$ha[$#ha][$last_index]}, \%dl_entry;
}

print "----- Benchmarks -----\n";
for(my $i=0; $i<@bm; $i++) {
  print "$bm[$i]: ";
  for(my $j=0; $j<@{$ha[$i]}; $j++) {
    for(my $k=0; $k<@{$ha[$i][$j]}; $k++) {
      my %dl_entry = (%{$ha[$i][$j][$k]});
      my ($dl, $sc, $as, $ar) = ($dl_entry{"dl"}, $dl_entry{"sc"}, $dl_entry{"as"}, $dl_entry{"ar"});
      print "$dl - ($sc, $as, $ar) ";
    }
  }
  print "\n";
}


for(my $i=0; $i<@bm; $i++) {
  for(my $j=0; $j<@{$ha[$i]}; $j++) {
    for(my $k=0; $k<@{$ha[$i][$j]}; $k++) {
      my %first_dl_entry = (%{$ha[$i][0][$k]});
      my ($fdl, $fsc, $fas, $far) = ($first_dl_entry{"dl"}, $first_dl_entry{"sc"}, $first_dl_entry{"as"}, $first_dl_entry{"ar"});

      my %dl_entry = (%{$ha[$i][$j][$k]});
      my ($dl, $sc, $as, $ar) = ($dl_entry{"dl"}, $dl_entry{"sc"}, $dl_entry{"as"}, $dl_entry{"ar"});

      #verify that each expconfig has the same sequence of datalabels.
      die "Datalabel sequence must be identical for expconfigs " . $ec[$i][0] . " and " . $ec[$i][$j] . "with $dl vs $fdl (j=$j)" if($dl ne $fdl);
    }
  }
}

#### Generate score graphdata files
for(my $i=0; $i<@bm; $i++) {
  my $score_plotfile = "score_${plotname}_graphdata_$bm[$i].gnuplot";
  open PF, ">", $score_plotfile or die "Could not open score_plotfile=$score_plotfile: $!\n";

  my $title = "config";
  for(my $k=0; $k<@{$ha[$i][0]}; $k++) {
    my $dl = $ha[$i][0][$k]{"dl"};
    $title .= " \"$dl\"";
  }
  
  print PF "$title\n";

  for(my $j=0; $j<@{$ha[$i]}; $j++) {
    print PF "\"$ec[$i][$j]\" ";
    for(my $k=0; $k<@{$ha[$i][$j]}; $k++) {
      my %dl_entry = (%{$ha[$i][$j][$k]});
      my $score = $dl_entry{"sc"};
      print PF "$score ";
    }
    print PF "\n";
  }
  close PF or die "Could not close score_plotfile=$score_plotfile: $!\n";
}

### Generate score script file
my $score_scriptfile = "score_${plotname}_script.gnuplot";
open SF, ">", $score_scriptfile or die "Could not open score_scriptfile=$score_scriptfile: $!\n";
my $command = "
#set terminal pngcairo size 1280,640 enhanced font 'Verdana,16'
#set output '${plotname}.png'
set title '${plotname}'
#set yrange [-0.2:1.0]
set ylabel 'Improvement to SR'
set style data histograms
set style histogram clustered gap 1
set key
unset xtics
set style fill solid #pattern
set boxwidth 0.95
set xtics out nomirror
";

$command .= "plot ";
for(my $i=0; $i<@bm; $i++) {
  $command .= "     " if($i > 0);
  my $score_plotfile = "score_${plotname}_graphdata_$bm[$i].gnuplot";
  $command .= "newhistogram '$bm[$i]' font 'Verdana,18', ";
  for(my $k=0; $k<@{$ha[$i][0]}; $k++) {
    my $label = $ha[$i][0][$k]{"dl"};
    $label = "X" if($label eq ""); #FIXME
    my $scpf = "''";
    $scpf = "'$score_plotfile'" if($k == 0);
    my $delim = ", ";
    $delim .= "\\\n" if($k == $#{$ha[$i][0]});
    my $dl = $ha[$i][0][$k]{"dl"};
    $command .= "$scpf using '$dl':xtic(1)$delim";
  }
}
$command .= "     0 lt black notitle\n";

print SF $command;

exit;

my @dl;
my @sc;
my @as;
my @ar;

### Generate comparison graphdata files
for(my $i=0; $i<@bm; $i++) {
  my $comp_plotfile = "comp_${plotname}_graphdata_$bm[$i].gnuplot";
  open PF, ">", $comp_plotfile or die "Could not open comp_plotfile=$comp_plotfile: $!\n";

  print PF "DS controller bestfixed\n";

  for(my $k=0; $k<@{$ha[$i][0]}; $k++) {
    my $label = $ha[$i][0][$k]{"dl"};
    $label = "X" if($label eq ""); #FIXME
    my $comp_s = $ha[$i][0][$k]{"as"};
    my $comp_r = $ha[$i][0][$k]{"ar"};
    print PF "$label $comp_s $comp_r\n";
  }
  close PF or die "Could not close comp_plotfile=$comp_plotfile: $!\n";
}

my @dn = ();

### Generate comparison script file
my $comp_scriptfile = "comp_${plotname}_script.gnuplot";
open SF, ">", $comp_scriptfile or die "Could not open comp_scriptfile=$comp_scriptfile: $!\n";
$command = "
#set terminal pngcairo size 1280,640 enhanced font 'Verdana,16'
#set output '${plotname}.png'
set title '${plotname}'
set yrange [0:1.0]
set ylabel 'SR'
set style data histograms
set style histogram clustered gap 1
#set nokey
unset xtics
set style fill solid noborder
set boxwidth 0.95
set xtics out nomirror rotate by -45
";

$command .= "plot ";
for(my $i=0; $i<@bm; $i++) {
  $command .= ", " if($i > 0);
  my $comp_plotfile = "comp_${plotname}_graphdata_$bm[$i].gnuplot";
  my ($controller_title, $bestfixed_title) = ("title ''", "title ''");
  ($controller_title, $bestfixed_title) = ("title 'controller'", "title 'bestfixed'") if($i == 0);
  $command .= "newhistogram '$dn[$i]' font 'Verdana,18', '$comp_plotfile' using 'controller':xtic(1) ls 1 $controller_title, '' using 'bestfixed' ls 2 $bestfixed_title";
}
$command .= "\n";

print SF $command;


