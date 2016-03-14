#!/usr/bin/perl -w

use strict;
use Cwd;
use File::Basename;

my $usage = "
Usage: score.perl [-q] <path1/logname1> [<tag> <tag> ...] [<path2/logname2> [<tag> <tag> ...]] ...
where:
  path1 is the directory containing the graphdata_{sr|mseq}_dsall-<logname1>_* file (assumed current dir if only logname1 specified).
  tag = {+|-}{srt | fixed}[_identifiers], to limit scoring to specific datasets in the graphdata file, if specified.
     A + indicates that a score should be computed for the specified dataset.
     A - indicates that the specified dataset should only be used in computation of scores for other datasets.

  path2/logname2, path3/logname3 --- additional lognames can be specified, including the particular dataset tags to limit to in each log.

  -q suppresses detailed output, only the computed scores are printed one per line.

Examples:
  ./score 43.3 +srt -fixed            # score all srt runs against all the fixed case runs
  ./score 43.3 +srt_43.3 -fixed_40.1  # score a specific srt run against a specific run of the fixed cases
  ./score mpeg2enc/dolbycanyon_encode/43.3 +srt_43.3 -srt_43.3 \
          2parm_mpeg2enc/dolbycanyon_encode/43.3 +srt_43.3 -srt_43.3 \
          4parm_mpeg2enc/dolbycanyon_encode/43.3 +srt_43.3 -srt_43.3 
      # score the specific runs of srt in the 1, 2 and 4 parameter runs of the mpeg2enc benchmark against each other
";

die $usage if (@ARGV < 1);

print "# args -> @ARGV\n";

my $VERBOSE = 1;
if($ARGV[0] eq "-q") {
  $VERBOSE = 0;
  shift @ARGV;
}

my @logs = ();
my @tags = ();

my $max_index = -1;
foreach my $p (@ARGV) {
  if($p =~ /^\+/ || $p =~ /^-/) {
    die "logname must be specified before tag." if($max_index == -1);
    push(@{$tags[$max_index]}, $p);
  }
  else { # start of a new log
    $max_index++;
    push(@logs, $p);
    push(@tags, []);
  }
}

for(my $i=0; $i <= $max_index; $i++) {
  print "log = ", $logs[$i], " tags = " if($VERBOSE == 1);
  foreach my $t (@{$tags[$i]}) {
    print "$t " if($VERBOSE == 1);
  }
  print "\n" if($VERBOSE == 1);
} 

# Read graphdata files
my @tablevalues = ();
for(my $i=0; $i <= $max_index; $i++) {
  my $logname = $logs[$i];
  my $dn = dirname($logname);
  my $ln = basename($logname);
  if($dn eq "") {
    $dn = ".";
  }
  my $pattern = "${dn}/graphdata_sr_dsall-${ln}_*";
  my @candidates = glob $pattern;
  die "Could not find candidate graphdata file for $logname using pattern $pattern" if(@candidates == 0);
  my $graphdatafilename = $candidates[0];
  print "Reading $graphdatafilename\n" if($VERBOSE == 1);
  open(GFP, "<", "$graphdatafilename") or die "Could not open file $graphdatafilename: $!";

  push(@tablevalues, {});

  my $titleline = <GFP>;
  print "titleline = $titleline\n" if($VERBOSE == 1);
  my @titles = split(/\s/, $titleline);
  my $h = shift @titles; #drop the leading '#'
  die "titleline of $graphdatafilename does not start with #" if($h ne "#");

  for(my $j=0; $j <= $#titles; $j++) {
    $tablevalues[$i]->{$titles[$j]} = [];
  }

  while(my $line = <GFP>) {
    my @values = split(/\s/, $line);
    die "Mismatch in number of columns in title = $titleline\n and line = $line\n" if ($#titles != $#values);
    for(my $j=0; $j <= $#values; $j++) {
      push(@{$tablevalues[$i]->{$titles[$j]}}, $values[$j]);
    }
  }
}

if($VERBOSE == 1) {
  print "----------- Table read --------------\n";
  for(my $i=0; $i <= $max_index; $i++) {
    print "Read data for ", $logs[$i], "\n";

    for my $k (sort (keys %{$tablevalues[$i]})) {
      print "   $k: ";
      foreach my $v (@{$tablevalues[$i]->{$k}}) {
        print "$v ";
      }
      print "\n";
    }
  }
}

# Verify that the means used in each log's graphdata are the same values
for(my $i=0; $i <= $max_index; $i++) {
  my $num_means_i = @{$tablevalues[$i]->{"mean"}};
  for(my $ii=$i+1; $ii <= $max_index; $ii++) {
    my $num_means_ii = @{$tablevalues[$ii]->{"mean"}};
    die "Mismatching number of means: log $i has $num_means_i, log $ii has $num_means_ii" if($num_means_i != $num_means_ii);

    for(my $j=0; $j < $num_means_i; $j++) {
      my $mean_i  = ${$tablevalues[$i]->{"mean"}}[$j];
      my $mean_ii = ${$tablevalues[$ii]->{"mean"}}[$j];
      die "Mismatched mean values on index $j in log $i vs $ii: $mean_i vs $mean_ii" if($mean_i ne $mean_ii);
    }
  }
}



# Identify reference columns -- scores will be computed against all reference columns of the various graphdata files.
my @ref_columns = ();

# Identify score columns -- scores will be computed only for the score columns.
my @score_columns = ();

die "Internal error" if($#logs != $max_index or $#tags != $max_index);
for(my $i=0; $i <= $max_index; $i++) {
  push(@ref_columns, []);
  push(@score_columns, []);
  my @all_log_columns = grep {$_ ne "mean"} (sort(keys %{$tablevalues[$i]}));
  if(@{$tags[$i]} == 0) { #no limiting tags specified, include all columns from i'th graphdata file
    @{$ref_columns[$i]} = (@all_log_columns);
    @{$score_columns[$i]} = (@all_log_columns);
  }
  else {
    my @selected_ref_columns = ();
    my @selected_score_columns = ();
    foreach my $k (@all_log_columns) {
      my $matches_a_ref_tag = 0;
      my $matches_a_score_tag = 0;
      foreach my $t (@{$tags[$i]}) {
        my $ts = substr $t, 1; #skip the initial + or -
        my $flag = substr $t, 0, 1; #the initial + or -
        die "$t is not a valid tag" if($flag ne "+" and $flag ne "-");
        if($k =~ /$ts/) {
          $matches_a_ref_tag = 1 if($flag eq "-");
          $matches_a_score_tag = 1 if($flag eq "+");
        }
      }
      push(@selected_ref_columns, $k) if($matches_a_ref_tag == 1);
      push(@selected_score_columns, $k) if($matches_a_score_tag == 1);
    }
    @{$ref_columns[$i]} = (@selected_ref_columns);
    @{$score_columns[$i]} = (@selected_score_columns);
  }
}

if($VERBOSE == 1) {
  print "----------- Ref columns --------------\n";
  for(my $i=0; $i <= $max_index; $i++) {
    print "Ref columns for ", $logs[$i], "\n";
    for my $k (@{$ref_columns[$i]}) {
      print "   $k\n";
    }
  }
}

my @score_values = ();

die "Internal error" if($#ref_columns != $max_index or $#score_columns != $max_index);
for(my $i=0; $i <= $max_index; $i++) {
  push(@score_values, []);

  # compute each score in current log's graphdata
  for my $k (@{$score_columns[$i]}) {
    my $min_score;
    my @closest_avg_metrics;
    my @closest_ref;
    # go over all ref columns
    for(my $ii=0; $ii <= $max_index; $ii++) {
      for my $r (@{$ref_columns[$ii]}) {
        next if($i eq $ii and $k eq $r); # skip comparison against self

        my $r_score = 0.0;
        my $num_means = @{$tablevalues[$i]->{"mean"}};
        die "Mismatched num means" if($num_means != @{$tablevalues[$ii]->{"mean"}});
        die "0 means" if($num_means == 0);
        my $avg_metric_s = 0;
        my $avg_metric_r = 0;
        for(my $mi=0; $mi < $num_means; $mi++) {
          my $vs = ${$tablevalues[$i]->{$k}}[$mi];
          my $vr = ${$tablevalues[$ii]->{$r}}[$mi];
          $r_score += $vs - $vr;
          $avg_metric_s += $vs;
          $avg_metric_r += $vr;
        }
        $r_score = $r_score / $num_means;
        $avg_metric_s = $avg_metric_s / $num_means;
        $avg_metric_r = $avg_metric_r / $num_means;
        if(not defined($min_score) or $min_score > $r_score) {
          $min_score = $r_score;
          @closest_avg_metrics = ($avg_metric_s, $avg_metric_r);
          @closest_ref = ($ii, $r);
        }
      }
    }
    die "Could not score" if(not defined($min_score));
    push(@{$score_values[$i]}, [$min_score, $closest_avg_metrics[0], $closest_avg_metrics[1], $closest_ref[0], $closest_ref[1]]);
  }
}


print "----------- Scores --------------\n" if($VERBOSE == 1);
for(my $i=0; $i <= $max_index; $i++) {
  print "Scores for ", $logs[$i], "\n" if($VERBOSE == 1);
  die "Mismatch" if($#{$score_values[$i]} != $#{$score_columns[$i]});
  for(my $j=0; $j < @{$score_values[$i]}; $j++) {
    my $k = ${$score_columns[$i]}[$j];
    my $info = ${$score_values[$i]}[$j];
    my ($score, $avg_metric_s, $avg_metric_r, $closest_ii, $closest_ref) = (@{$info});
    print "   score $k: $score, $avg_metric_s, $avg_metric_r, $closest_ii, $closest_ref\n";
  }
}
