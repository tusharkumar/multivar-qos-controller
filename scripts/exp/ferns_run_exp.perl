#!/usr/bin/perl -w

use strict;
use File::Basename;

my %explore_space;
$explore_space{"mean"} = [0.06, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09];
$explore_space{"window_frac"} = [0.20, 0.10, 0.40];
$explore_space{"sliding_window"} = [4];
$explore_space{"srt"} = ["onlysrt", "fixedcases"];
$explore_space{"test"} = ["testing", "running"];


my $usage = "
Usage: run_exp.perl <log_number> <num_objs> [srt] [mean] [window_frac] [sliding_window] [test]
	If one or more of mean, window_frac, sliding_window are specified, then only the default case for that aspect will be run.

	srt => no fixed-level cases will be run
	test => show lists of commands, but does not execute any
";

die $usage if (@ARGV < 2);

for(my $i=2; $i<@ARGV; $i++) {
	my $arg = $ARGV[$i];
	my $argchoice;
	($arg, $argchoice) = ($1, $2) if($arg =~ /(.*)=(.*)/);
	die "Unknown option: $arg\n$usage" if(not defined $explore_space{$arg});
	$explore_space{$arg}[0] = $argchoice if(defined $argchoice); #change the default if choice specified
	$explore_space{$arg} = [ $explore_space{$arg}[0] ]; #limit to just the default case
	print "## Limiting exploration of \"$arg\" to " . $explore_space{$arg}[0] . "\n";
}

die if(@{$explore_space{"srt"}} == 0);
die if(@{$explore_space{"srt"}} > 2);
my $needs_srt = ($explore_space{"srt"}[0] eq "onlysrt" || $explore_space{"srt"}[1] eq "onlysrt" ? 1 : 0);
my $needs_fixedcases = ($explore_space{"srt"}[0] eq "fixedcases" || $explore_space{"srt"}[1] eq "fixedcases" ? 1 : 0);
$explore_space{"srt"} = [];

if($needs_srt == 1) {
	$explore_space{"srt"} = [-1];
}

if($needs_fixedcases == 1) {
	my $N = 2;
	foreach my $i (0..(2*$N)) {
		push @{$explore_space{"srt"}}, $i;
	}
}
die if(@{$explore_space{"srt"}} == 0);


###
my $log_number = $ARGV[0];
my $num_objs   = $ARGV[1]+0; #force treatment as number
print "log_number = $log_number num_objs=$num_objs\n";
die "invalid num_objs=$ARGV[1]" if ($num_objs ne $ARGV[1]);

my $extract_time_sequence_script = dirname($0) . "/scripts/extract_time_sequence.pl";

system "killall -9 Xvfb";
system "Xvfb :1 -auth /dev/null &";

foreach my $level (@{$explore_space{"srt"}}) {
	foreach my $mean (@{$explore_space{"mean"}}) {
		foreach my $window_frac (@{$explore_space{"window_frac"}}) {
			foreach my $sliding_window (@{$explore_space{"sliding_window"}}) {

my $level_string;
if($level >= 0)
{ $level_string = "fixed$level"; }
else
{ $level_string = "srt"; }

my $suffix = "_${level_string}_${log_number}_mean${mean}_wf${window_frac}_sl${sliding_window}";

my $command="DISPLAY=:1 MEAN_EXEC_TIME=$mean DEST_DIR=. ADVANCED_CONTROLLER=false LOWER_FRAC=${window_frac} UPPER_FRAC=${window_frac} DEFAULT_VALUE=${level} NUM_OBJS=${num_objs} FILENAME=info${suffix} ../../../ferns_demo -v input.mp4 -k 100 > log${suffix} 2>&1";
print "--- Invoking --- $suffix\n";
print "$command\n";

if(@{$explore_space{"test"}} > 1) {
	system $command;

	system "perl ${extract_time_sequence_script} log${suffix}";
	system "rm -f *.jpg";
}

			}
		}
	}
}

system "killall -9 Xvfb";
