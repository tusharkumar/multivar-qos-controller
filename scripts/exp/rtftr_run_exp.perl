#!/usr/bin/perl -w

use strict;
use File::Basename;

my %explore_space;
#$explore_space{"mean"} = [0.20, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.21, 0.22];
$explore_space{"mean"} = [0.20, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.22, 0.24, 0.26, 0.28, 0.30];
$explore_space{"window_frac"} = [0.20, 0.10, 0.40];
$explore_space{"sliding_window"} = [4];
$explore_space{"srt"} = ["onlysrt", "fixedcases"];
$explore_space{"test"} = ["testing", "running"];


my $usage = "
Usage: run_exp.perl <log_number> <scalefactor|minrect|both> [srt] [mean] [window_frac] [sliding_window] [test]
	If one or more of mean, window_frac, sliding_window are specified, then only the default case for that aspect will be run.

	srt => no fixed-level cases will be run
	test => show lists of commands, but does not execute any
";

die $usage if (@ARGV < 2);
die $usage if ($ARGV[1] ne "minrect" && $ARGV[1] ne "scalefactor" && $ARGV[1] ne "both");
my $use_scale_factor = ($ARGV[1] eq "scalefactor" || $ARGV[1] eq "both" ? "true" : "false");
my $use_min_rect     = ($ARGV[1] eq "minrect"     || $ARGV[1] eq "both" ? "true" : "false");

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
	if($use_scale_factor eq "true" && $use_min_rect eq "true") {
		$explore_space{"srt"} = [[-1,-1]];
	}
	else {
		$explore_space{"srt"} = [-1];
	}
}

if($needs_fixedcases == 1) {
	my $N = 2;
	if($use_scale_factor eq "true" && $use_min_rect eq "true") {
		foreach my $i (0..(2*$N)) {
			foreach my $j (0..(2*$N)) {
				push @{$explore_space{"srt"}}, [$i, $j];
			}
		}
	}
	else {
		die if($use_scale_factor eq "false" && $use_min_rect eq "false");
		foreach my $i (0..(2*$N)) {
			push @{$explore_space{"srt"}}, $i;
		}
	}
}
die if(@{$explore_space{"srt"}} == 0);


###
my $log_number = $ARGV[0];
print "log_number = $log_number\n";

my $extract_time_sequence_script = dirname($0) . "/scripts/extract_time_sequence.pl";

system "killall -9 Xvfb";
system "Xvfb :1 -auth /dev/null &";

foreach my $level (@{$explore_space{"srt"}}) {
	foreach my $mean (@{$explore_space{"mean"}}) {
		foreach my $window_frac (@{$explore_space{"window_frac"}}) {
			foreach my $sliding_window (@{$explore_space{"sliding_window"}}) {

my $scale_factor_level;
my $min_rect_level;

my $level_string;
if($use_scale_factor eq "true" && $use_min_rect eq "true") {
	$scale_factor_level = ${$level}[0];
	$min_rect_level     = ${$level}[1];
	if($scale_factor_level == -1) {
		die if($min_rect_level != -1);
		$level_string = "srt";
	}
	else
	{ $level_string = "fixed$scale_factor_level:$min_rect_level"; }
}
else {
	if($level >= 0)
	{ $level_string = "fixed$level"; }
	else
	{ $level_string = "srt"; }
	$scale_factor_level = $level;
	$min_rect_level     = $level;
}

my $suffix = "_${level_string}_${log_number}_mean${mean}_wf${window_frac}_sl${sliding_window}";

my $command="DISPLAY=:1 USE_SCALE_FACTOR=${use_scale_factor} USE_MIN_RECT=${use_min_rect} mean_exec_time=$mean DEST_DIR=. ADVANCED_CONTROLLER=false lower_frac=${window_frac} upper_frac=${window_frac} scale_factor_default_value=${scale_factor_level} min_rect_default_value=${min_rect_level} FILE1=info${suffix} SOURCE_FILE=input.avi ../../../rtftr > log${suffix} 2>&1";
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
