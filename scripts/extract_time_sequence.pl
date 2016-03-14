#!/usr/bin/perl -w

use strict;

my $usage =
"extract_time_sequence.pl <log-file-name>
   <log-file-name> must begin with 'log_'.
   Two output files will be produced with same name as log-file,
     except 'script_' and 'graphdata_' will replace the 'log_' prefix.
";

die $usage if(@ARGV != 1);

my $logfilename = $ARGV[0];
print "logfilename=$logfilename\n";

open LFP, $logfilename or die "Could not open log-file: $logfilename: $!";

my $scriptfilename = $logfilename;
$scriptfilename =~ s/log_/timeline_script_/ or die "extract_time_sequence.pl: ERROR: log-file-name must start with 'log_'\n";
open SFP, "> $scriptfilename" or die "extract_time_sequence.pl: ERROR: Could not open script-file for writing: $!\n";

my $graphdatafilename = $logfilename;
$graphdatafilename =~ s/log_/timeline_graphdata_/ or die "extract_time_sequence.pl: ERROR: log-file-name must start with 'log_'\n";
open GFP, "> $graphdatafilename" or die "extract_time_sequence.pl: ERROR: Could not open graphdata-file for writing: $!\n";


my %event_index;
my $index_count = 0;
$event_index{"APPLIED M(new)"} = $index_count++;
$event_index{"REFINED Mp(balance)"} = $index_count++;
$event_index{"REFINED Mp(better)"} = $index_count++;
$event_index{"LLSE NO SOL"} = $index_count++;
$event_index{"APPLIED M(balanced)"} = $index_count++;
$event_index{"APPLIED M(better)"} = $index_count++;
$event_index{"LQR SOLUTION"} = $index_count++;
$event_index{"RESET M,C"} = $index_count++;
$event_index{"RESET C(newM)"} = $index_count++;
$event_index{"FORCED EXPLORATION"} = $index_count++;
$event_index{"RESCALING for Oscillation"} = $index_count++;
$event_index{"RESCALING for Sluggishness"} = $index_count++;

$event_index{"q"} = $index_count++;
$event_index{"advMpM"} = $index_count++;
$event_index{"l_tr"} = $index_count++;
$event_index{"SR_M"} = $index_count++;
$event_index{"coverage_met"} = $index_count++;
$event_index{"log2 |H|"} = $index_count++;

my @event_name_from_index = map {"undefined"} 0..($index_count-1);
for my $e (keys %event_index) {
	@event_name_from_index[ $event_index{$e} ] = $e;
}

my $frame_number;
my $xs;
my $ys;
my $num_x;
my $num_y;
my $num_frames;
my @eventflags = map { -1 } 1..$index_count;
while(<LFP>) {
	if(/write_controlparameters\(\): frame_number=(\d+)/) {
		$frame_number = $1;
	}
	if(/APPLYING frame-choices=\[(.*)\]/) {
		$num_x = split(/,/, $1);
		$xs = join(" ", split(/,/, $1));
	}

	if(/----/) { #for efficiency
		for my $e (keys %event_index) {
			if(/\Q$e\E/) {
				my $index = $event_index{$e};
				$eventflags[ $index ] = $index+4;
			}
		}
	}
	elsif(/ q=(\d*(\.)?(\d+))/) {
		$eventflags[ $event_index{"q"} ] = $1;
	}
	elsif(/advantage_Mp_over_M=((-)?\d*(\.)?(\d+))/) {
		$eventflags[ $event_index{"advMpM"} ] = $1;
	}
	elsif(/l_tr=(\d+) .*estimated_SR_M=(\d*(\.)?(\d+))/) {
		$eventflags[ $event_index{"l_tr"} ] = $1;
		$eventflags[ $event_index{"SR_M"} ] = $2;
	}
	elsif(/coverage_met=1/) {
		$eventflags[ $event_index{"coverage_met"} ] = 2.5;
	}
	elsif(/history_len=(\d+)/ and $1 > 0) {
		$eventflags[ $event_index{"log2 |H|"} ] = log($1)/log(2);
	}

	if(/READ frame-normalized-observations=\[(.*)\]/) {
		$num_y = split(/,/, $1);
		$ys = join(" ", split(/,/, $1));
		my $flags = join(" ", @eventflags);
		print GFP "$frame_number\t\t$xs\t\t$ys\t\t$flags\n";
		$num_frames = $frame_number;
		@eventflags = map { -1 } 1..$index_count;
	}
}
print "num_x=$num_x\n";
print "num_y=$num_y\n";

print SFP "set tmargin 1\n";
print SFP "set bmargin 1\n";
print SFP "set lmargin 3\n";
print SFP "set rmargin 3\n";
print SFP "set multiplot layout 3,1\n";

my $plot_command;

print SFP "set yrange [0:*]\n";
$plot_command = "plot ";
for(my $i=1; $i<=$index_count; $i++) {
	$plot_command .= ", " if($i > 1);
	my $column = $i+1+$num_x+$num_y;
	my $event_name = $event_name_from_index[$i-1];
	$plot_command .= "\"${graphdatafilename}\" using 1:${column} title '${event_name}' with points ps 0.8";
}
print SFP "$plot_command\n\n";

print SFP "set yrange [*:*]\n";
$plot_command = "plot ";
for(my $i=1; $i<=$num_y; $i++) {
	$plot_command .= ", " if($i > 1);
	my $column = $i+1+$num_x;
	$plot_command .= "\"${graphdatafilename}\" using 1:${column} title 'y${i}' with linespoints ps 0.8";
}
print SFP "$plot_command\n\n";

$plot_command = "plot ";
for(my $i=1; $i<=$num_x; $i++) {
	$plot_command .= ", " if($i > 1);
	my $column = $i+1;
	$plot_command .= "\"${graphdatafilename}\" using 1:${column} title 'x${i}' with linespoints ps 0.8";
}
print SFP "$plot_command\n\n";

print SFP "unset multiplot\n";
