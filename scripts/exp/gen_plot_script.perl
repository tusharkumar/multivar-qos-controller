#!/usr/bin/perl -w

use strict;
use Cwd;
use File::Basename;

my $usage = "
Usage: gen_plot_script <data-set1> <data-set2>
   or: gen_plot_script -f <file-listing-data-sets>

Example: ./gen_plot_script 3 4 5
Example: ./gen_plot_script -f datasets.txt
where datasets.txt has following contents:
3
4
5
";

die $usage if (@ARGV < 1);

my $pwd = cwd();
my $expname = basename($pwd); 

my @data_sets;
if($ARGV[0] eq "-f") {
	my $flds = $ARGV[1];
	print "Reading data sets from file $flds\n";
	open FLDS, $flds or die "Could not open file listing datasets=$flds: $!\n";
	while(<FLDS>) {
		chomp $_;
		push(@data_sets, $_);
	}
}
else {
	@data_sets = @ARGV;
}

print "data_sets=";
foreach my $ds (@data_sets) {
	print "$ds, ";
}
print "\n";

my @allfiles = ();
foreach my $ds (@data_sets) {
	my @files = <log*_${ds}_*>;
	push(@allfiles, @files);
}

my %xset = ();
my %meanset = ();
my %wfset = ();
my %slset = ();
foreach my $f (@allfiles) {
	$xset{$1} = 1 if($f =~ /^log_([^_]+)/);
	$meanset{$1} = 1 if($f =~ /mean([\d.]+)/);
	$wfset{$1} = 1 if($f =~ /wf([\d.]+)/);
	$slset{$1} = 1 if($f =~ /sl([\d.]+)/);
}
my @xs = sort keys %xset;
my @means = sort keys %meanset;
my @wfs = sort keys %wfset;
my @sls = sort keys %slset;

foreach my $wf (@wfs) {
	foreach my $sl (@sls) {

		my @headers = ();
		my %header_hash = ();
		my %MSEQ_hash = ();
		my %SR_hash = ();
		foreach my $mean (@means) {
			foreach my $x (@xs) {
				foreach my $ds (@data_sets) {
					my $header = "${x}_${ds}";
					if(not defined $header_hash{$header}) {
						$header_hash{$header} = 1;
						push(@headers, $header);
					}
					my @candidates = <log_${x}_${ds}_mean${mean}_wf${wf}_sl${sl}*>;
					die "There should have been atmost one match: log_${x}_${ds}_mean${mean}_wf${wf}_sl${sl}*\n" if (@candidates > 1);
					my $MSEQ = -0.1;
					my $SR = -0.1;
					if(@candidates == 1) {
						my $logfile = $candidates[0];
						open LF, $logfile or die "Could not open logfile=${logfile}: $!\n";
						while(<LF>) {
							($MSEQ, $SR)=($1, $2) if(/STATS:.* cumulative_MSEQ=([\d.]+).* cumulative_satisfaction_ratio=([\d.]+)/);
						}
						print "logfile=${logfile} has MSEQ=${MSEQ} SR=${SR}\n";
					}
					$MSEQ_hash{$mean}{$header} = $MSEQ;
					$SR_hash{$mean}{$header} = $SR;
				}
			}
			print "\n";
		}

		## Filter
		my @keep_header_indices = ();
		for(my $i=0; $i<@headers; $i++) {
			my $h = $headers[$i];
			my $count = 0;
			foreach my $mean (@means) {
				$count++ if $SR_hash{$mean}{$h} >= 0;
			}
			push(@keep_header_indices, $i) if $count > 0;
		}
		print "all headers=@headers\n";
		print "keep_header_indices=@keep_header_indices\n";
		@headers = @headers[@keep_header_indices];
		print "filtered headers=@headers\n";

		## Output
		my @metrics = ("sr", "mseq");
		foreach my $m (@metrics) {
			#my $plotfile = "graphdata_ds" . join("-",@data_sets) . "_" . "wf${wf}_sl${sl}.gnuplot";
			my $plotfile = "graphdata_" . $m . "_ds" . "all-" . $data_sets[-1] . "_" . "wf${wf}_sl${sl}.gnuplot";
			open PF, "> $plotfile" or die "Could not open plotfile=$plotfile: $!\n";
			my $header_string = "# mean\t" . join("\t", @headers);
			print PF $header_string . "\n";
			foreach my $mean (@means) {
				my $outstring = "${mean}";
				foreach my $h (@headers) {
					my $MSEQ = $MSEQ_hash{$mean}{$h};
					my $SR = $SR_hash{$mean}{$h};
					my $mval = ($m eq "sr" ? $SR : log($MSEQ)/log(10));
					$outstring = $outstring . "\t${mval}";
				}
				print PF $outstring . "\n";
			}

			#my $scriptfile = "script_ds" . join("-",@data_sets) . "_" . "wf${wf}_sl${sl}.gnuplot";
			my $scriptfile = "script_" . $m . "_ds" . "all-" . $data_sets[-1] . "_" . "wf${wf}_sl${sl}.gnuplot";
			open SF, "> $scriptfile" or die "Could not open scriptfile=$scriptfile: $!\n";
			my $command = "set title '${m} ${expname}'\n";
			$command .= "plot ";
			for(my $i=0; $i<@headers; $i++) {
				$command .= ", " if($i > 0);
				my $column = $i+2;
				my $h = $headers[$i];
				$command .= "\"${plotfile}\" using 1:${column} title '${h}' with linespoints"
			}
			$command .= " lw 3"; #make the last line thicker
			print SF $command;
		}
	}
}
