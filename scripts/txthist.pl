#!/usr/bin/perl -w

use strict;

my %hist;

while(<>) {
	chomp;
	$hist{$_}++;
}

my $firstkey = (keys %hist)[0];
my $maxfrequency = $hist{$firstkey};
foreach my $k (keys %hist) {
	$maxfrequency = $hist{$k} if($maxfrequency < $hist{$k});
}
print "maxfrequency=$maxfrequency\n";

foreach my $k (sort {$a <=> $b} keys %hist) {
	print "$k : " . ("X" x ($hist{$k} * 50.0 / $maxfrequency)) . " $hist{$k}\n";
}
