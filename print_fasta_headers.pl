#!/usr/bin/perl

use 5.010;
use strict;
use warnings;

open (QUERY,"<$ARGV[0]") or die $!;

#my $negregex = hypothetical;

while (my $line = <QUERY>) {
	if (($line =~m/^>/) # String to match. In this case, matches all fasta headers as they begin the line (^) with ">".

#&& ($line =~ m/^(?:(?!$negregex).)*$/) # String not to match. Excludes any lines containing "hypothetical" as this is the negative regular expression specified in the earlier code negregex variable.

)
	 { print $line ; }
	
}
