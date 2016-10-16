#!/usr/bin/perl
use strict;
use warnings;


	open(INFILE, $ARGV[0]) or die "Cannot open $ARGV[0]: $!.\n";	
	open(GENES, $ARGV[1]) or die "Cannot open $ARGV[1]: $!.\n";
	open(OUTFILE, ">$ARGV[2]") or die "Cannot write $ARGV[2]: $!.\n";
	
	
	chomp(my @genes = <GENES>);
	close GENES;

	while(my $l = <INFILE>) {
		 if($l =~ m/^(\w+)\t/){
		 	if (($1 ~~ @genes) or ($1 eq "Gene" )) {
				 print OUTFILE "$l";
			 }
		 }
	 }
	 close INFILE;
	 close OUTFILE;
	
