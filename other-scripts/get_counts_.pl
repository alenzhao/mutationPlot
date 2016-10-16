#!/usr/bin/perl
use strict;
use warnings;


	open(INFILE, $ARGV[0]) or die "Cannot open $ARGV[0]: $!.\n";
	open(OUTFILE, ">$ARGV[1]") or die "Cannot write $ARGV[1]: $!.\n";
	open(GENES, $ARGV[2]) or die "Cannot open $ARGV[2]: $!.\n";
	
	chomp(my @genes = <GENES>);
	
	print OUTFILE "gene\ttranscript\tposition\tchromosome\tcounts\tannotation\tLOF\n";
	while(my $l = <INFILE>) {
		next if $l =~ m/^#/;
		my @fields = split("\t", $l);
		if($l =~ m/AN_NFE=(\d+).*CSQ=(.*(ENSG\d+).*(ENST\d+).*)/){
			my $csq = $2;
        	my  @a =split(/\|/,$csq);
			my $gene=$3;
			my $transcript=$4;
			my $counts=$1;
			
			if($gene ~~ @genes){
				print OUTFILE "$gene\t$transcript\t$fields[1]\t$fields[0]\t$counts\t$a[1]\t$a[56]\n";
			}
		}
	}
	close INFILE;
	close OUTFILE;
	
