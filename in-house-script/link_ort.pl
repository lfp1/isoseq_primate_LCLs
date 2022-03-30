#!/usr/bin/perl -w
use strict;
die "Usage: <ort files> <contig>\n" unless @ARGV >= 2;
my $contig = pop;

## Get reference species and the output order of species.
my @order;
my $ref_sp;
&contig_parser($contig, \@order);
$ref_sp = shift @order;

## Get ortholog pairs.
my %ort;
foreach my $ort_file (@ARGV) {
	&ort_parser($ort_file, \%ort);
}

## Output ortholog table of multiple species.
my $ref_file = $ARGV[0];
open IN, $ref_file;
while (<IN>) {
	my @info = split /\s+/;
	my $ref_gene = $info[0];
	my $out = join "\t", @info[0..5];
	foreach my $hit_sp (@order) {
		@{$ort{"$ref_sp\t$hit_sp"}{$ref_gene}} = ("NA", "NA", "NA", "NA", "NA", "NA") unless $ort{"$ref_sp\t$hit_sp"}{$ref_gene};
		$out .= "\t" . join "\t", @{$ort{"$ref_sp\t$hit_sp"}{$ref_gene}};
	}
	print "$out\n";
}
close IN;


######################################           subroutine          ####################################
sub ort_parser {
	my ($in_file, $ref) = @_;	
	open IN, $in_file;
	while (<IN>) {
		my @info = split /\s+/;
		my $ref_sp = (split /_/, $info[0])[0];
		my $hit_sp = (split /_/, $info[6])[0];
		$ref->{"$ref_sp\t$hit_sp"}{$info[0]} = [@info[6..11]];
	}
	close IN;
}


sub contig_parser {
my ($in_file, $a_ref) = @_;
open IN, $in_file;
$/ = "#--END--";
while (<IN>) {
	s/^\s+|\s+$//g;
	my @lines = split /\n/;
	foreach my $line (@lines) {
		next if $line =~ /^#/;
		my @info = split /\s+/, $line;
		push @$a_ref, $info[0];
	}
	last;
}
$/ = "\n";
close IN;
}
