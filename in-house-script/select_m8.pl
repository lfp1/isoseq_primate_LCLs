#!/usr/bin/perl -w
use strict;
die "Usage: <m8 file> <config> <outdir>\n" unless @ARGV == 3;
my $outdir = $ARGV[2];
mkdir $outdir unless -e $outdir;

my ($ref_sp, %spNeed) = &getContig($ARGV[1]);

my %result;
if ($ARGV[0] =~ /\.gz$/) {
	open (IN,"zcat $ARGV[0] |") 
} else {
	open IN, $ARGV[0];
}
while (<IN>) {
	my @info = split /\s+/;
	my $sp1 = (split /_/, $info[0])[0];
	my $sp2 = (split /_/, $info[1])[0];
	next unless $spNeed{$sp1} && $spNeed{$sp2};
	next if $sp1 eq $sp2;
	next unless $sp1 eq $ref_sp || $sp2 eq $ref_sp;
	($sp1, $sp2) = ($sp2, $sp1) unless $sp1 eq $ref_sp;
	push @{$result{$sp1}{$sp2}}, $_;
}
close IN;

foreach my $sp1 (keys %result) {
	foreach my $sp2 (keys %{$result{$sp1}}) {
		open OUT, ">$outdir/${sp1}_$sp2.m8";
		foreach my $line (@{$result{$sp1}{$sp2}}) {
			print OUT $line;
		}
		close OUT;
	}
}


## subroutine
######################
sub getContig {
my $contig = shift;
open IN, $contig;
my @data = <IN>; 
close IN;
my $ref_sp;
foreach my $line (@data) {
	next if $line =~ /^#/ || $line =~ /^\s+$/;
	my @info = split /\s+/, $line;
	if ($info[1] eq "ref") {
		$ref_sp = $info[0];
		last;
	}
}
my %spNeed;
foreach my $line (@data) {
	next if $line =~ /^#/ || $line =~ /^\s+$/;
	my @info = split /\s+/, $line;
	$spNeed{$info[0]} ++;
	last if $line =~ /#--END--/;
}
return ($ref_sp, %spNeed);
}
