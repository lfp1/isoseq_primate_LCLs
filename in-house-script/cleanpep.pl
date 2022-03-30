#!/usr/bin/perl -w
use strict;

my %hash;

open IN, "<$ARGV[0]";
while(<IN>){
	chomp;
	next if ($_=~/^#/);
	my @a=split/\t+/;
	next unless ($a[1] !~ "pseudogene" && $a[2] eq "CDS");
	my $gid=$1 if ($a[-1]=~/gene_id "(\S+)";/);
	my $pid=$1 if ($a[-1]=~/protein_id "(\S+)";/);
	###print "$gid\t$pid\n";
	if ($gid && $pid){
		$hash{$pid}=$gid;
	}
}
close IN;

open IN, "<$ARGV[1]";
while(<IN>){
	chomp;
	if ($_=~/^>/){
		s/>//;
		my @a=split/\s+/;
		my $pid=$a[0];
		die if (!exists$hash{$pid});
		print ">$pid gene:$hash{$pid}\n";
	}
	else{
		print "$_\n";
	}
}
close IN;
