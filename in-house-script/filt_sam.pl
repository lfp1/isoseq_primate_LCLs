#!/usr/bin/perl -w
use strict;

if ($ARGV[0]=~/.bam$/){
	open IN, "samtools view $ARGV[0] |";
}
else{
	open IN, "<$ARGV[0]";
}
while(<IN>){
	chomp;
	if ($_=~/^transcript/){
		my @a=split/\s+/;
		next if ($a[2] eq 'chrM');
		my $tag=0;
		foreach my $i(@a){
			$tag=1 if ($i eq 'NH:i:1');
		}
		next unless ($tag==1);
		print "$_\n";
	}
	else{
		next if ($_=~/chrM/);
		print "$_\n";
	}
}
close IN;
