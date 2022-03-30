#!/usr/bin/perl -w
use strict;

my %del;

if ($ARGV[0]=~/.sam$/){
	open IN, "<$ARGV[0]";
}
else{
	open IN, "samtools view $ARGV[0] |";
}
while(<IN>){
	chomp;
	next unless ($_=~/^transcript/);
	my @a=split/\s+/;
	my $tag=0;
	foreach my $i(@a){
		$tag=1 if ($i eq 'NH:i:1');
	}
	$tag=0 if ($a[2] eq 'chrM');
	next unless ($tag==0);
	$del{$a[0]}=1;
}
close IN;

if ($ARGV[1]=~/.gz$/){
	open IN, "gzip -dc $ARGV[1] |";
}
else{
	open IN, "<$ARGV[1]";
}
local $/ = ">";
while(<IN>){
	chomp;
	next if ($_ eq '');
	my @a=split/\n/,$_;
	my $info=shift@a;
	my $seq=join"",@a;

	my @b=split/\s+/,$info;
	my $id=$b[0];
	next if (exists$del{$id});
	
	chomp($id,$seq);
	print ">$id\n$seq\n";
}
close IN;

