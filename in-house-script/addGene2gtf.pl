#!/usr/bin/perl -w
use strict;
use List::Util qw/sum/;

my %hash;
open IN, "<$ARGV[0]";
while(<IN>){
	chomp;
	next if ($_=~/^isoform/);

	my @a=split/\s+/;
	my ($gid,$tid)=($a[6],$a[0]);

	$hash{$tid}=$gid;
}
close IN;

my (%chr,%std,%bg,%ed,%exo);
open IN, "<$ARGV[1]";
while(<IN>){
	chomp;
	my $tid=$1 if ($_=~/transcript_id "(\S+)";/);

	my @a=split/\t+/;
	my ($chr,$bg,$ed,$std)=($a[0],$a[3],$a[4],$a[6]);

	$chr{$tid}=$chr;
	$std{$tid}=$std;
	$bg{$tid}=$bg if (!exists$bg{$tid} || $bg{$tid}>$bg);
	$ed{$tid}=$ed if (!exists$ed{$tid} || $ed{$tid}<$ed);

	my $elo="$chr,$std,$bg,$ed";
	$exo{$tid}{$elo}=[$chr,$std,$bg,$ed];
}
close IN;

foreach my $tid(keys%exo){
	my ($chr,$bg,$ed,$std)=($chr{$tid},$bg{$tid},$ed{$tid},$std{$tid});
	my $gid=$hash{$tid};

	print "$chr\tPacBio\ttranscript\t$bg\t$ed\t\.\t$std\t\.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";

	if ($std eq '+'){
		foreach my $elo (sort {$exo{$tid}{$a}[2] <=> $exo{$tid}{$b}[2]} keys %{$exo{$tid}}){
			my ($bgt,$edt)=($exo{$tid}{$elo}[2],$exo{$tid}{$elo}[3]);
			print "$chr\tPacBio\texon\t$bgt\t$edt\t\.\t$std\t\.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
		}
	}
	else{
		foreach my $elo (sort {$exo{$tid}{$b}[3] <=> $exo{$tid}{$a}[3]} keys %{$exo{$tid}}){
			my ($bgt,$edt)=($exo{$tid}{$elo}[2],$exo{$tid}{$elo}[3]);
			print "$chr\tPacBio\texon\t$bgt\t$edt\t\.\t$std\t\.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
		}
	}
}
