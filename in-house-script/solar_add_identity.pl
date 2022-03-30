#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my ($solarFIle, $tabFile, $m8File, $best);
my ($Verbose,$Help);
GetOptions(
	"solar:s"=>\$solarFIle,
	"tab:s"=>\$tabFile,
	"m8:s"=>\$m8File,
	"best"=>\$best,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
if (!$solarFIle || (!$solarFIle && $tabFile) || (!$solarFIle && $m8File) || $Help) {
	 print <<"Usage End.";

Usage: 
perl solar_add_identity.pl --solar Foc1_miss.pep.blast.solar --m8 ../Foc1_miss.pep.blast.m8 >Foc1_miss.pep.blast.solar.idAdd
perl solar_add_identity.pl --solar Foc1_miss.pep.blast.solar --m8 ../Foc1_miss.pep.blast.m8 --best >Foc1_miss.pep.blast.solar.idAdd.best

Usage End.
	exit;
}

my %identity;
###############################################################################
## get the identity of each record
if ($tabFile) {
	open IN, $tabFile;
	while (<IN>) {
		next if /Query_id/;
		my @info = split /\s+/;
		my $key = "$info[0]\t$info[4]\t$info[2]\t$info[3]\t$info[6]\t$info[7]";
		$identity{$key} = $info[8];
	}
	close IN;
} elsif ($m8File) {
	open IN, $m8File;
	while (<IN>) {
		my @info = split /\s+/;
		($info[6], $info[7]) = ($info[7], $info[6]) if $info[6] > $info[7];
		($info[8], $info[9]) = ($info[9], $info[8]) if $info[8] > $info[9];
		my $key = "$info[0]\t$info[1]\t$info[6]\t$info[7]\t$info[8]\t$info[9]";
		$identity{$key} = $info[2];
	}
	close IN;
}
##############################################################################


################################################################################
##
open IN, $solarFIle;
my @data = <IN>;
close IN;
# Qname Qlen Qstart Qstop strand Sname Slen Sstart Sstop #blocks total_score \
#                Qstart,Qstop;...; Sstart,Sstop;...; score;...;
my %line_to_score;
my %line_to_query;
foreach (@data) {
	my @info = split /\s+/;
	#my $score = &get_score($info[-1]);
	$line_to_query{$_} = $info[0];
	#$line_to_score{$_} = $score;
	$line_to_score{$_} = $info[10];
}

my %printTime;
print "#Query_id\tQuery_length\tQuery_start\tQuery_end\tQ_align_ratio\tStrand\tSubject_id\tSubject_length\tSubject_start\tSubject_end\tS_align_ratio\tScore\tIdentity\n" if $tabFile || $m8File;
foreach (sort by_QS_score @data) {
	my @info = split /\s+/;
	if ($best) {
		next if $printTime{$info[0]};
	}
	if ($tabFile || $m8File) {
		my @queryPos = &get_pos($info[11], $info[4]);
		my @subjectPos = &get_pos($info[12], $info[4]);
		my @identity;
		my ($query_total_match_len, $query_total_len, $average_identity, $query_aligning_ratio);
		my ($subject_total_len, $subject_aligning_ratio);
		for (my $i =0; $i < @queryPos; $i+=2) {
			my $key = "$info[0]\t$info[5]\t$queryPos[$i]\t" . $queryPos[$i+1] . "\t$subjectPos[$i]\t" . $subjectPos[$i+1];
			#print "$key\n";
			push @identity, $identity{$key};
			$query_total_match_len += (abs($queryPos[$i+1] - $queryPos[$i]) + 1) * $identity{$key};
			$query_total_len += abs($queryPos[$i+1] - $queryPos[$i]) + 1;
			$subject_total_len += abs($subjectPos[$i+1] - $subjectPos[$i]) + 1;
		}
		$average_identity = sprintf "%.2f", $query_total_match_len / $query_total_len;
		#$query_aligning_ratio = sprintf "%.2f", $query_total_len / $info[1];
		$query_aligning_ratio = sprintf "%.2f", ($info[3]-$info[2]+1) / $info[1];
		#$subject_aligning_ratio = sprintf "%.2f", $subject_total_len / $info[6];
		$subject_aligning_ratio = sprintf "%.2f", ($info[8]-$info[7]+1) / $info[6];
		#print $outPut;
		print "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$query_aligning_ratio\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t$subject_aligning_ratio\t$info[10]\t$average_identity\n";
	} else {
		print;
	}

	$printTime{$info[0]} ++;
}
##################################################################################



## subroutine
#################
sub get_score {
my $input = shift;
$input =~ s/\+|-//g;
my @scores = split /;/, $input;
my $total_score;
foreach my $score (@scores) {
	$total_score += $score;
}
return $total_score;
}
################

#################
sub get_pos {
my $input = shift;
my $strand = shift;
my @info = split /;/, $input;
my @pos;
foreach my $bg_ed (@info) {
	my ($bg, $ed) = split /,/, $bg_ed;
	if ($strand eq "-") {
		($bg, $ed) = ($ed, $bg) if $bg > $ed;
	}
	push @pos, ($bg, $ed);
}
return @pos;
}
###############

##################
sub by_QS_score {
	$line_to_query{$a} cmp $line_to_query{$b}
			or
	$line_to_score{$b} <=> $line_to_score{$a}
}
##################
