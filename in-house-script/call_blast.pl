#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
die "\nUsage: <query_dir> <database_file> <blast type> <out_dir> <suffix of query> <memory for qsub(G)> <-q for qsub (e.g st.q)> <-P for qsub (e.g P17Z10200N0101)>\n\n" unless @ARGV == 8;
my $query_dir = shift;
my $database_file = shift;
my $blast_type = shift;
my $out_dir = shift;
my $suffix = shift;
my $mem = shift;
my $qsub_q = shift;
my $qsub_p = shift;
#my $program = "/share/project002/liqiye/bin/genBlastA/blastall";
#my $program = "/opt/blc/genome/biosoft/blast-2.2.23/bin/blastall";
my $program = "/share/app/blast-2.2.26/bin/blastall";
mkdir $out_dir unless -e $out_dir;

#foreach my $p (\$query_dir, \$database_file, \$out_dir) {
#	$$p = abs_path($$p);
#}

opendir IN, $query_dir;
while (my $file = readdir IN) {
	next unless $file =~ /.+\.$suffix$/;
	my $query_file = "$query_dir/$file";
	my $out_file = "$out_dir/$file.m8";
	my $sh_file = $file =~ /^\d+/ ? "$out_dir/blast_$file.sh" : "$out_dir/$file.sh";
	open SH, ">$sh_file";
	#print SH "date\n$program -p $blast_type -i $query_file -d $database_file -o $out_file -F F -a 4 -e 1e-5 -m 8\ndate\n";
	print SH "date\n$program -p $blast_type -i $query_file -d $database_file -o $out_file -F F -a 4 -e 1e-2 -m 8\ndate\n";
    close SH;
	system "qsub -S /bin/sh -cwd -l vf=${mem}G,num_proc=1 -q $qsub_q -P $qsub_p -binding linear:4 $sh_file";
}
closedir IN;
