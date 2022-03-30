#!usr/bin/perl -w
use strict;
=pod
unless(@ARGV==3)
{
	print "input file list:
	ref_gff        the reference gff file needed to locate the genes(including locating in which scaffold and get the order in this scaffold)
	sub_gff        the gff file of subject species
	blast_file     blast result of above species (the result should be dealed as idAdd file)\n";
	exit;
}
=cut

die "\nUsage: perl $0 ref.gff tar.gff blast.idAdd > out.ort\n\n" unless @ARGV == 3;

my $ref_gff_file=shift;
my $sub_gff_file=shift;
my $blast_file=shift;

my(%ref_sca,%ref_gene,%sub_sca,%sub_gene);
&read_gff($ref_gff_file,\%ref_sca,\%ref_gene);
&read_gff($sub_gff_file,\%sub_sca,\%sub_gene);

my %ref_gene_bg_ed;
my %sub_gene_bg_ed;
&getGenePos($ref_gff_file, \%ref_gene_bg_ed);
&getGenePos($sub_gff_file, \%sub_gene_bg_ed);

my%blast_pair=&read_blast($blast_file);

my(%seed,%seed_extend,%sub0,%sub1);
foreach my$sca_id(keys%ref_sca)
{
	my$gene_num=keys%{$ref_sca{$sca_id}};
	next unless$gene_num>=3;
	for(my$i=1;$i<=$gene_num-2;$i++)
	{
		my$ref_gene1=$ref_sca{$sca_id}{$i};
		my$ref_gene2=$ref_sca{$sca_id}{$i+1};
		my$ref_gene3=$ref_sca{$sca_id}{$i+2};

		my($flag1,$sub_gene1)=&best_check($ref_gene1,\%blast_pair);
		my($flag2,$sub_gene2)=&best_check($ref_gene2,\%blast_pair);
		my($flag3,$sub_gene3)=&best_check($ref_gene3,\%blast_pair);

		next unless $flag1 && $flag2 && $flag3;
		die "$sub_gene1" unless $sub_gene{$sub_gene1};
		my($s_sca1,$s_seq1)=@{$sub_gene{$sub_gene1}};
		my($s_sca2,$s_seq2)=@{$sub_gene{$sub_gene2}};
		my($s_sca3,$s_seq3)=@{$sub_gene{$sub_gene3}};
		
		if($s_sca1 eq $s_sca2 && $s_sca2 eq $s_sca3)
		{
			next unless $s_seq1+$s_seq3==2*$s_seq2 && ($s_seq1==$s_seq2+1||$s_seq1==$s_seq2-1);
			my$strand=($s_seq1>$s_seq2)?'-':'+';
			$seed{$sca_id}{$i}=[$s_sca1,$s_seq1,$strand];
			$sub0{$s_sca1}{$s_seq1}++;
			$sub0{$s_sca2}{$s_seq2}++;
			$sub0{$s_sca3}{$s_seq3}++;
			
			$i+=2;
		}
	}
}
foreach my$sca_id(keys%seed)
{
	my@seqs=sort{$a<=>$b}keys%{$seed{$sca_id}};
	for(my$i=0;$i<@seqs;$i++)
	{
		## the range of the seed extending
		my$start=($i==0)?1:$seqs[$i-1]+3;
		my$end=($i==@seqs-1)?(keys%{$ref_sca{$sca_id}}):$seqs[$i+1]-1;

		my$seq=$seqs[$i];
		my($sca_,$seq_,$strand)=@{$seed{$sca_id}{$seq}};

		## record the gene order in sub gene sequence
		my($bac_seq,$for_seq)=($strand eq '+')?($seq_,$seq_+2):($seq_,$seq_-2);

		## extend backward
		for(my$j=$seq-1;$j>=$start;$j--)
		{
			## slide across the ref gene sequence
			my$gene=$ref_sca{$sca_id}{$j};
			next unless exists$blast_pair{$gene};
			my@sub_genes=reverse sort{$blast_pair{$gene}{$a}->[0]<=>$blast_pair{$gene}{$b}->[0]}keys%{$blast_pair{$gene}};
			my$top_gene=$sub_genes[0];
			foreach my$sub_gene(@sub_genes)
			{
				my($s_sca,$s_seq)=@{$sub_gene{$sub_gene}};
				next unless$s_sca eq $sca_;
				if($strand eq '+'){last unless$s_seq==$bac_seq-1;}
				if($strand eq '-'){last unless$s_seq==$bac_seq+1;}
				my$lev=($sub_gene eq $top_gene)?'L1+':'L1';
				$seed_extend{$sca_id}{$j}=[$s_sca,$s_seq,$lev]unless exists$sub0{$s_sca}{$s_seq};
				$sub0{$s_sca}{$s_seq}++if$lev eq 'L1+';
				$sub1{$s_sca}{$s_seq}{$gene}++if$lev eq 'L1';
				$bac_seq-- if $strand eq '+';
				$bac_seq++ if $strand eq '-';

				last;
			}
		}

		## extend forward
		for(my $j=$seq+3;$j<=$end;$j++)
		{
			## slide across the ref gene sequence
			my $gene=$ref_sca{$sca_id}{$j};

			next unless exists$blast_pair{$gene};
			my@sub_genes=reverse sort{$blast_pair{$gene}{$a}->[0]<=>$blast_pair{$gene}{$b}->[0]}keys%{$blast_pair{$gene}};
			my$top_gene=$sub_genes[0];

			foreach my$sub_gene(@sub_genes)
			{
				my($s_sca,$s_seq)=@{$sub_gene{$sub_gene}};
				next unless$s_sca eq $sca_;
				if($strand eq '+'){last unless$s_seq==$for_seq+1;}
				if($strand eq '-'){last unless$s_seq==$for_seq-1;}
				my$lev=($sub_gene eq $top_gene)?'L1+':'L1';
				$seed_extend{$sca_id}{$j}=[$s_sca,$s_seq,$lev]unless exists$sub0{$s_sca}{$s_seq};
				$sub0{$s_sca}{$s_seq}++if$lev eq 'L1+';
				$sub1{$s_sca}{$s_seq}{$gene}++if$lev eq 'L1';
				$for_seq++ if $strand eq '+';
				$for_seq-- if $strand eq '-';

				last;
			}
		}

	}
}

foreach my$sca_id(keys%seed_extend)
{
	foreach my$seq(keys%{$seed_extend{$sca_id}})
	{
		my($s_sca,$s_seq,$lev)=@{$seed_extend{$sca_id}{$seq}};
		my$gene=$ref_sca{$sca_id}{$seq};
		next if$lev eq 'L1+';
		if($sub0{$s_sca}{$s_seq})
		{
			delete$seed_extend{$sca_id}{$seq};
		}elsif(keys %{$sub1{$s_sca}{$s_seq}}>=2){
			delete$seed_extend{$sca_id}{$seq};
			delete$sub1{$s_sca}{$s_seq}{$gene};
		}
	}
}

foreach my$sca_id(keys%seed)
{
	my$lev='L1++';
	foreach my$seq(keys%{$seed{$sca_id}})
	{
		my($s_sca,$s_seq,$s_strand)=@{$seed{$sca_id}{$seq}};
		for(my$i=0;$i<3;$i++)
		{
			$seed_extend{$sca_id}{$seq++}=[$s_sca,$s_seq,$lev];
			$s_seq++ if $s_strand eq '+';
			$s_seq-- if $s_strand eq '-';
		}
	}
}

my(%add,%sca_ref);
foreach my$sca_id(keys%ref_sca)
{
	foreach my$seq(keys %{$ref_sca{$sca_id}})
	{
		if(exists$seed_extend{$sca_id}{$seq})
		{
			my($s_sca,$s_seq,$lev)=@{$seed_extend{$sca_id}{$seq}};
			$sub0{$s_sca}{$s_seq}++;
			next;
		}
		$sca_ref{$sca_id}{$seq}++;
	}
}
foreach my$sca_id(keys%sca_ref)
{
	my@blanks;
	for(my$i=1;$i<=(keys%{$ref_sca{$sca_id}});$i++)
	{
		next unless$sca_ref{$sca_id}{$i};
		my($flag,$seq_);
		for(my$j=$i+1;;$j++)
		{
			last unless exists$sca_ref{$sca_id}{$j};
			$flag++;
			$seq_=$j;
			last if$j==(keys %{$ref_sca{$sca_id}});
		}
		
		next unless$flag;
		push@blanks,[$i,$seq_];
		$i=$seq_;
	}
	
	for(my$i=0;$i<@blanks;$i++)
	{
		my($start,$end)=@{$blanks[$i]};
		my$last_ed=$start-1;
		my$last_st=($i==0)?1:$blanks[$i-1]->[1]+1;
		my$next_st=$end+1;
		my$next_ed=($i==@blanks-1)?(keys%{$ref_sca{$sca_id}}):$blanks[$i+1]->[0]-1;
		my($last_sca1,$last_seq1)=($start==1)?('NA','NA'):@{$seed_extend{$sca_id}{$last_ed}}[0,1];
		#my($last_sca2,$last_seq2)=($start==1)?('NA','NA'):@{$seed_extend{$sca_id}{$last_st}}[0,1];
		my($next_sca1,$next_seq1)=($end==(keys%{$ref_sca{$sca_id}}))?('NA','NA'):@{$seed_extend{$sca_id}{$next_st}}[0,1];
		#my($next_sca2,$next_seq2)=($end==(keys%{$ref_sca{$sca_id}}))?('NA','NA'):@{$seed_extend{$sca_id}{$next_ed}}[0,1];
		my($last_sca,$last_seq,$next_sca,$next_seq);
		if($last_sca1 eq 'NA'){$last_sca1=$next_sca1;}
		if($next_sca1 eq 'NA'){$next_sca1=$last_sca1;}
		#my($last_seq,$next_seq)=(sort{$a<=>$b}($last_seq1,$last_seq2,$next_seq1,$next_seq2))[0,-1];
		($last_sca,$next_sca)=($last_sca1,$next_sca1);          ## $last_sca1==$last_sca2&&$next_sca1==$next_sca2;
		($last_seq,$next_seq)=($last_seq1,$next_seq1);
		
		my(%pair1,%pair2);
		for(my$j=$start;$j<=$end;$j++)
		{
			my$gene=$ref_sca{$sca_id}{$j};
			my($flag,$s_gene)=&best_check($gene,\%blast_pair);
			next unless$flag;
			my($s_sca,$s_seq)=@{$sub_gene{$s_gene}};
			if($last_sca eq $next_sca)
			{
				next unless$s_sca eq $last_sca;
				#if($last_seq eq 'NA'){next unless($s_seq>=$next_seq-10 && $s_seq<=$next_seq+10);}
				#elsif($next_seq eq 'NA'){next unless($s_seq<=$last_seq+10 && $s_seq>=$last_seq-10);}
				#else{next unless($s_seq>=$last_seq-5 && $s_seq<=$next_seq+5)||($s_seq<=$last_seq+5 && $s_seq>=$next_seq-5);}
			}
			push@{$pair1{$s_sca}},$s_seq;
			$pair2{$s_sca}{$s_seq}=$j;
		}
		if($last_sca eq $next_sca)
		{
			next unless exists$pair1{$last_sca};
			my@s_seqs=@{$pair1{$last_sca}};
			my@locs=&array_test(\@s_seqs);
			next unless@locs;
			for(my$j=0;$j<@locs;$j++)
			{
				my$s_seq=$s_seqs[$locs[$j]];
				my$seq=$pair2{$last_sca}{$s_seq};
				$add{$sca_id}{$seq}=[$last_sca,$s_seq]unless exists$sub0{$last_sca}{$s_seq};
				$sub0{$last_sca}{$s_seq}++;
			}
		}else{
			foreach my$s_sca(keys%pair1)
			{
				my@s_seqs=@{$pair1{$s_sca}};
				next unless@s_seqs>=2;
				my@locs=&array_test(\@s_seqs);
				next unless@locs>=2;
				for(my$k=0;$k<@locs;$k++)
				{
					my$s_seq=$s_seqs[$locs[$k]];
					my$seq=$pair2{$s_sca}{$s_seq};
					$add{$sca_id}{$seq}=[$s_sca,$s_seq]unless exists$sub0{$s_sca}{$s_seq};
					$sub0{$s_sca}{$s_seq}++;
				}
			}
		}
	}
}

foreach my$sca_id(sort keys %ref_sca)
{
	foreach my$seq(sort{$a<=>$b}keys%{$ref_sca{$sca_id}})
	{
		my$gene=$ref_sca{$sca_id}{$seq};
		my($s_gene,$s_sca,$s_seq);
		my($align_rate1,$align_rate2,$identity,$lev);
		if(exists$seed_extend{$sca_id}{$seq})
		{
			($s_sca,$s_seq,$lev)=@{$seed_extend{$sca_id}{$seq}};
			$s_gene=$sub_sca{$s_sca}{$s_seq};
		}elsif(exists$add{$sca_id}{$seq}){
			$lev='L1--';
			($s_sca,$s_seq)=@{$add{$sca_id}{$seq}};
			$s_gene=$sub_sca{$s_sca}{$s_seq};
		}else{
			$lev='L2';
			#$s_gene=(reverse sort{$blast_pair{$gene}{$a}->[0]<=>$blast_pair{$gene}{$b}->[0]}keys%{$blast_pair{$gene}})[0];
			my$flag;
			($flag,$s_gene)=&best_check($gene,\%blast_pair);
			if($s_gene ne 'NA')
			{
				($s_sca,$s_seq)=@{$sub_gene{$s_gene}};
				$s_gene='NA'if$sub0{$s_sca}{$s_seq};
			}
		}
		($align_rate1,$align_rate2,$identity)=@{$blast_pair{$gene}{$s_gene}}[1..3]if$s_gene ne 'NA';
		($s_gene,$s_sca,$s_seq,$align_rate1,$align_rate2,$identity,$lev)=('NA','NA','NA','NA','NA','NA','NA')if$s_gene eq 'NA';
		
		$sub_gene_bg_ed{$s_gene} = ["NA", "NA", "NA", "NA"] unless defined($sub_gene_bg_ed{$s_gene});
		print"$gene\t$sca_id\t$seq\t$ref_gene_bg_ed{$gene}->[1]\t$ref_gene_bg_ed{$gene}->[2]\t$ref_gene_bg_ed{$gene}->[3]\t$s_gene\t$s_sca\t$s_seq\t$sub_gene_bg_ed{$s_gene}->[1]\t$sub_gene_bg_ed{$s_gene}->[2]\t$sub_gene_bg_ed{$s_gene}->[3]\t$align_rate1\t$align_rate2\t$identity\t$lev\n";
	}
}

##**********************************
# read the gff file to get the gene order in the whole genome.
sub read_gff
{
	my($gff_file,$sca,$gene)=@_;

	my%gff;
	open IN,$gff_file;
	while(<IN>)
	{
		chomp;
		my@data=split /\t+/;

		next unless $data[2]eq'mRNA';

		die unless $data[8] =~ /^ID=(\S+?);/;
		my$gene_id=$1;
		$gff{$data[0]}{$gene_id}=[$data[3],$data[4]];
	}
	close IN;

	foreach my$sca_id(keys%gff)
	{
		my@sort_id=sort{$gff{$sca_id}{$a}->[0]<=>$gff{$sca_id}{$b}->[0]}keys%{$gff{$sca_id}};
		for(my$i=0;$i<@sort_id;$i++)
		{
			$$gene{$sort_id[$i]}=[$sca_id,$i+1];
			$$sca{$sca_id}{$i+1}=$sort_id[$i];
		}
	}
}

# read the blast file (idAdd format).
sub read_blast
{
	my($blast_file)=@_;

	my%blast;
	open IN,$blast_file;
	while(<IN>)
	{
		chomp;
		my@data=split /\t+/;

		$blast{$data[0]}{$data[6]}=[$data[11],$data[4],$data[10],$data[12]];
	}
	close IN;

	return %blast;	
}

# check a gene whether have a best hit gene.
sub best_check
{
	my($gene,$blast)=@_;

	my($flag,$best_gene);
	if($$blast{$gene})
	{
		$best_gene=(reverse sort{$$blast{$gene}{$a}->[0]<=>$$blast{$gene}{$b}->[0]}keys%{$$blast{$gene}})[0];
		if($$blast{$best_gene})
		{
			my$media=(reverse sort{$$blast{$best_gene}{$a}->[0]<=>$$blast{$best_gene}{$b}->[0]}keys%{$$blast{$best_gene}})[0];
			($flag,$best_gene)=(1,$best_gene)if$media eq $gene;
		}
	}
	($flag,$best_gene)=(0,'NA')unless$flag;

	return($flag,$best_gene);
}

# check a array whether arrangeed in order.
sub array_test
{
	my($array)=@_;
	my(@locs1,@locs2);
	for(my$i=0;$i<@$array;$i++)
	{
		my($flag1_1,$flag1_2,$flag2_1,$flag2_2)=(0,0,0,0);
		for(my$j=$i+1;$j<@$array;$j++)
		{
			$flag1_1++if$$array[$i]>$$array[$j];
			$flag2_1++if$$array[$i]<$$array[$j];
		}
		for(my$j=0;$j<$i;$j++)
		{
			$flag1_2++if$$array[$i]<$$array[$j];
			$flag2_2++if$$array[$i]>$$array[$j];
		}
		push@locs1,$i if($flag1_1+$flag1_2)/@$array<=0.25;
		push@locs2,$i if($flag2_1+$flag2_2)/@$array<=0.25;
	}

	return unless @locs1/@$array>=0.6||@locs2/@$array>=0.6;
	return@locs1 if@locs1/@$array>=0.6;
	return@locs2 if@locs2/@$array>=0.6;
}

sub getGenePos {
    my ($in_file, $ref) = @_;
	open IN, $in_file;
	while (<IN>) {
		my @info = split /\t/;
		#my @info = split /\s+/;
		next unless $info[2] eq "mRNA";
		die unless $info[8] =~ /^ID=(\S+?);/;
		my $id = $1;
		$ref->{$id} = [$info[0], $info[6], $info[3], $info[4]]; ## id chr strand bg ed
	}
	close IN;
}

