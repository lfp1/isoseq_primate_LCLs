#Correcting Iso-seq transcripts with all available pair-end RNA-seq data for identical species by LoRDEC, which can be installed with conda using the command 'conda create -n loredc env lordec'
#Because high-QV and low-QV isoforms differ only in predicted accuracy, both will be used in the following analysis
less Macaca_mulatta.polished.*.fasta.gz > Macaca_mulatta.polished.fa
/home/bin/lordec-correct -t 5 -b 200 -e 0.4 -s 3 -k 18 -T 50 -i Macaca_mulatta.polished.fa -2 /data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.1.fq.gz,/data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.2.fq.gz,/data/rnaseq/Macaca_mulatta/R02027/rep2/Macaca_mulatta.R02027.rep2.1.fq.gz,/data/rnaseq/Macaca_mulatta/R02027/rep2/Macaca_mulatta.R02027.rep2.2.fq.gz,/data/rnaseq/Macaca_mulatta/R02027/rep3/Macaca_mulatta.R02027.rep3.1.fq.gz,/data/rnaseq/Macaca_mulatta/R02027/rep3/Macaca_mulatta.R02027.rep3.2.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep1/Macaca_mulatta.R05040.rep1.1.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep1/Macaca_mulatta.R05040.rep1.2.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep2/Macaca_mulatta.R05040.rep2.1.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep2/Macaca_mulatta.R05040.rep2.2.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep3/Macaca_mulatta.R05040.rep3.1.fq.gz,/data/rnaseq/Macaca_mulatta/R05040/rep3/Macaca_mulatta.R05040.rep3.2.fq.gz -o Macaca_mulatta.correct.fa

#Trimming ployA tail remains in corrected transcript using HMM
#Installation instructions for trim_isoseq_polyA programs can be found at https://github.com/bowhan/trim_isoseq_polyA
/home/bin/trim_isoseq_polyA-master/bin/trim_isoseq_polyA -i Macaca_mulatta/Macaca_mulatta.correct.fa -t 1 -G > Macaca_mulatta.correct.trim.fa

#Aligning transcripts to reference genome by GMAP and sorting the output
/home/bin/gmap-20170215/bin/gmap -D /ref/genome/gmapindex -d Macaca_mulatta -f samse -t 50 -n 0 Macaca_mulatta.correct.trim.fa > Macaca_mulatta.GMAP.tmp.sam
less Macaca_mulatta.GMAP.tmp.sam | awk '$_!~/^transcript/' > header
less Macaca_mulatta.GMAP.tmp.sam | awk '$_~/^transcript/' > sam
cat header | sort | uniq | sort -k 3,3 -k 4,4n > header2
sort -k 3,3 -k 4,4n sam > sam2
cat header2 sam2 > Macaca_mulatta.GMAP.sam
rm Macaca_mulatta.GMAP.tmp.sam header header2 sam sam2

#Correcting Iso-seq transcripts by genome sequence and all splice junctions (SJ) identified from RNA-seq mapping with TranscriptClean
#Installation instructions for TranscriptClean can be found at https://github.com/mortazavilab/TranscriptClean
python /home/bin/TranscriptClean-master/TranscriptClean.py --sam Macaca_mulatta.GMAP.sam --genome /ref/genome/Macaca_mulatta.fa --spliceJns /analysis_RNA/01.align-genome/Macaca_mulatta/all.SJ.out.tab --outprefix ./

#Retainning the uniquely-mapped reads and remove the reads mapped to chrM by in-house scirpts
perl /home/in-house-script/filt_fa.pl TC_clean.sam TC_clean.fa > Macaca_mulatta.correct.trim.rm.fa
perl /home/in-house-script/filt_sam.pl TC_clean.sam > Macaca_mulatta.GMAP.rm.sam

#TOFU cupcake was used to collapse transcripts
#Installation instructions for cDNA_Cupcake can be found at https://github.com/Magdoll/cDNA_Cupcake
#Collapsing HQ isoform results to unique isoforms (based on genome alignment)
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/collapse_isoforms_by_sam.py --input Macaca_mulatta.correct.trim.rm.fa -s Macaca_mulatta.GMAP.rm.sam -o Macaca_mulatta.all
#Converting fasta format to fastq format
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/sequence/fa2fq.py Macaca_mulatta.all.collapsed.rep.fa; mv Macaca_mulatta.all.collapsed.rep.fastq Macaca_mulatta.all.collapsed.rep.fq
#Obtaining count information post collapse to unique isoforms
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/get_abundance_post_collapse.py Macaca_mulatta.all.collapsed Macaca_mulatta.polished.cluster_report.csv
#Filtering away 5' degraded isoforms
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/filter_away_subset.py Macaca_mulatta.all.collapsed
#Preparing file for running subsampling (rarefaction curve)
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/annotation/make_file_for_subsampling_from_collapsed.py -i Macaca_mulatta.all.collapsed.filtered -o Macaca_mulatta.all.collapse.filtered.for_subsampling.txt
#Finding fusion genes
/home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/fusion_finder.py --input Macaca_mulatta.correct.trim.fa -s TC_clean.sam -o Macaca_mulatta.fusion -t 0.99

#SAQNTI was used for quality-control for transcript isoform
#Installation instructions for SAQNTI can be found at https://github.com/ConesaLab/SQANTI
#Performing the in-depth characterization of transcripts
/home/bin/python-2.7.13/bin/python /home/bin/ConesaLab-sqanti-6927e53e56d2/sqanti_qc.py Macaca_mulatta.all.collapsed.filtered.rep.fasta /ref/gtf/Macaca_mulatta.gtf /ref/genome/Macaca_mulatta.fa -fl Macaca_mulatta.all.collapsed.filtered.abundance.txt -c /data/rnaseq/01.align-genome/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1SJ.out.tab,/data/rnaseq/01.align-genome/Macaca_mulatta/R02027/rep2/Macaca_mulatta.R02027.rep2SJ.out.tab,/data/rnaseq/01.align-genome/Macaca_mulatta/R02027/rep3/Macaca_mulatta.R02027.rep3SJ.out.tab,/data/rnaseq/01.align-genome/Macaca_mulatta/R05040/rep1/Macaca_mulatta.R05040.rep1SJ.out.tab,/data/rnaseq/01.align-genome/Macaca_mulatta/R05040/rep2/Macaca_mulatta.R05040.rep2SJ.out.tab,/data/rnaseq/01.align-genome/Macaca_mulatta/R05040/rep3/Macaca_mulatta.R05040.rep3SJ.out.tab -e /analysis_RNA/02.align-isoform/Macaca_mulatta/transcript_tpms_all_samples.tsv -x /ref/genome/gmapindex/Macaca_mulatta/Macaca_mulatta -t 5
#Applying matching learning methods to filter transcripts that are likely to be artifacts
/home/bin/python-2.7.13/bin/python /home/bin/ConesaLab-sqanti-6927e53e56d2/sqanti_filter.py Macaca_mulatta.all.collapsed.filtered.rep_classification.txt -i Macaca_mulatta.all.collapsed.filtered.rep_corrected.fasta

#SUPPA was used for identifying alternative splicing (AS) event
#Installation instructions for SUPPA can be found at https://github.com/comprna/SUPPA
#Producing gtf containing gene information from SQANTI
perl /home/in-house-script/addGene2gtf.pl Macaca_mulatta.all.collapsed.filtered.rep_classification.txt Macaca_mulatta.all.collapsed.filtered.rep_corrected.gtf > Mmul.isof.gtf
#Identifying AS event
python /home/bin/SUPPA-master/suppa.py generateEvents -i Mmul.isof.gtf -o Mmul.isof.AS.suppa -f ioe -e {SE,SS,MX,RI,FL}
