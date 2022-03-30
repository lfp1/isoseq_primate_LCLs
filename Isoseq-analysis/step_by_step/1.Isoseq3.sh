# Iso-Seq3 (SMRTLINK v6) was used to process all PacBio Sequel subreads, which can be obtained by running the following command:
#wget -P /home/bin/ https://downloads.pacbcloud.com/public/software/installers/smrtlink_6.0.0.47841.zip

#The same analysis was carried out on five primate species. As an example, we'll look at data from the Macaca mulatta.
#First, all subreads from each bio-replicate, library, and tech-replicate were combined using the bamtools merge function, which requires an input for all above subread file directories (refer to 'config.in.lst').
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/bamtools merge -list config.in.lst -out Macaca_mulatta.subreads.bam; date

#Sorting the combined subread bam by samtools sort function.
date; /home/bin/samtools-0.1.19/bin/samtools sort -n -@ 20 Macaca_mulatta.subreads.bam Macaca_mulatta.subreads.sort; rm Macaca_mulatta.subreads.bam; date

#Building pacbio index for the bam file
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/pbindex Macaca_mulatta.subreads.sort.bam; date

#Producing circular consensus sequences (CCS)
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/ccs Macaca_mulatta.subreads.sort.bam Macaca_mulatta.ccs.bam --noPolish --minPasses 1 --numThreads 50 --reportFile Macaca_mulatta.ccs_report.txt; date

#Identifying full-length non-chimeric (FLNC) CCS utilizing demultiplex barcoding algorithm (LIMA) with special `--isoseq` mode with known primer sequence (refer to 'primer.fa')
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/lima --isoseq --dump-clips --no-pbi Macaca_mulatta.ccs.bam primer.fa Macaca_mulatta.demux.bam; date

#Predicting de novo consensus isoforms from classified FLNC CCS using the ICE (Iterative Clustering and Error Correction) algorithm
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/isoseq3 cluster Macaca_mulatta.demux.primer_5p--primer_3p.bam Macaca_mulatta.unpolished.bam -j 50 --require-polya; date

#Polishing predicted consensus isoforms using Quiver and classifying the polished isoforms into high-QV or low-QV isoforms based on predicted accuracy
date; /home/bin/smrtlink_6.0.0/smrtcmds/bin/isoseq3 polish -j 50 Macaca_mulatta.unpolished.bam Macaca_mulatta.subreads.sort.bam Macaca_mulatta.polished.bam; date
