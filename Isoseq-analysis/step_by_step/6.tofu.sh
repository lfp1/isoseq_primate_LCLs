#retain the uniquely-mapped reads and remove the reads mapped to chrM
date; perl step1.filt_fa.pl TC_clean.sam TC_clean.fa > Macaca_mulatta.correct.trim.rm.fa; perl step1.filt_sam.pl TC_clean.sam > Macaca_mulatta.GMAP.rm.sam; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/collapse_isoforms_by_sam.py --input Macaca_mulatta.correct.trim.rm.fa -s Macaca_mulatta.GMAP.rm.sam -o Macaca_mulatta.all; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/sequence/fa2fq.py Macaca_mulatta.all.collapsed.rep.fa; date
date; mv Macaca_mulatta.all.collapsed.rep.fastq Macaca_mulatta.all.collapsed.rep.fq; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/get_abundance_post_collapse.py Macaca_mulatta.all.collapsed Macaca_mulatta.polished.cluster_report.csv; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/filter_away_subset.py Macaca_mulatta.all.collapsed; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/annotation/make_file_for_subsampling_from_collapsed.py -i Macaca_mulatta.all.collapsed.filtered -o Macaca_mulatta.all.collapse.filtered.for_subsampling.txt; date
date; /home/bin/Python-2.7.15/bin/python2 /home/bin/cDNA_Cupcake-master.20181122/cDNA_Cupcake-master/cupcake/tofu/fusion_finder.py --input Macaca_mulatta.correct.trim.fa -s TC_clean.sam -o Macaca_mulatta.fusion -t 0.99; date
