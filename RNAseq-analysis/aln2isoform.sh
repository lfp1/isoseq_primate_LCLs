#Aligning pair-end RNA-seq data to transcript isoform for quantification by kallisto
#Installation instructions for kallisto can be found at https://github.com/pachterlab/kallisto
#The same analysis was carried out on five primate species. As an example, we'll look at data from one replicate of the Macaca mulatta

#Building an index from a FASTA formatted file of target sequences
kallisto index -i Macaca_mulatta.kidx Macaca_mulatta.all.collapsed.filtered.rep.fa

#Running the quantification algorithm
kallisto quant -i Macaca_mulatta.kidx --rf-stranded -t 20 -o ./ /data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.1.fq.gz /data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.2.fq.gz
