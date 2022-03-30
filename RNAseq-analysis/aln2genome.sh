#Aligning pair-end RNA-seq data to reference genome by STAR
#Installation instructions for STAR can be found at https://github.com/alexdobin/STAR
#The same analysis was carried out on five primate species. As an example, we'll look at data from one replicate of the Macaca mulatta

#Building index for mapping
/home/bin/STAR-2.7.0a/source/STAR --runThreadN 50 --runMode genomeGenerate --limitGenomeGenerateRAM 1000000000000 --outFileNamePrefix Macaca_mulatta --genomeDir /ref/genome/starindex/Macaca_mulatta/ --genomeFastaFiles /ref/genome/Macaca_mulatta.fa --sjdbGTFfile /ref/gtf/Macaca_mulatta.gtf --sjdbOverhang 100

#Running mapping
/home/bin/STAR-2.7.0a/source/STAR --runThreadN 10 --readFilesCommand zcat --genomeDir /ref/genome/starindex/Macaca_mulatta/ --readFilesIn /data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.1.fq.gz /data/rnaseq/Macaca_mulatta/R02027/rep1/Macaca_mulatta.R02027.rep1.2.fq.gz --readStrand Reverse --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Macaca_mulatta.R05040.rep1
