date; /home/bin/gmap-20170215/bin/gmap -D /ref/genome/gmapindex -d Macaca_mulatta -f samse -t 50 -n 0 Macaca_mulatta.correct.trim.fa > Macaca_mulatta.GMAP.tmp.sam; date
date; perl /home/bin/RNA_PacBio_ISOseq_Reference_2017a/RNA_module/bin/split_sam.pl Macaca_mulatta.GMAP.tmp.sam ./; date
date; cat header | sort | uniq | sort -k 3,3 -k 4,4n > header2; sort -k 3,3 -k 4,4n sam > sam2; cat header2 sam2 > Macaca_mulatta.GMAP.sam; rm Macaca_mulatta.GMAP.tmp.sam header header2 sam sam2; date
date; /home/bin/samtools-0.1.19/bin/samtools view -bhS Macaca_mulatta.GMAP.sam -o Macaca_mulatta.GMAP.bam; date
