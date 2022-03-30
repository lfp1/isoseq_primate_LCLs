date
conda activate transcript_clean
python /home/bin/TranscriptClean-master/TranscriptClean.py --sam Macaca_mulatta.GMAP.sam --genome /ref/genome/Macaca_mulatta.fa --spliceJns /analysis_RNA/01.align-genome/Macaca_mulatta/all.SJ.out.tab --outprefix ./
conda deactivate
date
