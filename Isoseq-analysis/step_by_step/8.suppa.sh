#get the suitable gtf format for suppa input
date; perl step2.gtf.pl Macaca_mulatta.all.collapsed.filtered.rep_classification.txt Macaca_mulatta.all.collapsed.filtered.rep_corrected.gtf > Mmul.isof.gtf; date
date; python /home/bin/SUPPA-master/suppa.py generateEvents -i Mmul.isof.gtf -o ./Result_Event/Mmul.isof.AS.suppa -f ioe -e {SE,SS,MX,RI,FL}; date
