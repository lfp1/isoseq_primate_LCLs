#Selecting representative ORF for each species
#The same analysis was carried out on five primate species. As an example, we'll look at data from the Macaca mulatta.
perl /home/in-house-script/gtf2gff.pl /ref/gtf/Macaca_mulatta.gtf > Macaca_mulatta.all.gff
perl /home/in-house-script/check_orf_for_gff.pl Macaca_mulatta.all.gff /ref/genome/Macaca_mulatta.fa > Macaca_mulatta.all.gff.orf
perl /home/in-house-script/cleanpep.pl /ref/gtf/Macaca_mulatta.gtf /ref/pep/Macaca_mulatta.pep.fa > Macaca_mulatta.pep
perl /home/in-house-script/bestProteinFromEnsembl.pl Macaca_mulatta.pep Macaca_mulatta.all.gff.orf 1 > Macaca_mulatta.pep.best
perl /home/in-house-script/select_gff.pl Macaca_mulatta.pep.best Macaca_mulatta.all.gff > Macaca_mulatta.best.gff
perl /home/in-house-script/getGene.pl Macaca_mulatta.best.gff /ref/genome/Macaca_mulatta.fa > Macaca_mulatta.best.cds.fa
perl /home/in-house-script/check_orf_for_cds.pl Macaca_mulatta.best.cds.fa | awk '$4>0' > Macaca_mulatta.best.cds.prestop
perl /home/in-house-script/cds2aa.pl Macaca_mulatta.best.cds.fa > Macaca_mulatta.best.pep.fa
perl /home/in-house-script/select_orf.pl Macaca_mulatta.best.pep.fa Macaca_mulatta.all.gff.orf > Macaca_mulatta.best.gff.orf
awk '$8==1 && $9==1' Macaca_mulatta.best.gff.orf > Macaca_mulatta.best.gff.orf.intact
sed "s/>/>Mmul/" Macaca_mulatta.best.pep.fa > Mmul.best.pep.fa
sed "s/=/=Mmul/" Macaca_mulatta.best.gff > Mmul.best.pep.gff

#Combining all representative protein sequence in five species and running the mapping with BLASTP
cat Mmul.best.pep.fa Pabe.best.pep.fa Ggor.best.pep.fa Ptro.best.pep.fa Hsap.best.pep.fa > 5species.best.pep.fa
#Formatting database sequences for mapping
/home/bin/formatdb -i 5species.best.pep.fa -p T -o T
#Spliting fasta and running blastp independtly
mkdir split; cd split; perl /home/in-house-script/Split_fasta.pl -seq ../5species.best.pep.fa -nf 200 -od ./; perl /home/in-house-script/call_blast.pl ./ ../5species.best.pep.fa blastp ./ best.pep.fa 0.01 queue_ID project_ID; cd ..;
cat split/5species_*.best.pep.fa.m8 > 5species.best.pep.fa.m8; rm -r split

#Obtaining pair-wise blastp result
perl /home/in-house-script/select_m8.pl 5species.best.pep.fa.m8 config_for_ortholog ./
for i in *.m8; do perl /home/bin/solar-0.9.6/solar.pl -a prot2prot -f m8 $i > $i.solar; perl /home/in-house-script/solar_add_realLen.pl $i.solar 5species.best.pep.fa > $i.solar.cor; perl /home/in-house-script/solar_add_identity.pl --solar $i.solar.cor --m8 $i > $i.solar.cor.idAdd; done

#Defining orthologous genes among five species by using human as reference
perl /home/in-house-script/ort.pl Hsap.best.pep.gff Mmul.best.pep.gff Hsap_Mmul.m8.solar.cor.idAdd > Hsap_Mmul.m8.solar.cor.idAdd.ort
perl /home/in-house-script/ort.pl Hsap.best.pep.gff Pabe.best.pep.gff Hsap_Pabe.m8.solar.cor.idAdd > Hsap_Pabe.m8.solar.cor.idAdd.ort
perl /home/in-house-script/ort.pl Hsap.best.pep.gff Ggor.best.pep.gff Hsap_Ggor.m8.solar.cor.idAdd > Hsap_Ggor.m8.solar.cor.idAdd.ort
perl /home/in-house-script/ort.pl Hsap.best.pep.gff Ptro.best.pep.gff Hsap_Ptro.m8.solar.cor.idAdd > Hsap_Ptro.m8.solar.cor.idAdd.ort
perl /home/in-house-script/link_ort.pl Hsap_Mmul.m8.solar.cor.idAdd.ort Hsap_Pabe.m8.solar.cor.idAdd.ort Hsap_Ggor.m8.solar.cor.idAdd.ort Hsap_Ptro.m8.solar.cor.idAdd.ort config_for_ortholog > 5species.ort
