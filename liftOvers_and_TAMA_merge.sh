#!/bin/bash

# We restrict to real 'Isoform' (not 'Artifacts' according to SQANTI filtering) located in orthologous 1:1 genes in the 5 species.
# Get their GTFs to be lifted afterwards (BED12)

for species in Homo_sapiens Pan_troglodytes Gorilla_gorilla Pongo_abelii Macaca_mulatta
do
	sqanti_decision="/media/luis/Data/PhD/Isoseq_LCLs/frozen_filtered_isoseq_and_epigenomics/${species}.all.collapsed.filtered.rep_classification.txt_filterResults.txt"

	gtf_file="/media/luis/Data/PhD/Isoseq_LCLs/Result_QC.sqanti/${species}.all.collapsed.filtered.rep_corrected.gtf.gz"

	ort_genes="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/5species.ort.overlap.filt"

	# grep -w will not take any antisense or fusion
	# Get GTF for those isoforms
	grep -w -Ff <(awk -F'\t' '{print $1}' 5species.ort.overlap.filt | awk -F'_' '{print $3}') <(awk -F'\t' '$38=="Isoform"' ${sqanti_decision}) | awk -F'\t'  '{print $1}' | grep -w -Ff - <(zcat ${gtf_file}) > ${species}_isoform_ort_one2one.gtf

	# GTF to BED12
	tools/gtfToGenePred ${species}_isoform_ort_one2one.gtf ${species}_isoform_ort_one2one.genepred
	tools/genePredToBed ${species}_isoform_ort_one2one.genepred ${species}_isoform_ort_one2one.bed12

	# Add species to BED12	
	awk -F'\t' -v OFS="\t" -v species=$species '$4=species"_"$4' ${species}_isoform_ort_one2one.bed12 > tmp && mv tmp ${species}_isoform_ort_one2one.bed12

done

# Chain files from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenpath/panTro5/liftOver/panTro5ToPanTro6.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/gorGor4/liftOver/gorGor4ToGorGor6.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/ponAbe2/liftOver/ponAbe2ToPonAbe3.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/rheMac8/liftOver/rheMac8ToRheMac10.over.chain.gz

wget http://hgdownload.cse.ucsc.edu/goldenpath/panTro6/liftOver/panTro6ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/gorGor6/liftOver/gorGor6ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/ponAbe3/liftOver/ponAbe3ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/rheMac10/liftOver/rheMac10ToHg38.over.chain.gz

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToPanTro6.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToGorGor6.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToPonAbe3.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToRheMac10.over.chain.gz

# LiftOvers (upgrade to new assemblies)
liftOver Pan_troglodytes_isoform_ort_one2one.bed12 panTro5ToPanTro6.over.chain.gz panTro5ToPanTro6.bed12 panTro5ToPanTro6.unmapped
liftOver Gorilla_gorilla_isoform_ort_one2one.bed12 gorGor4ToGorGor6.over.chain.gz gorGor4ToGorGor6.bed12 gorGor4ToGorGor6.unmapped
liftOver Pongo_abelii_isoform_ort_one2one.bed12 ponAbe2ToPonAbe3.over.chain.gz ponAbe2ToPonAbe3.bed12 ponAbe2ToPonAbe3.unmapped
liftOver Macaca_mulatta_isoform_ort_one2one.bed12 rheMac8ToRheMac10.over.chain.gz rheMac8ToRheMac10.bed12 rheMac8ToRheMac10.unmapped

# LiftOvers from NHP to human (less conservative minMatch for cross-sp liftOver)
liftOver -minMatch=0.5 panTro5ToPanTro6.bed12 panTro6ToHg38.over.chain.gz panTro6ToHg38.bed12 panTro6ToHg38.unmapped
liftOver -minMatch=0.5 gorGor4ToGorGor6.bed12 gorGor6ToHg38.over.chain.gz gorGor6ToHg38.bed12 gorGor6ToHg38.unmapped
liftOver -minMatch=0.5 ponAbe2ToPonAbe3.bed12 ponAbe3ToHg38.over.chain.gz ponAbe3ToHg38.bed12 ponAbe3ToHg38.unmapped
liftOver -minMatch=0.5 rheMac8ToRheMac10.bed12 rheMac10ToHg38.over.chain.gz rheMac10ToHg38.bed12 rheMac10ToHg38.unmapped

# Collapse isoforms and create merged models (TAMA) in hg38 genome
# Modify BED12 for TAMA merge
arr=`find -name "*ToHg38.bed12" -o -name "Homo_sapiens_isoform_ort_one2one.bed12"`
for file in $arr; do awk -F'\t' -v OFS="\t" '$4=$4";"$4' $file > $file.TAMA && mv $file.TAMA TAMA/ ; done

# Create 'filenames_TAMA.txt'
find "/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA" -name "*TAMA" -type f | awk -F'\t' -v OFS="\t" '{print $1,"no_cap","1,1,1","PB"}' > TAMA/filenames_TAMA.txt

# TAMA merge (collapse exon ends as well)
source activate python2
python /home/luis/tama/tama_merge.py -f TAMA/filenames_TAMA.txt -p TAMA/merged_models_in_human -a 100 -m 100 -z 100 -e longest_ends

###############################################################################################
# The merged model in human will be projected to NHP. Keep only isoforms lifted in the 5 species
###############################################################################################

merged_human="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/merged_models_in_human.bed"
outPath="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models"

liftOver -minMatch=0.5 ${merged_human} hg38ToPanTro6.over.chain.gz ${outPath}/merged_models_in_chimpanzee.bed ${outPath}/merged_models_in_chimpanzee.unmapped
liftOver -minMatch=0.5 ${merged_human} hg38ToGorGor6.over.chain.gz ${outPath}/merged_models_in_gorilla.bed ${outPath}/merged_models_in_gorilla.unmapped
liftOver -minMatch=0.5 ${merged_human} hg38ToPonAbe3.over.chain.gz ${outPath}/merged_models_in_orangutan.bed ${outPath}/merged_models_in_orangutan.unmapped
liftOver -minMatch=0.5 ${merged_human} hg38ToRheMac10.over.chain.gz ${outPath}/merged_models_in_macaque.bed ${outPath}/merged_models_in_macaque.unmapped

# Select isoform merged models that are present in the 5 species
find "/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA" -name "*bed" -exec cat {} \; | awk -F'\t' '{print $4}' | sort | uniq -c | sed -r 's/^\s+//g' | awk '$1=="5" {print $2}' > ${outPath}/ids_present_in_5_species.txt

# Filter BED12 to 49.797 transcript models in each species (comparative)
grep -w -Ff ${outPath}/ids_present_in_5_species.txt ${merged_human} > ${outPath}/present_in_5_species/Homo_sapiens_5_sp.bed12
grep -w -Ff ${outPath}/ids_present_in_5_species.txt ${outPath}/merged_models_in_chimpanzee.bed > ${outPath}/present_in_5_species/Pan_troglodytes_5_sp.bed12
grep -w -Ff ${outPath}/ids_present_in_5_species.txt ${outPath}/merged_models_in_gorilla.bed > ${outPath}/present_in_5_species/Gorilla_gorilla_5_sp.bed12
grep -w -Ff ${outPath}/ids_present_in_5_species.txt ${outPath}/merged_models_in_orangutan.bed > ${outPath}/present_in_5_species/Pongo_abelii_5_sp.bed12
grep -w -Ff ${outPath}/ids_present_in_5_species.txt ${outPath}/merged_models_in_macaque.bed > ${outPath}/present_in_5_species/Macaca_mulatta_5_sp.bed12

############################################################################
# Get FASTA files of transcript models for Kallisto pseudomapping (Methods)
############################################################################

assemblies="/media/luis/Data/PhD/Isoseq_LCLs/new_assemblies_fasta"

/bin/bedtools2/bin/bedtools getfasta -name -s -split -fi ${assemblies}/hg38.fa -bed ${outPath}/present_in_5_species/Homo_sapiens_5_sp.bed12 > ${outPath}/present_in_5_species/Homo_sapiens_5_sp.fasta

/bin/bedtools2/bin/bedtools getfasta -name -s -split -fi ${assemblies}/panTro6.fa -bed ${outPath}/present_in_5_species/Pan_troglodytes_5_sp.bed12 > ${outPath}/present_in_5_species/Pan_troglodytes_5_sp.fasta

/bin/bedtools2/bin/bedtools getfasta -name -s -split -fi ${assemblies}/gorGor6.fa -bed ${outPath}/present_in_5_species/Gorilla_gorilla_5_sp.bed12 > ${outPath}/present_in_5_species/Gorilla_gorilla_5_sp.fasta

/bin/bedtools2/bin/bedtools getfasta -name -s -split -fi ${assemblies}/ponAbe3.fa -bed ${outPath}/present_in_5_species/Pongo_abelii_5_sp.bed12 > ${outPath}/present_in_5_species/Pongo_abelii_5_sp.fasta

/bin/bedtools2/bin/bedtools getfasta -name -s -split -fi ${assemblies}/rheMac10.fa -bed ${outPath}/present_in_5_species/Macaca_mulatta_5_sp.bed12 > ${outPath}/present_in_5_species/Macaca_mulatta_5_sp.fasta


