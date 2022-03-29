##################################################
# Main code for differential exon usage analyses #
##################################################

# List of exon count files resulting from htseq (order in NHP has been set in reference to the order in the human reference)
countFiles<-c(list.files(path="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results",pattern="Homo_sapiens.+exonic_counts.txt",full.names=TRUE),list.files(path="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order",pattern="exonic_counts_final.txt",full.names=TRUE))

# Create the design dataframe from the previous vector. Include experimental batches
sampleTable = data.frame(
row.names = gsub(".+?\\.(.+)_exonic_counts[_final]*.txt","\\1",basename(countFiles)),
condition = gsub("(.+?)\\..+","\\1",basename(countFiles)),
libType= rep("paired_end",16),
batch= c("snyder","snyder","raquel","luis","raquel","luis","raquel","raquel","raquel","luis","raquel","luis","raquel","luis","luis","raquel"))

# Fix the '-' in CRL-1850 to avoid warning message
rownames(sampleTable)[rownames(sampleTable)=="CRL-1850.rep1"]<-"CRL1850.rep1"

# While all RNA-seq has been loaded, we removed GM12878.rep1 as in the rest of analyses (technical replicate) #
countFiles<-countFiles[-1]
sampleTable<-sampleTable[-1,]

# Iterate over pairwise comparatives while subsetting the above dfs
list_species<-c("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta")
pairwise_comparisons<-as.data.frame(t(combn(list_species, 2)))
colnames(pairwise_comparisons)<-c("species1","species2")

# Run DEXSeq pairwise comparisons for differential exon usage (10) #
library("DEXSeq")
BPPARAM = MulticoreParam(workers=4)

for (i in 1:nrow(pairwise_comparisons)) {
	
	species1<-as.character(pairwise_comparisons[i,"species1"])
	species2<-as.character(pairwise_comparisons[i,"species2"])
	toMatch <- c(species1,species2)

	# sampleTable and countFiles are subsetted and order of samples is kept in the grepping
	sampleTable_subset<-sampleTable[grep(paste(toMatch,collapse="|"),sampleTable$condition),]
	sampleTable_subset<-droplevels(sampleTable_subset) # drop unused levels to prevent DEXSeq warnings
	countFiles_subset<-countFiles[grep(paste(toMatch,collapse="|"),countFiles)]
	
	# Read the counts and design for each pairwise comparative (subset).
	print(paste("Reading files for",paste(species1,species2,sep="-")))

	dxd = DEXSeqDataSetFromHTSeq(
	countFiles_subset,
	sampleData=sampleTable_subset,
	design= ~ sample + exon + batch:exon + condition:exon,
	flattenedfile="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/Homo_sapiens_5_sp_ens_exonic_parts_filtered.gff" )

	# Estimate size factors
	dxd = estimateSizeFactors( dxd )

	# Include batch effect in estimateDispersions
	print(paste("Estimating dispersions for",paste(species1,species2,sep="-")))

	formulaFullModel    =  ~ sample + exon + batch:exon + condition:exon
	formulaReducedModel =  ~ sample + exon + batch:exon 
	dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=BPPARAM )

	# DEU
	print(paste("Testing DEU for",paste(species1,species2,sep="-")))

	dxd = testForDEU ( dxd, 
	reducedModel = formulaReducedModel, 
	fullModel = formulaFullModel ,
	BPPARAM=BPPARAM )

	# Estimate fold-changes
	print(paste("Estimating fold changes for",paste(species1,species2,sep="-")))
	dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM )

	# Generate final DEXSeq results
	dxr1 = DEXSeqResults( dxd )
	
	# Create two objects per comparative
	assign(x=paste(species1,species2,"dxd",sep="_"),value=dxd)
	assign(x=paste(species1,species2,"dxr1",sep="_"),value=dxr1)

	print(paste("Finished ",paste(species1,species2,sep="-")))
}


# Also, DEXSeq run with all samples together #

# Read the counts and design.
print("Reading files")
dxd = DEXSeqDataSetFromHTSeq(
countFiles,
sampleData=sampleTable,
design= ~ sample + exon + batch:exon + condition:exon,
flattenedfile="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/Homo_sapiens_5_sp_ens_exonic_parts_filtered.gff" )
# Estimate size factors
dxd = estimateSizeFactors( dxd )
# Include batch effect in estimateDispersions
print("Estimating dispersions")
formulaFullModel    =  ~ sample + exon + batch:exon + condition:exon
formulaReducedModel =  ~ sample + exon + batch:exon 
dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=BPPARAM )
# DEU
print("Testing DEU")
dxd = testForDEU ( dxd, 
reducedModel = formulaReducedModel, 
fullModel = formulaFullModel ,
BPPARAM=BPPARAM )
# Estimate fold-changes
print("Estimating fold changes")
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM )
# Generate final DEXSeq results
dxr1 = DEXSeqResults( dxd )
print("Finished")

###########################################
# Retrieve significant exon usage changes #
###########################################

library("dplyr")
library("ggplot2")
library("ggsci")

# Read file with the average phastCons score for each exonic part
phastCons_exonic_parts<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/phastCons_Homo_sapiens_exonic_parts_filtered.txt",sep="\t",header=FALSE)
phastCons_exonic_parts<-phastCons_exonic_parts %>% dplyr::select(V1,V6)
colnames(phastCons_exonic_parts)<-c("exonic_part","mean_phastCon")

# Get DEXSeq result objects
DEXSeq_result_objects<-ls(pattern="*_dxr1") 
# Put DEXSeq objects into into a list
DEXSeq_result_list<-do.call("list",mget(DEXSeq_result_objects)) 
# Convert them to dataframe
DEXSeq_result_list <- lapply(DEXSeq_result_list, function(x) { x<-as.data.frame(x) } ) 

# Explore exonic parts that are significant in each pairwise comparative through upset plot
library("UpSetR")
# Significant exonic parts: exonBaseMean>10 & padj<0.05 & abs(.[[10]])>1.2
significant_exonic_parts_10_cmp<-lapply(DEXSeq_result_list, function(x) { x %>% dplyr::filter(exonBaseMean>10 & padj<0.05 & abs(.[[10]])>1.2) %>% rownames(.) } )
upset(fromList(significant_exonic_parts_10_cmp),nsets = 10,order.by="freq")

# Define 2 functions to easily retrieve UpSetR intersections (from @docmanny in https://github.com/hms-dbmi/UpSetR/issues/85)
# This function will tell you the original rownames (default fromList() in UpSetR is not returning the names)
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
      }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
  }

# To get desired intersection subsets
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

# Get species-specific exon usage changes
hs<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1',
'Homo_sapiens_Gorilla_gorilla_dxr1',
'Homo_sapiens_Pongo_abelii_dxr1',
'Homo_sapiens_Macaca_mulatta_dxr1') 

cs<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1',
'Pan_troglodytes_Gorilla_gorilla_dxr1',
'Pan_troglodytes_Pongo_abelii_dxr1',
'Pan_troglodytes_Macaca_mulatta_dxr1') 

gs<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Gorilla_gorilla_dxr1',
'Pan_troglodytes_Gorilla_gorilla_dxr1',
'Gorilla_gorilla_Pongo_abelii_dxr1',
'Gorilla_gorilla_Macaca_mulatta_dxr1') 

os<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pongo_abelii_dxr1',
'Pan_troglodytes_Pongo_abelii_dxr1',
'Gorilla_gorilla_Pongo_abelii_dxr1',
'Pongo_abelii_Macaca_mulatta_dxr1') 

ms<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Macaca_mulatta_dxr1',
'Pan_troglodytes_Macaca_mulatta_dxr1',
'Gorilla_gorilla_Macaca_mulatta_dxr1',
'Pongo_abelii_Macaca_mulatta_dxr1') 

# Concatenate all exonic parts showing species-specific usage changes together
all_sp_specific_exonic_parts<-rbind(hs,cs,gs,os,ms)
# For species-specific exonic parts, represent the usage coefficient from dxr1
library("pheatmap")
dxr1 %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% rownames(all_sp_specific_exonic_parts)) %>% dplyr::select("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta") %>% pheatmap(.,scale = "row",show_rownames = FALSE)

###########################################################################
### Classify exonic parts and load annotation features for each of them ###
###########################################################################

# Convert the global analyses to data frame
dxr1_df<-as.data.frame(dxr1)
# Remove genes with only one exonic part! 147 genes were removed in this step. Keep rownames
dxr1_df<-dxr1_df %>% rownames_to_column('my_rows') %>% group_by(groupID) %>% filter(n()>1) %>% ungroup() %>% column_to_rownames('my_rows')

# Add phastCons to the dataframe of exonic parts
dxr1_df<-merge(dxr1_df,phastCons_exonic_parts,by.x="row.names",by.y="exonic_part")

# Classify exonic parts by the changes they display
# The exonic parts that are not changing are those NOT in fromList(significant_exonic_parts_10_cmp) (no significant change in any of the 10 comparatives)
dxr1_df<-dxr1_df %>% mutate(category = case_when(Row.names %in% rownames(hs) ~ "HS",Row.names %in% rownames(cs) ~ "CS",Row.names %in% rownames(gs) ~ "GS",Row.names %in% rownames(os) ~ "OS",Row.names %in% rownames(ms) ~ "MS",!dxr1_df$Row.names %in% rownames(fromList(significant_exonic_parts_10_cmp))~"Unchanged",TRUE ~ "Other_changing"))

# Order of levels
dxr1_df$category <- factor(dxr1_df$category, levels = c("Unchanged","HS","CS","GS","OS","MS","Other_changing"))

#########################################
# CDS intersection (Linux command-line) #
#########################################

# gff="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/SQANTI/Homo_sapiens_5_sp_SQANTI3_corrected.gtf.cds.gff"
# awk -F'\t' '$3=="CDS"' $gff | awk -F'\t' -v OFS="\t" '{print $1,$4-1,$5,"CDS",".",$7}' > hg38_cds.bed
# exonic_parts="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/Homo_sapiens_5_sp_ens_exonic_parts_filtered.gff"
# grep -w "exonic_part" ${exonic_parts} | awk -v OFS="\t" '{print $1,$4-1,$5,$10,$14,$6,$7}' | sed 's/"//g' | sed 's/;//g' | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4":E"$5,$6,$7}' > Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed
# Annotate number of CDS intersections in each exonic part
# /bin/bedtools2/bin/bedtools intersect -c -s -a Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed -b hg38_cds.bed > CDS_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed

##################################################
# Repeat masker intersection (Linux command-line) #
##################################################

# Download repeatmasker BED file from UCSC tables (hg38)
# Not strand-specific intersection
# /bin/bedtools2/bin/bedtools intersect -c -a Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed -b hg38_repeatmasker.bed > repeatmasker_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed

##############################################
# Seg dups intersection (Linux command-line) #
##############################################

# Download segmental duplications BED file from UCSC tables (hg38)
# Not strand-specific intersection
# /bin/bedtools2/bin/bedtools intersect -c -a Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed -b hg38_segdups.bed > segdups_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed

################
# Annotate CDS #
################

CDS<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/CDS_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed",sep="\t",header=FALSE)

# Define if there is CDS overlap or not
CDS$cds<-ifelse(CDS$V7>0,"CDS","no_CDS")
CDS<-CDS %>% dplyr::select(V4,cds)
colnames(CDS)<-c("exonic_part","cds")

# Add CDS intersection to exonic parts
dxr1_df<-merge(dxr1_df,CDS,by.x="Row.names",by.y="exonic_part")

############################
# Annotate repeat elements #
############################

repeat_masker<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/repeatmasker_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed",sep="\t",header=FALSE)

# Define if there is repeat elements overlap or not
repeat_masker$repeats<-ifelse(repeat_masker$V7>0,"Repeats","No_repeats")
repeat_masker<-repeat_masker %>% dplyr::select(V4,repeats)
colnames(repeat_masker)<-c("exonic_part","repeats")

# Add overlap of repeat elements with exonic parts
dxr1_df<-merge(dxr1_df,repeat_masker,by.x="Row.names",by.y="exonic_part")

###################################
# Annotate segmental duplications #
###################################

segdups<-read.table("/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/segdups_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed",sep="\t",header=FALSE)

# Define if there is overlap with segmental duplications
segdups$segdup<-ifelse(segdups$V7>0,"Segdup","No_segdup")
segdups<-segdups %>% dplyr::select(V4,segdup)
colnames(segdups)<-c("exonic_part","segdup")

# Add overlap with segmental duplications with exonic parts
dxr1_df<-merge(dxr1_df,segdups,by.x="Row.names",by.y="exonic_part")

####################################################################
#### APPRIS principal isoforms intersection (Linux command-line) ###
####################################################################

# exonic_parts="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed"
# hg38_ensembl="/media/luis/Data/PhD/Isoseq_LCLs/new_assemblies_fasta/annotations/hg38.ensembl.gtf"
# From http://appris-tools.org/#/downloads, download Gencode27/Ensembl90 (APPRIS scores txt). Note it is not identical to V91 but very close.
# appris="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/APPRIS/appris_data.appris.txt"

# For principal isoforms, get transcript coordinates and intersect with our exonic parts
# /bin/bedtools2/bin/bedtools intersect -s -c -a $exonic_parts -b <(grep -w "PRINCIPAL" $appris | awk -F'\t' '{print $3}' | grep -w -Ff - $hg38_ensembl | awk -F'\t' '$3=="exon"') > appris_intersection.txt
# Do the same for the alternative/minor/other isoforms. 
# /bin/bedtools2/bin/bedtools intersect -s -c -a appris_intersection.txt -b <(grep -v -w "PRINCIPAL" $appris | awk -F'\t' '{print $3}' | grep -w -Ff - $hg38_ensembl | awk -F'\t' '$3=="exon"') > tmp && mv tmp appris_intersection.txt
# 'appris_intersection.txt': 7th column is the overlapping with 'principal' exons, 8th column is the overlapping with 'minor' exons.
# Some will be overlapping both 'principal' and 'minor'.

###########################
# Annotate APPRIS overlap #
###########################

# Read file containing the overlaps with exons that are part of 'principal' or 'minor' isoforms in APPRIS
appris<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/APPRIS/appris_intersection.txt",sep="\t",header=FALSE)
# Classify each exonic part
appris<-appris %>% mutate(appris_overlap = case_when(V7>0 & V8>0 ~ "Overlap_principal_and_minor",V7>0 & V8==0 ~ "Only_principal",V7==0 & V8>0 ~ "Only_minor", V7==0 & V8==0 ~ "Unclassified"))
appris<-appris %>% dplyr::select(V4,appris_overlap)
colnames(appris)<-c("exonic_part","appris_overlap")
# Add APPRIS information to exonic parts
dxr1_df<-merge(dxr1_df,appris,by.x="Row.names",by.y="exonic_part")

# Set order of levels
dxr1_df$appris_overlap <- factor(dxr1_df$appris_overlap, levels = c("Overlap_principal_and_minor","Only_principal","Only_minor", "Unclassified"))

########################################################
### Get genes showing species-specific DEU changes #####
########################################################

hs_genes<-gsub(":.+","",rownames(hs)) %>% unique()
cs_genes<-gsub(":.+","",rownames(cs)) %>% unique()
gs_genes<-gsub(":.+","",rownames(gs)) %>% unique()
os_genes<-gsub(":.+","",rownames(os)) %>% unique()
ms_genes<-gsub(":.+","",rownames(ms)) %>% unique()


##############################################################
# Load the Percent-Spliced-In (PSI) per each exonic part #####
##############################################################

library(tidyverse)
# Filenames with matched PSI
PSI_filenames<-c(list.files("/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/exclusion_counts/compute_PSI",pattern = "Homo_sapiens.+psi",full.names = TRUE), list.files("/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/exclusion_counts/compute_PSI/matching_order",pattern = "*psi.final.txt",full.names = TRUE) )
# Remove GM12878.rep1
PSI_filenames<-PSI_filenames[-1]
# Load into list of DFs
PSI_dfs<-lapply(PSI_filenames,read.table,sep="\t")
# Second column will be the sample name
for (i in 1:15) {
	colnames(PSI_dfs[[i]])[2] <- gsub(".final.txt","",basename(PSI_filenames))[i]
}
# Merge the list of PSI DFs into a single DF
PSI<-PSI_dfs %>% reduce(left_join, by = "V1")
# Exonic parts to row.names
rownames(PSI)<-PSI$V1
PSI<-PSI %>% dplyr::select(-V1)
# Set NA (0 inclusion and 0 exclusion reads) as 0
PSI[is.na(PSI)]<-0
# Get mean in all samples
PSI$meanPSI<-rowMeans(PSI)
rownames(PSI)<-gsub(":",":E",rownames(PSI))
# Classify into PSI bins
PSI<-PSI %>% mutate(PSI_class = case_when(meanPSI>=0.8 ~ "Highly included", meanPSI<0.8 & meanPSI>=0.4 ~ "Mid included",meanPSI<0.4 & meanPSI>=0.2   ~ "Mid excluded",meanPSI<0.2 ~ "Highly excluded"))
# Add PSI to exonic parts
PSI_classes<-PSI %>% dplyr::select(PSI_class)
dxr1_df<-merge(dxr1_df,PSI_classes,by.x="Row.names",by.y="row.names")
dxr1_df$PSI_class<-as.factor(dxr1_df$PSI_class)
dxr1_df$PSI_class <- factor(dxr1_df$PSI_class, levels = c("Highly included", "Mid included", "Mid excluded", "Highly excluded"))

########################################
### Add AS classification (SUPPA) ######
########################################

# Add if exonic part is SE, RI, A5, A3, AF, AL or MX.
SE_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_SE.bed",sep="\t",header=FALSE)
SE_df$SE<-ifelse(SE_df$V7>0,"Yes","No")
SE_df<-SE_df %>% dplyr::select(V4,SE)
colnames(SE_df)<-c("exonic_part","SE")
dxr1_df<-merge(dxr1_df,SE_df,by.x="Row.names",by.y="exonic_part")

RI_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_RI.bed",sep="\t",header=FALSE)
RI_df$RI<-ifelse(RI_df$V7>0,"Yes","No")
RI_df<-RI_df %>% dplyr::select(V4,RI)
colnames(RI_df)<-c("exonic_part","RI")
dxr1_df<-merge(dxr1_df,RI_df,by.x="Row.names",by.y="exonic_part")

A5_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_A5.bed",sep="\t",header=FALSE)
A5_df$A5<-ifelse(A5_df$V7>0,"Yes","No")
A5_df<-A5_df %>% dplyr::select(V4,A5)
colnames(A5_df)<-c("exonic_part","A5")
dxr1_df<-merge(dxr1_df,A5_df,by.x="Row.names",by.y="exonic_part")

A3_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_A3.bed",sep="\t",header=FALSE)
A3_df$A3<-ifelse(A3_df$V7>0,"Yes","No")
A3_df<-A3_df %>% dplyr::select(V4,A3)
colnames(A3_df)<-c("exonic_part","A3")
dxr1_df<-merge(dxr1_df,A3_df,by.x="Row.names",by.y="exonic_part")

AF_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_AF.bed",sep="\t",header=FALSE)
AF_df$AF<-ifelse(AF_df$V7>0,"Yes","No")
AF_df<-AF_df %>% dplyr::select(V4,AF)
colnames(AF_df)<-c("exonic_part","AF")
dxr1_df<-merge(dxr1_df,AF_df,by.x="Row.names",by.y="exonic_part")

AL_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_AL.bed",sep="\t",header=FALSE)
AL_df$AL<-ifelse(AL_df$V7>0,"Yes","No")
AL_df<-AL_df %>% dplyr::select(V4,AL)
colnames(AL_df)<-c("exonic_part","AL")
dxr1_df<-merge(dxr1_df,AL_df,by.x="Row.names",by.y="exonic_part")

MX_df<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/SUPPA2_all_splicing_events/intersect_with_exonic_parts/EP_vs_MX.bed",sep="\t",header=FALSE)
MX_df$MX<-ifelse(MX_df$V7>0,"Yes","No")
MX_df<-MX_df %>% dplyr::select(V4,MX)
colnames(MX_df)<-c("exonic_part","MX")
dxr1_df<-merge(dxr1_df,MX_df,by.x="Row.names",by.y="exonic_part")

# Alternative splicing classes summarized in a single column:
dxr1_df<-dxr1_df %>% mutate(AS_class = case_when(SE=="Yes" & RI=="No" & A5=="No" & A3=="No" & AF=="No" & AL=="No" & MX=="No" ~ "SE",
SE=="No" & RI=="Yes" & A5=="No" & A3=="No" & AF=="No" & AL=="No" & MX=="No" ~ "RI",
SE=="No" & RI=="No" & A5=="Yes" & A3=="No" & AF=="No" & AL=="No" & MX=="No" ~ "A5",
SE=="No" & RI=="No" & A5=="No" & A3=="Yes" & AF=="No" & AL=="No" & MX=="No" ~ "A3",
SE=="No" & RI=="No" & A5=="No" & A3=="No" & AF=="Yes" & AL=="No" & MX=="No" ~ "AF",
SE=="No" & RI=="No" & A5=="No" & A3=="No" & AF=="No" & AL=="Yes" & MX=="No" ~ "AL",
SE=="No" & RI=="No" & A5=="No" & A3=="No" & AF=="No" & AL=="No" & MX=="Yes" ~ "MX",
SE=="No" & RI=="No" & A5=="No" & A3=="No" & AF=="No" & AL=="No" & MX=="No" ~ "None",
TRUE ~ "Complex_AS"))

# Set order of levels
dxr1_df$AS_class <- factor(dxr1_df$AS_class, levels = c("None","SE","RI","A5","A3","AF","AL","MX","Complex_AS")) 

#################################
## Annotate UCSC Pfam domains ###
#################################

# Linux command line:
# We will add the overlap with Pfam domains from UCSC (download BED12 from UCSC). Note bedtools intersect must use the --split option for BED12
# bedtools intersect -split -c -s -a Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed -b hg38_Pfam_UCSC.bed > Pfam_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed

Pfam_in_exonic_parts<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/Pfam_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed",sep="\t",header=FALSE)
Pfam_in_exonic_parts<-Pfam_in_exonic_parts %>% dplyr::select(V4,V7)
colnames(Pfam_in_exonic_parts)<-c("Row.names","Pfam_UCSC")
Pfam_in_exonic_parts$Pfam_UCSC<-ifelse(Pfam_in_exonic_parts$Pfam_UCSC>0,"Yes","No")

# Add Pfam domains overlap
dxr1_df<-merge(dxr1_df,Pfam_in_exonic_parts)

##################################################################
# Heatmap of species-specific DEU change (only skipping exons) ###
##################################################################

library("ggpubr")
# Select species-specific SE changes.
# Note 'direction' of change is only used for sorting the df
spsp_DEU_SE<-dxr1_df %>% dplyr::filter((category=="HS" | category=="CS" | category=="GS" | category=="OS" | category=="MS") & AS_class=="SE")
spsp_DEU_SE<-spsp_DEU_SE %>% mutate(direction = case_when(category=="HS" & Homo_sapiens>Pan_troglodytes ~ "up",category=="CS" & Pan_troglodytes>Homo_sapiens ~ "up", category=="GS" & Gorilla_gorilla>Homo_sapiens ~ "up",category=="OS" & Pongo_abelii>Homo_sapiens ~ "up",category=="MS" & Macaca_mulatta>Homo_sapiens ~ "up", TRUE ~ "down")) # if we ask up or down in all pairwise comp., the numbers are the same

# Sort the exons for the representation
spsp_DEU_SE_plot <- spsp_DEU_SE %>% arrange(category,direction) %>% dplyr::select("Row.names","Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta") 
colnames(spsp_DEU_SE_plot)<-c("Row.names","Human","Chimpanzee","Gorilla","Orangutan","Rhesus macaque")
rownames(spsp_DEU_SE_plot)<-spsp_DEU_SE_plot$Row.names
spsp_DEU_SE_plot<-spsp_DEU_SE_plot %>% dplyr::select(-Row.names)

# Generate an equivalent df for annotation of exons
spsp_DEU_SE$divisible_by_three<-ifelse(spsp_DEU_SE$genomicData.width %% 3 == 0,"Yes","No")
spsp_DEU_SE_row_annotation<-spsp_DEU_SE %>% dplyr::select(Row.names,groupID,PSI_class,mean_phastCon,Pfam_UCSC,divisible_by_three,cds)

# Clean annotation df
rownames(spsp_DEU_SE_row_annotation)<-spsp_DEU_SE_row_annotation$Row.names
spsp_DEU_SE_row_annotation<-spsp_DEU_SE_row_annotation %>% dplyr::select(-Row.names,-groupID)
spsp_DEU_SE_row_annotation$PSI_class<-gsub(" ","_",spsp_DEU_SE_row_annotation$PSI_class)
spsp_DEU_SE_row_annotation$cds<-gsub("no_CDS","No",spsp_DEU_SE_row_annotation$cds)
spsp_DEU_SE_row_annotation$cds<-gsub("CDS","Yes",spsp_DEU_SE_row_annotation$cds)
colnames(spsp_DEU_SE_row_annotation)<-c("PSI","phastCons","Pfam_UCSC","divisible_by_3","CDS")
spsp_DEU_SE_row_annotation <- spsp_DEU_SE_row_annotation %>% mutate(phastCons = case_when(phastCons<=0.3 ~ "Low",phastCons>=0.8 ~ "High",TRUE ~ "Mid"))
spsp_DEU_SE_row_annotation$phastCons<-as.factor(spsp_DEU_SE_row_annotation$phastCons)

# Set colors for annotation df
my_colour<-list(
PSI = c(Highly_excluded = "#F2F0F7", Mid_excluded = "#CBC9E2", Mid_included = "#9E9AC8" , Highly_included = "#6A51A3"),
CDS = c(Yes = "deeppink4", No = "#FFF5F0"),
phastCons = c(Low = "#FDE0DD", Mid="#FA9FB5", High="#C51B8A"),
divisible_by_3 = c(Yes = "darkslategray3", No = "azure1"),
Pfam_UCSC=c(Yes="lightblue",No="lightcyan" ) )

# Represent heatmap (scaled by row) for SE events
pheatmap(spsp_DEU_SE_plot,scale = "row",show_rownames = FALSE,cellwidth=80,cellheight = 2,cluster_cols=TRUE,cluster_rows=FALSE,color=colorRampPalette(c("#9eb3c2","#065a82","#21295c"))(100),annotation_row = spsp_DEU_SE_row_annotation,annotation_colors = my_colour,angle_col = 45,fontsize_col = 15,fontsize = 15)

#############################
### Annotate Alu elements ###
#############################

# Alu elements in SE exons.
# Allow intersection in any strand, since Alus can modify splicing in any orientation (especially if antisense!)
# Linux command-line
# grep "Alu" hg38_repeatmasker.bed  > hg38_Alu.bed
# bedtools intersect -c -a Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand.bed -b hg38_Alu.bed > Alu_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed

Alu<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/Alu_intersection_Homo_sapiens_5_sp_ens_exonic_parts_filtered.bed",sep="\t")
Alu<-Alu %>% dplyr::select(V4,V7)
Alu$V7<-ifelse(Alu$V7>0,"Yes","No")
colnames(Alu)<-c("Row.names","Alu")
# Add Alu overlap for each exonic part
dxr1_df<-merge(dxr1_df,Alu)

##########################################################################
# Functional enrichment for genes showing species-specific DEU changes ###
##########################################################################

# Concatenate all genes with species-specific DEU changes and define the universe of genes evaluated in DEU
spsp_DEU<-c(hs_genes,cs_genes,gs_genes,os_genes,ms_genes) %>% unique()
universe_DEU<- dxr1_df$groupID %>% unique()

WebGestaltR::WebGestaltR(projectName="WebGestalt_spsp_DEU_GO",interestGene = spsp_DEU,referenceGene = universe_DEU,interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")
WebGestaltR::WebGestaltR(projectName="WebGestalt_spsp_DEU_Panther",interestGene = spsp_DEU,referenceGene = universe_DEU,interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

####################################################################################################
# For the following analyses, conserved exonic parts are also required to pass expression cutoff ###
####################################################################################################

# We will remove the 'unchanged' that are exonBaseMean le than 10 (i.e. keeping exonBaseMean>10) leaving the rest as they are
dxr1_df<-dxr1_df %>% dplyr::filter(!(category=="Unchanged" & exonBaseMean<=10))

# After this filtering, there are still 61 exonic parts with NAs. They are classified as 'Unchanged' (and 9 of them show spurious differences in only 1 comp)
# Also remove them
dxr1_df<-dxr1_df %>% filter(!is.na(Homo_sapiens)) # 143,430 exonic parts remaining

#######################################################
### Gene-level classification based on DEU patterns ###
#######################################################

# Perform a more accurate classification of exonic parts depending on DEU. Based on these classes, provide a gene-level classification (Methods).

# Group 1: genes where all exonic parts are non-significant in all species (all_conserved)
# Group 2: genes where all exonic parts are either human specific up (single species) or conserved or other (only_human_sp)
# Group 3: genes where all isoforms are either NHP specific up (single species) or conserved or other (only_NHP_sp)
# Group 4: genes with isoforms with spsp DEU up in several species or conserved or other (convergence_spsp)
# Group 5: genes with isoforms with 3+2 patterns/spsp down together by conserved (other)

# Start from global results matrix
classify_DEU<-dxr1_df %>% dplyr::select(Row.names,groupID,Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Pongo_abelii,Macaca_mulatta,category)

# The 'other' category will have 3+2 changes and spsp-down
# Define 3+2 categories (10 in total) by hand for each exonic part
h_c<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Gorilla_gorilla_dxr1','Homo_sapiens_Pongo_abelii_dxr1',
'Homo_sapiens_Macaca_mulatta_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1','Pan_troglodytes_Pongo_abelii_dxr1','Pan_troglodytes_Macaca_mulatta_dxr1') 
h_g<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Homo_sapiens_Pongo_abelii_dxr1',
'Homo_sapiens_Macaca_mulatta_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1','Gorilla_gorilla_Pongo_abelii_dxr1','Gorilla_gorilla_Macaca_mulatta_dxr1') 
h_o<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Homo_sapiens_Gorilla_gorilla_dxr1',
'Homo_sapiens_Macaca_mulatta_dxr1','Pan_troglodytes_Pongo_abelii_dxr1','Gorilla_gorilla_Pongo_abelii_dxr1','Pongo_abelii_Macaca_mulatta_dxr1')
h_m<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Homo_sapiens_Gorilla_gorilla_dxr1',
'Homo_sapiens_Pongo_abelii_dxr1','Pan_troglodytes_Macaca_mulatta_dxr1','Gorilla_gorilla_Macaca_mulatta_dxr1','Pongo_abelii_Macaca_mulatta_dxr1')
c_g<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Pan_troglodytes_Pongo_abelii_dxr1',
'Pan_troglodytes_Macaca_mulatta_dxr1','Homo_sapiens_Gorilla_gorilla_dxr1','Gorilla_gorilla_Pongo_abelii_dxr1','Gorilla_gorilla_Macaca_mulatta_dxr1')
c_o<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1',
'Pan_troglodytes_Macaca_mulatta_dxr1','Homo_sapiens_Pongo_abelii_dxr1','Gorilla_gorilla_Pongo_abelii_dxr1','Pongo_abelii_Macaca_mulatta_dxr1') 
c_m<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pan_troglodytes_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1',
'Pan_troglodytes_Pongo_abelii_dxr1','Homo_sapiens_Macaca_mulatta_dxr1','Gorilla_gorilla_Macaca_mulatta_dxr1','Pongo_abelii_Macaca_mulatta_dxr1')
g_o<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Gorilla_gorilla_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1',
'Gorilla_gorilla_Macaca_mulatta_dxr1','Homo_sapiens_Pongo_abelii_dxr1','Pan_troglodytes_Pongo_abelii_dxr1','Pongo_abelii_Macaca_mulatta_dxr1')
g_m<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Gorilla_gorilla_dxr1','Pan_troglodytes_Gorilla_gorilla_dxr1',
'Gorilla_gorilla_Pongo_abelii_dxr1','Homo_sapiens_Macaca_mulatta_dxr1','Pan_troglodytes_Macaca_mulatta_dxr1','Pongo_abelii_Macaca_mulatta_dxr1')
o_m<-get_intersect_members(fromList(significant_exonic_parts_10_cmp),'Homo_sapiens_Pongo_abelii_dxr1','Pan_troglodytes_Pongo_abelii_dxr1',
'Gorilla_gorilla_Pongo_abelii_dxr1','Homo_sapiens_Macaca_mulatta_dxr1','Pan_troglodytes_Macaca_mulatta_dxr1','Gorilla_gorilla_Macaca_mulatta_dxr1')

# Groups can be checked using upset and pheatmap
#ok<-fromList(significant_exonic_parts_10_cmp)
#ok<-ok[rowSums(ok)=="6",]
#upset(ok,nsets=10,nintersects=1000,order.by="freq")
#classify_DEU[classify_DEU$Row.names %in% rownames(h_c),] %>% dplyr::select(Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Pongo_abelii,Macaca_mulatta) %>% pheatmap()

# Define spsp-DEU that are down-regulated
hs_down<-classify_DEU %>% dplyr::filter(category=="HS" & Homo_sapiens<Pan_troglodytes & Homo_sapiens<Gorilla_gorilla & Homo_sapiens<Pongo_abelii & Homo_sapiens<Macaca_mulatta)
cs_down<-classify_DEU %>% dplyr::filter(category=="CS" & Pan_troglodytes<Homo_sapiens & Pan_troglodytes<Gorilla_gorilla & Pan_troglodytes<Pongo_abelii & Pan_troglodytes<Macaca_mulatta)
gs_down<-classify_DEU %>% dplyr::filter(category=="GS" & Gorilla_gorilla<Homo_sapiens & Gorilla_gorilla<Pan_troglodytes & Gorilla_gorilla<Pongo_abelii & Gorilla_gorilla<Macaca_mulatta)
os_down<-classify_DEU %>% dplyr::filter(category=="OS" & Pongo_abelii<Homo_sapiens & Pongo_abelii<Pan_troglodytes & Pongo_abelii<Gorilla_gorilla & Pongo_abelii<Macaca_mulatta)
ms_down<-classify_DEU %>% dplyr::filter(category=="MS" & Macaca_mulatta<Homo_sapiens & Macaca_mulatta<Pan_troglodytes & Macaca_mulatta<Gorilla_gorilla & Macaca_mulatta<Pongo_abelii)

# Exonic parts with 3+2 pattern or spsp-down make group 'other'
other<-c(rownames(h_c),rownames(h_g),rownames(h_o),rownames(h_m),rownames(c_g),rownames(c_o),rownames(c_m),rownames(g_o),rownames(g_m),rownames(o_m),
hs_down$Row.names,cs_down$Row.names,gs_down$Row.names,os_down$Row.names,ms_down$Row.names)

# Classify exonic parts:
classify_DEU <- classify_DEU %>% mutate(isoform_class = case_when(category=="HS" & Homo_sapiens>Pan_troglodytes & Homo_sapiens>Gorilla_gorilla & Homo_sapiens>Pongo_abelii & Homo_sapiens>Macaca_mulatta ~ "human_up",
category=="CS" & Pan_troglodytes>Homo_sapiens & Pan_troglodytes>Gorilla_gorilla & Pan_troglodytes>Pongo_abelii & Pan_troglodytes>Macaca_mulatta ~ "chimpanzee_up",
category=="GS" & Gorilla_gorilla>Homo_sapiens & Gorilla_gorilla>Pan_troglodytes & Gorilla_gorilla>Pongo_abelii & Gorilla_gorilla>Macaca_mulatta ~ "gorilla_up",
category=="OS" & Pongo_abelii>Homo_sapiens & Pongo_abelii>Pan_troglodytes & Pongo_abelii>Gorilla_gorilla & Pongo_abelii>Macaca_mulatta ~ "orangutan_up",
category=="MS" & Macaca_mulatta>Homo_sapiens & Macaca_mulatta>Pan_troglodytes & Macaca_mulatta>Gorilla_gorilla & Macaca_mulatta>Pongo_abelii ~ "macaque_up",
category=="Unchanged" ~ "conserved",
Row.names %in% other ~ "other"
))

# At this point, 'NA' correspond to exonic parts that are neither conserved, or spsp-up or other (3+2 or spsp-downs). Ignore these changes (variability, see Methods).
classify_DEU<-classify_DEU[!is.na(classify_DEU$isoform_class),]
# Simplify df
classify_DEU<-classify_DEU %>% dplyr::select(groupID,isoform_class)

# Aggregate exonic part classification at the level of genes (alphabetical sorting)
library("data.table")
classify_DEU<-as.data.table(classify_DEU)[, toString(sort(unique(isoform_class))), by = "groupID"]
colnames(classify_DEU)<-c("gene_id","isoform_class")

# Based on exonic part classes for each gene, perform gene-level classification
classify_DEU<-classify_DEU %>% dplyr::mutate(
gene_class = case_when(isoform_class=="conserved" ~ "all_conserved", 
isoform_class=="human_up" | isoform_class=="conserved, human_up" | isoform_class=="conserved, human_up, other" | isoform_class=="human_up, other"   ~ "only_human_up",
isoform_class=="chimpanzee_up" | isoform_class=="chimpanzee_up, conserved" | isoform_class=="chimpanzee_up, conserved, other" | isoform_class=="chimpanzee_up, other"  ~ "only_NHP_up",
isoform_class=="gorilla_up" | isoform_class=="conserved, gorilla_up" | isoform_class=="conserved, gorilla_up, other" | isoform_class=="gorilla_up, other" ~ "only_NHP_up",
isoform_class=="orangutan_up" | isoform_class=="conserved, orangutan_up" | isoform_class=="conserved, orangutan_up, other" | isoform_class=="orangutan_up, other" ~ "only_NHP_up",
isoform_class=="macaque_up" | isoform_class=="conserved, macaque_up" | isoform_class=="conserved, macaque_up, other" | isoform_class=="macaque_up, other" ~ "only_NHP_up",
isoform_class=="other" | isoform_class=="conserved, other" ~ "other",
TRUE ~ "convergence_spsp"
))

# Inspect human population-based dN/dS from https://pubmed.ncbi.nlm.nih.gov/28977405/ across gene classes defined by DEU
dN_dS_human<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/dNdS_primates/ExAc_human_dNdS.txt",sep="\t",header=TRUE)
dN_dS_human<-dN_dS_human %>% dplyr::select(Ensembl.Gene.ID,ExAC.unique_variants.dN_dS)
classify_DEU_with_dN_dS<-merge(classify_DEU,dN_dS_human,by.x="gene_id",by.y="Ensembl.Gene.ID")
classify_DEU_with_dN_dS$gene_class <- factor(classify_DEU_with_dN_dS$gene_class, levels = c("all_conserved", "only_human_up", "only_NHP_up","other","convergence_spsp"))

ggboxplot(classify_DEU_with_dN_dS,x="gene_class",y="ExAC.unique_variants.dN_dS",fill="gene_class") + coord_cartesian(ylim=c(0,1)) + geom_violin(alpha=0.1)

#############################################################################################
# Inspect phastCons scores and pi diversity based on exonic parts classification
# Also inspect overlap with Alu elements, CDS and APPRIS principal isoforms (focus on human-specific changes)
#############################################################################################

# For this, we will use 3 groups of exonic parts: conserved, spsp-up and -down and 'other' (3+2 changes)
# We exclude non-grouped changes since they are due to variability (Methods)

# Define spsp-DEU-upregulated
hs_up<- dxr1_df %>% dplyr::filter(category=="HS" & Homo_sapiens>Pan_troglodytes & Homo_sapiens>Gorilla_gorilla & Homo_sapiens>Pongo_abelii & Homo_sapiens>Macaca_mulatta)
cs_up<- dxr1_df %>% dplyr::filter(category=="CS" & Pan_troglodytes>Homo_sapiens & Pan_troglodytes>Gorilla_gorilla & Pan_troglodytes>Pongo_abelii & Pan_troglodytes>Macaca_mulatta)
gs_up<- dxr1_df %>% dplyr::filter(category=="GS" & Gorilla_gorilla>Homo_sapiens & Gorilla_gorilla>Pan_troglodytes & Gorilla_gorilla>Pongo_abelii & Gorilla_gorilla>Macaca_mulatta)
os_up<- dxr1_df %>% dplyr::filter(category=="OS" & Pongo_abelii>Homo_sapiens & Pongo_abelii>Pan_troglodytes & Pongo_abelii>Gorilla_gorilla & Pongo_abelii>Macaca_mulatta)
ms_up<- dxr1_df %>% dplyr::filter(category=="MS" & Macaca_mulatta>Homo_sapiens & Macaca_mulatta>Pan_troglodytes & Macaca_mulatta>Gorilla_gorilla & Macaca_mulatta>Pongo_abelii)

# Add new column with more detailed exonic part classes (phastCons and pi diversity analyses)
dxr1_df<-dxr1_df %>% mutate(fine_category = case_when(
Row.names %in% rownames(h_c) ~ "other",
Row.names %in% rownames(h_g) ~ "other",
Row.names %in% rownames(h_o) ~ "other",
Row.names %in% rownames(h_m) ~ "other",
Row.names %in% rownames(c_g) ~ "other",
Row.names %in% rownames(c_o) ~ "other",
Row.names %in% rownames(c_m) ~ "other",
Row.names %in% rownames(g_o) ~ "other",
Row.names %in% rownames(g_m) ~ "other",
Row.names %in% rownames(o_m) ~ "other",
Row.names %in% hs_up$Row.names ~ "spsp_up",
Row.names %in% cs_up$Row.names ~ "spsp_up",
Row.names %in% gs_up$Row.names ~ "spsp_up",
Row.names %in% os_up$Row.names ~ "spsp_up",
Row.names %in% ms_up$Row.names ~ "spsp_up",
Row.names %in% hs_down$Row.names ~ "spsp_down",
Row.names %in% cs_down$Row.names ~ "spsp_down",
Row.names %in% gs_down$Row.names ~ "spsp_down",
Row.names %in% os_down$Row.names ~ "spsp_down",
Row.names %in% ms_down$Row.names ~ "spsp_down",
category=="Unchanged" ~ "conserved",
TRUE ~ "non_group_changes"))

# Here 'non_group_changes' are exonic parts whose usage differences doesn't form consistent groups. Also sp-sp that are not consistently up or down (only 5)
# Not included in the plot
dxr1_df2<-dxr1_df %>% dplyr::filter(fine_category!="non_group_changes")

#######################
# 1000Genomes masking #
#######################

# How many exonic parts are confidently called in 1000Genomes? Select those where strict masking overlaps all bases within the exonic part
mask<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/pi_diversity_GRCh38/Homo_sapiens_5_sp_ens_exonic_parts_filtered_w_strand_overlap_with_passed_bases.txt",sep="\t")
length_bases<-mask %>% dplyr::mutate(length=V3-V2) %>% dplyr::select(V4,length) %>% unique()
colnames(length_bases)<-c("exonic_part","length")
sum_passed_bases<-aggregate(mask$V11~mask$V4,FUN=sum)
colnames(sum_passed_bases)<-c("exonic_part","passed_bases")
passed_bases_over_total_length<-merge(sum_passed_bases,length_bases)
passed_bases_over_total_length<-passed_bases_over_total_length %>% dplyr::mutate(percentage_passed=(passed_bases/length)*100)
# hist(passed_bases_over_total_length$percentage_passed) # most of them are fully covered by passed bases in strict mask (1000 Genomes)

# Exonic parts to be represented must pass 100% coverage by passed bases (high-confidence).
passed_bases_OK<-passed_bases_over_total_length %>% dplyr::filter(percentage_passed==100) %>% dplyr::select(exonic_part)

# Average pi=sum divided by length. Since we are restricting to exonic parts fully covered by 1000 Genomes, length==length of passed bases
pi<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DEU/liftOver_merged_models/present_in_5_species/defining_genes_by_ensembl/quantifying_filtered_exonic_parts/results/matching_order/DEXSeq/pi_diversity_GRCh38/my_regions_sum_diversity_and_length.txt", sep="\t",header=FALSE)
colnames(pi)<-c("exonic_part","sum_pi","length")
pi<-pi %>% dplyr::mutate(pi_estimation=sum_pi/length) %>% dplyr::select(exonic_part,pi_estimation)
# Add pi with 0 if there is no pi estimation (no variation)
dxr1_df2<-left_join(dxr1_df2, pi, by = c("Row.names" = "exonic_part"))
dxr1_df2[is.na(dxr1_df2$pi_estimation),"pi_estimation"]<-0
# Restrict to exonic parts in passed_bases_OK, thus keeping exonic parts completely callable in 1000Genomes (100% of passed bases)
dxr1_df2<-dxr1_df2 %>% dplyr::filter(Row.names %in% passed_bases_OK$exonic_part)
# Sort levels
dxr1_df2$fine_category <- factor(dxr1_df2$fine_category, levels = c("conserved","spsp_up","spsp_down","other"))
# Define comparisons for statistical testing
my_comparisons <- list( c("conserved", "spsp_up"), c("conserved", "spsp_down"), c("conserved", "other"), c("spsp_up", "spsp_down"), c("spsp_up", "other"), c("spsp_down", "other"))
# Limit to exonic parts with >5 bases, save it as dxr1_df3.
dxr1_df3<-dxr1_df2 %>% dplyr::filter(genomicData.width>5)

################################
# PhastCons boxplot (>5 bases) #
################################

ggboxplot(dxr1_df3,x="fine_category",y="mean_phastCon",fill="fine_category",outlier.size=0.001) + coord_cartesian(ylim=c(0,1.7)) + theme(aspect.ratio=0.8) + geom_violin(alpha=0.2,aes(fill=fine_category)) + scale_fill_nejm() + stat_compare_means(comparisons = my_comparisons,aes(label = ..p.signif..)) + xlab("Exonic parts") + ylab("Average phastCons (17 primates)") +  theme(text = element_text(size=20)) + theme(legend.title=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

compare_means(data=dxr1_df3,mean_phastCon~fine_category) 

###########################
# Pi boxplot (>5 bases) ###
###########################

ggboxplot(dxr1_df3,x="fine_category",y="pi_estimation",fill="fine_category",outlier.shape=NA) + coord_cartesian(ylim=c(0,0.004)) + theme(aspect.ratio=0.8)+ scale_fill_nejm() + stat_compare_means(comparisons = my_comparisons,aes(label = ..p.signif..)) + xlab("Exonic parts") + ylab("Pi diversity (1000 Genomes)") +  theme(text = element_text(size=20)) + theme(legend.title=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

compare_means(data=dxr1_df3,pi_estimation~fine_category)

###################################
### Alu, CDS and APPRIS ###########
###################################

ggthemr::ggthemr_reset()

# Alu, CDS and APPRIS overlap will be compared between human-specific DEU (up or down) and the conserved exonic parts
# Specify species and change of direction for exonic parts showing species-specific DEU.

dxr1_df<-dxr1_df %>% mutate(fine_category = case_when(
Row.names %in% rownames(h_c) ~ "other",
Row.names %in% rownames(h_g) ~ "other",
Row.names %in% rownames(h_o) ~ "other",
Row.names %in% rownames(h_m) ~ "other",
Row.names %in% rownames(c_g) ~ "other",
Row.names %in% rownames(c_o) ~ "other",
Row.names %in% rownames(c_m) ~ "other",
Row.names %in% rownames(g_o) ~ "other",
Row.names %in% rownames(g_m) ~ "other",
Row.names %in% rownames(o_m) ~ "other",
Row.names %in% hs_up$Row.names ~ "human_up",
Row.names %in% cs_up$Row.names ~ "chimpanzee_up",
Row.names %in% gs_up$Row.names ~ "gorilla_up",
Row.names %in% os_up$Row.names ~ "orangutan_up",
Row.names %in% ms_up$Row.names ~ "macaque_up",
Row.names %in% hs_down$Row.names ~ "human_down",
Row.names %in% cs_down$Row.names ~ "chimpanzee_down",
Row.names %in% gs_down$Row.names ~ "gorilla_down",
Row.names %in% os_down$Row.names ~ "orangutan_down",
Row.names %in% ms_down$Row.names ~ "macaque_down",
category=="Unchanged" ~ "conserved",
TRUE ~ "non_group_changes"))

# Alu enrichment #
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down" | fine_category=="human_up") %>% droplevels() %>% group_by(fine_category) %>% dplyr::count(Alu) %>% mutate(prop=n/(sum(n))) %>% ggplot(data=., aes(x=fine_category,y=prop,fill=Alu)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab("Exonic parts") + ylab("Proportion") + scale_fill_manual(values=c("grey60", "red2"))
# Fisher tests
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value = 1.386e-10 --> x3 = 4.16e-10
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="conserved") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value=1.093231e-18 --> x3 = 3.28e-18
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value = 0.8922 --> x3 = 1
# Alu enrichment for skipping exons #
dxr1_df %>% dplyr::filter(AS_class=="SE") %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down" | fine_category=="human_up") %>% droplevels() %>% group_by(fine_category,AS_class) %>% dplyr::count(Alu) %>% mutate(prop=n/(sum(n))) %>% ggplot(data=., aes(x=fine_category,y=prop,fill=Alu)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab("Skipping exons (SE)") + ylab("Proportion") + scale_fill_manual(values=c("grey60", "red2"))
# Fisher tests
dxr1_df %>% dplyr::filter(AS_class=="SE") %>% dplyr::filter(fine_category=="human_up" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value = 0.003889 --> x3 = 0.012
dxr1_df %>% dplyr::filter(AS_class=="SE") %>% dplyr::filter(fine_category=="human_up" | fine_category=="conserved") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value = 1.88e-15 --> x3 = 5.64e-15
dxr1_df %>% dplyr::filter(AS_class=="SE") %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(Alu,fine_category) %>% table() %>% fisher.test() # p-value = 1 --> x3 = 1

# APPRIS principal isoforms' enrichment #
dxr1_df<-dxr1_df %>% mutate(overlap_principal_transcript = case_when(appris_overlap== "Overlap_principal_and_minor" ~ "Yes",appris_overlap=="Only_principal" ~ "Yes",TRUE ~ "No")) # overlaps principal transcript or not
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down" | fine_category=="human_up") %>% droplevels() %>% group_by(fine_category) %>% dplyr::count(overlap_principal_transcript) %>% mutate(prop=n/(sum(n))) %>% ggplot(data=., aes(x=fine_category,y=prop,fill=overlap_principal_transcript)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab("Exonic parts") + ylab("Proportion") + scale_fill_manual(values=c("grey60", "red2"))
# Fisher tests
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(overlap_principal_transcript,fine_category) %>% table() %>% fisher.test() #  p-value = 1.229525e-28 --> x3 = 3.69e-28
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="conserved") %>% droplevels() %>% dplyr::select(overlap_principal_transcript,fine_category) %>% table() %>% fisher.test() # p-value = 1.002362e-47 --> x3 = 3.01e-47
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(overlap_principal_transcript,fine_category) %>% table() %>% fisher.test() # p-value = 0.574 --> x3 = 1

# CDS enrichment #
dxr1_df$cds <- factor(dxr1_df$cds, levels = c("no_CDS","CDS")) # set order of levels
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down" | fine_category=="human_up") %>% droplevels() %>% group_by(fine_category) %>% dplyr::count(cds) %>% mutate(prop=n/(sum(n))) %>% ggplot(data=., aes(x=fine_category,y=prop,fill=cds)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab("Exonic parts") + ylab("Proportion") + scale_fill_manual(values=c("grey60", "red2"))
# Fisher tests
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(cds,fine_category) %>% table() %>% fisher.test() # p-value = 1.481e-12 --> x3 = 4.44e-12
dxr1_df %>% dplyr::filter(fine_category=="human_up" | fine_category=="conserved") %>% droplevels() %>% dplyr::select(cds,fine_category) %>% table() %>% fisher.test() # p-value = 1.115269e-19 --> x3 = 3.35e-19
dxr1_df %>% dplyr::filter(fine_category=="conserved" | fine_category=="human_down") %>% droplevels() %>% dplyr::select(cds,fine_category) %>% table() %>% fisher.test() # p-value = 0.8699 --> x3 = 1
# All p-values*3 (3 comparatives) (Bonferroni correction).


