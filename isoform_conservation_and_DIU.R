#############################################################################################
# Main code for isoform expression gains and losses and differential isoform usage analyses #
#############################################################################################

library("data.table")
library("dplyr")
library("sva")
library("vegan")
library("tweeDEseq")
library("stringr")
library("UpSetR")
library("ggsci")

###############################
### Isoform TPM calculation ###
###############################

# Read Kallisto abundance files per LCL
myfiles<-list.files(path="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/kallisto",pattern = "abundance.tsv",recursive = TRUE,full.names = TRUE) 
samples<-gsub(".+/(.+)/.+","\\1",myfiles)
list.DFs <- lapply(myfiles,fread)
est_counts <- lapply(list.DFs, function(x) x%>% select(target_id,est_counts))
for(i in 1:16) { colnames(est_counts[[i]])[2]<-samples[i] }
est_counts<-Reduce(function(x,y)merge(x,y,by="target_id",all=TRUE), est_counts)
est_counts<-as.data.frame(est_counts)
rownames(est_counts)<-est_counts$target_id
est_counts<-est_counts %>% select(-target_id)
# All 'NA' (isoform was not quantified by Kallisto in that transcriptome) will be assigned 0 counts
est_counts[is.na(est_counts)] <- 0
# Define experimental batches
batch<-c("R","L","R","L","R","L","L","S","S","R","L","R","R","R","R","L")
# ComBat batch-effect correction
corrected_counts<-sva::ComBat_seq(as.matrix(est_counts),batch = batch, group=NULL)
corrected_counts<-as.data.frame(corrected_counts)
# TMM-normalisation of counts here
TMM_counts <- normalizeCounts(corrected_counts,method="TMM")
TMM_counts<-as.data.frame(TMM_counts)

# TPM calculation #
eff_length <- lapply(list.DFs, function(x) x%>% select(target_id,eff_length))
for(i in 1:16) { colnames(eff_length[[i]])[2]<-samples[i] }
eff_length<-Reduce(function(x,y)merge(x,y,by="target_id",all=TRUE), eff_length)
eff_length<-as.data.frame(eff_length)
rownames(eff_length)<-eff_length$target_id
eff_length<-eff_length %>% select(-target_id)
# All 'NA' (isoform was not quantified by Kallisto in that transcriptome) will be assigned eff_lenght=1 to avoid 0/0 division)
eff_length[is.na(eff_length)] <- 1
# TPM calculation function
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# Compute TPM from batch-effect corrected counts and effective lengths of isoform in each species
TPM<-tpm3(corrected_counts,eff_length)
TPM<-as.data.frame(TPM)
# Remove GM12878 technical replicate 1 (only used for batch effect correction)
TPM<-TPM %>% select(-GM12878.rep1)

###########################################################
###### PCA before and after batch-effect correction #######
###########################################################

# TPM for est_counts (pre-ComBat), corrected_counts (after ComBat) and TMM_counts (after ComBat and TMM)
library("PCAtools")

# 1) est_counts
TPM_est_counts<-tpm3(est_counts,eff_length)
TPM_est_counts<-as.data.frame(TPM_est_counts)
TPM_est_counts<-TPM_est_counts %>% select(-GM12878.rep1)

metadata<-data.frame(row.names=colnames(TPM_est_counts))
metadata$samples<-gsub("(.+).rep[123]","\\1",rownames(metadata))
metadata$species<-c("Chimpanzee","Chimpanzee","Chimpanzee","Orangutan","Gorilla","Orangutan","Gorilla","Human","Human","Human","Gorilla","Orangutan","Macaque","Macaque","Macaque")
metadata$species <- factor(metadata$species, levels = c("Human","Chimpanzee","Gorilla","Orangutan","Macaque"))
metadata$batch<-c("B1","B2","B1","B2","B1","B2","B2","B3","B1","B2","B1","B1","B1","B1","B2")
# log(TPM+1) PCA
p <- pca(as.matrix(log(TPM_est_counts+1)), metadata = metadata, removeVar = 0.1)
biplot(p,legendPosition='right',colby = 'species',pointSize = 5,title='Isoform expression (log TPM + 1) without batch effect correction',lab=p$metadata$samples,shape = 'batch') + theme_bw() + theme(aspect.ratio=0.9) + scale_color_npg(alpha=0.9) + theme(plot.title = element_text(hjust = 0.5,vjust=0.5)) + coord_cartesian(xlim=c(-200,200)) + labs(colour="Species",shape="Batch") + theme(text=element_text(size=15))

# 2) corrected counts
TPM_corrected_counts<-tpm3(corrected_counts,eff_length)
TPM_corrected_counts<-as.data.frame(TPM_corrected_counts)
TPM_corrected_counts<-TPM_corrected_counts %>% select(-GM12878.rep1)
# log(TPM+1) PCA
p <- pca(as.matrix(log(TPM_corrected_counts+1)), metadata = metadata, removeVar = 0.1)
biplot(p,legendPosition='right',colby = 'species',pointSize = 5,title='Isoform expression (log TPM + 1) after ComBat correction',lab=p$metadata$samples,shape = 'batch') + theme_bw() + theme(aspect.ratio=0.9) + scale_color_npg(alpha=0.9) + theme(plot.title = element_text(hjust = 0.5,vjust=0.5)) + coord_cartesian(xlim=c(-200,200)) + labs(colour="Species",shape="Batch") + theme(text=element_text(size=15))

# 3) TMM counts
TPM_TMM_counts<-tpm3(TMM_counts,eff_length)
TPM_TMM_counts<-as.data.frame(TPM_TMM_counts)
TPM_TMM_counts<-TPM_TMM_counts %>% select(-GM12878.rep1)
# log(TPM+1) PCA
p <- pca(as.matrix(log(TPM_TMM_counts+1)), metadata = metadata, removeVar = 0.1)
biplot(p,legendPosition='right',colby = 'species',pointSize = 5,title='Isoform expression (log TPM + 1) after ComBat and TMM normalization',lab=p$metadata$samples,shape = 'batch') + theme_bw() + theme(aspect.ratio=0.9) + scale_color_npg(alpha=0.9) + theme(plot.title = element_text(hjust = 0.5,vjust=0.5)) + coord_cartesian(xlim=c(-200,200)) + labs(colour="Species",shape="Batch") + theme(text=element_text(size=15)) 

######################################################
### Keep only intra-species consistent transcripts ###
######################################################

# Format columns
colnames(TPM)<-c("Pan_troglodytes_CH114","Pan_troglodytes_CH170","Pan_troglodytes_CH391","Pongo_abelii_CRL.1850","Gorilla_gorilla_DIAN","Pongo_abelii_EB185","Gorilla_gorilla_GG05","Homo_sapiens_GM12878","Homo_sapiens_GM19150","Homo_sapiens_GM19238","Gorilla_gorilla_OMOYE","Pongo_abelii_PPY6_1","Macaca_mulatta_R02027","Macaca_mulatta_R05040","Macaca_mulatta_R94011")
# TPM cutoff is 0
my_threshold<-TPM>0
my_threshold<-as.data.frame(my_threshold)
# Remove intra-species inconsistencies(example: one human sample is TRUE and two human are FALSE)
copy<-my_threshold
copy$isoform<-rownames(copy)
copy<-reshape2::melt(copy,id.vars="isoform")
copy<-copy %>% dplyr::mutate(species = str_extract(variable, ".+[^_[.0-9A-Z_]+]"))
colnames(copy)<-c("isoform","sample","expressed","species")
# If TRUE/FALSE are both present in isoform AND species, then discard
inconsistent_isoforms<-copy %>% group_by(isoform,species) %>% mutate(n_distinct=n_distinct(expressed)) %>% ungroup() %>% filter(n_distinct=="2")
# Remove inconsistencies from 'my_threshold' matrix
my_threshold<-my_threshold %>% filter(!rownames(.) %in% inconsistent_isoforms$isoform)
# Save expressed isoforms in each sample in different objects
for (name in colnames(my_threshold)) { 
	assign(x=name,value=row.names(my_threshold)[my_threshold[[name]]])
}

# Upset plot for intra-sp consistent transcripts	
listInput <- list(Homo_sapiens_GM19150 = Homo_sapiens_GM19150, 
Homo_sapiens_GM19238 = Homo_sapiens_GM19238,
Homo_sapiens_GM12878 = Homo_sapiens_GM12878,
Pan_troglodytes_CH114 = Pan_troglodytes_CH114,
Pan_troglodytes_CH391 = Pan_troglodytes_CH391,
Pan_troglodytes_CH170 = Pan_troglodytes_CH170,
Gorilla_gorilla_DIAN = Gorilla_gorilla_DIAN,
Gorilla_gorilla_OMOYE = Gorilla_gorilla_OMOYE,
Gorilla_gorilla_GG05 = Gorilla_gorilla_GG05,
Pongo_abelii_PPY6_1 = Pongo_abelii_PPY6_1,
Pongo_abelii_EB185 = Pongo_abelii_EB185,
Pongo_abelii_CRL.1850 = Pongo_abelii_CRL.1850,
Macaca_mulatta_R02027 = Macaca_mulatta_R02027,
Macaca_mulatta_R05040 = Macaca_mulatta_R05040,
Macaca_mulatta_R94011 = Macaca_mulatta_R94011)
# Set colors
colors<-list(rep("grey80",3),rep("#BDD7E7",3),rep("#6BAED6",3),rep("#3182BD",3),rep("#08519C",3))
colors<-unlist(colors)
# Plotting isoform intersection sets
upset(fromList(listInput),nintersects=62,nsets = 15,keep.order=TRUE,sets.bar.color = colors,order.by="freq",sets.x.label="Transcripts above 0 TPM",
sets = c("Macaca_mulatta_R02027",
"Macaca_mulatta_R05040",
"Macaca_mulatta_R94011",
"Pongo_abelii_PPY6_1",
"Pongo_abelii_EB185",
"Pongo_abelii_CRL.1850",
"Gorilla_gorilla_DIAN",
"Gorilla_gorilla_OMOYE",
"Gorilla_gorilla_GG05",
"Pan_troglodytes_CH114",
"Pan_troglodytes_CH391",
"Pan_troglodytes_CH170",
"Homo_sapiens_GM19150",
"Homo_sapiens_GM12878",
"Homo_sapiens_GM19238"))

####################################
### Species-specific transcripts ###
####################################

library("reshape2")
library("ggplot2")

list_species<-c("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta")
TPM_2<-TPM
# Get binary matrix of isoform expression
TPM_2$nb_samples<-rowSums(TPM_2>0)

# Dataframe with all orthologous transcripts classified by SQANTI in each genome. Add label for species-specific vs background.
for (species in list_species) {

	subdf<-TPM_2 %>% filter(nb_samples=="3") %>% select(matches(species))
	subdf$nb_samples<-rowSums(subdf>0)
	spsp<-subdf %>% filter(nb_samples=="3") %>% rownames()

	# Load SQANTI df computed from each species and its genome/annotation
	SQANTI<-paste("/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/",species,"/QC/",species,"_5_sp_SQANTI3_classification.txt",sep="")
	SQANTI<-read.table(file=SQANTI,sep="\t",header=TRUE)
	# Restrict to only intra-sp consistent transcripts and add condition of spsp for each species	
	SQANTI<-SQANTI %>% filter(isoform %in% rownames(my_threshold)) 
	SQANTI$spsp<-ifelse(SQANTI$isoform %in% spsp,"Species-specific","Background")
	SQANTI$associated_transcript2<-ifelse(SQANTI$associated_transcript=="novel","novel","annotated")
	SQANTI$species<-species

	SQANTI_name<-paste(species,"_SQANTI",sep="")
	assign(x=SQANTI_name,value=SQANTI)	
}

# Concatenate all SQANTI's and set levels
all_species_SQANTI<-rbind(Homo_sapiens_SQANTI,Pan_troglodytes_SQANTI,Gorilla_gorilla_SQANTI,Pongo_abelii_SQANTI,Macaca_mulatta_SQANTI)
all_species_SQANTI$species <- factor(all_species_SQANTI$species,levels=c("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta"))

# Some plots based on SQANTI features in each genome+reference annotation. Species-specific vs rest.

# Annotated versus novel
all_species_SQANTI %>% group_by(species,spsp) %>% count(associated_transcript2) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=associated_transcript2)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + guides(fill=guide_legend(title="Species reference annotation")) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# SQANTI subcategories
all_species_SQANTI %>% group_by(species,spsp) %>% count(subcategory) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=subcategory)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Splice junction canonicity
all_species_SQANTI %>% group_by(species,spsp) %>% count(all_canonical) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=all_canonical)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Coding potential
all_species_SQANTI %>% group_by(species,spsp) %>% count(coding) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=coding)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# RTS evaluation
all_species_SQANTI %>% group_by(species,spsp) %>% count(RTS_stage) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=RTS_stage)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# NMD prediction
all_species_SQANTI %>% group_by(species,spsp) %>% count(predicted_NMD) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=spsp,y=prop,fill=predicted_NMD)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Isoform lengths
ggplot(all_species_SQANTI, aes(x=spsp,y=length)) + geom_boxplot() + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + coord_cartesian(ylim=c(0,8000))
# Number of exons
ggplot(all_species_SQANTI, aes(x=spsp,y=exons)) + geom_boxplot() + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~species,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + coord_cartesian(ylim=c(0,20))
# Enrichment in novel isoforms in the set of spsp-transcripts compared to background
table(all_species_SQANTI$spsp,all_species_SQANTI$associated_transcript2) %>% fisher.test() %>% str() # OR=3.86, p-value=1.4e-117

# Evaluate species-specific isoforms in the context of other species' genome+annotation
spsp_df<-data.frame()
for (my_species in list_species) {
    my_spsp<-all_species_SQANTI %>% filter(species==my_species & spsp=="Species-specific") %>% select(isoform)
    spsp_df<-rbind(spsp_df,all_species_SQANTI %>% filter(isoform %in% my_spsp$isoform) %>% mutate(spsp_all=paste(my_species,"specific",sep="_")))
}
# Set order of levels
spsp_df$species <- factor(spsp_df$species,levels=c("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta"))
spsp_df$spsp_all <- factor(spsp_df$spsp_all,levels=c("Homo_sapiens_specific","Pan_troglodytes_specific","Gorilla_gorilla_specific","Pongo_abelii_specific","Macaca_mulatta_specific"))
# Species-specific transcripts increase proportion of non-canonical in the genome of other species
spsp_df %>% group_by(species,spsp_all) %>% count(all_canonical) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=species,y=prop,fill=all_canonical)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~spsp_all,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Genome and annotation")
# Also, some cases of species-specific that are annotated in other species
spsp_df %>% group_by(species,spsp_all) %>% count(associated_transcript2) %>% mutate(prop=n/sum(n)) %>% ggplot(., aes(x=species,y=prop,fill=associated_transcript2)) + geom_bar(position="stack", stat="identity",color="black",width = 0.5) + xlab("") + theme(aspect.ratio=0.8) + facet_wrap(~spsp_all,nrow=1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + guides(fill=guide_legend(title="Reference annotation")) + xlab("Genome and annotation")

# Add new column with species-specific label independent of the species
all_species_SQANTI$all_spsp<-ifelse(all_species_SQANTI$isoform %in% spsp_df$isoform,"Species-specific","Background")

######################################################################
### GO enrichment for genes producing species-specific transcripts ###
######################################################################

spsp<-all_species_SQANTI %>% filter(spsp=="Species-specific") %>% select(isoform)
# Get genes producing species-specific transcripts (human Ensembl gene IDs)
spsp_genes<-Homo_sapiens_SQANTI[Homo_sapiens_SQANTI$isoform %in% spsp$isoform,"associated_gene"] %>% unique()
# 5,627 genes evaluated (universe)
universe_genes<-Homo_sapiens_SQANTI$associated_gene %>% unique()

library("WebGestaltR")
WebGestaltR(projectName="WebGestalt_spsp_gains_GO",interestGene = as.vector(spsp_genes),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")
WebGestaltR(projectName="WebGestalt_spsp_gains_Panther",interestGene = as.vector(spsp_genes),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

# Same test removing spsp_DGE_up genes
# spsp_DGE_up<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/DGE/send_to_Luis/result/spsp_DGE_up.txt",sep="\t",header=FALSE)
#WebGestaltR(projectName="WebGestalt_spsp_gains_removing_spsp_DGE_up_GO",interestGene =  as.vector(setdiff(spsp_genes,spsp_DGE_up$V1)),referenceGene = as.vector(setdiff(universe_genes,spsp_DGE_up$V1)),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")
#WebGestaltR(projectName="WebGestalt_spsp_gains_removing_spsp_DGE_up_Panther",interestGene =  as.vector(setdiff(spsp_genes,spsp_DGE_up$V1)),referenceGene = as.vector(setdiff(universe_genes,spsp_DGE_up$V1)),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

####################################################################
### Tissue-specific expression for genes producing spsp-isoforms ###
####################################################################

library(TissueEnrich)
# Genes producing species-specific transcripts
gs<-GeneSet(geneIds=as.character(spsp_genes),organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())
# Genes in the background (all evaluated here)
bg<-GeneSet(geneIds=as.character(universe_genes),organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())
te_output<-teEnrichment(inputGenes = gs,backgroundGenes = bg,rnaSeqDataset = 2) # with GTEx dataset

seEnrichmentOutput<-te_output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
print(enrichmentOutput)

# Log2 FC by tissues
ggthemr::ggthemr_reset()
ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
geom_bar(stat = 'identity')+
labs(x='', y = '-LOG10(P-Adjusted)')+
theme_bw()+
theme(legend.position="none")+
theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

# How is the expression of spleen specific/enhanced genes with sp-specific isoforms across tissues?
library(tidyr)
seExp<-te_output[[2]][["Spleen"]]
exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
exp$Gene<-row.names(exp)
exp<-exp %>% gather(key = "Tissue", value = "expression",1:(ncol(exp)-1))

ggplot(exp, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+
labs(x='', y = '')+
theme_bw()+
guides(fill = guide_legend(title = "Log2(TPM)"))+
theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

# Against background, genes with sp-specific isoforms are ~2x enriched (p-value 0.00073) in spleen-specific genes (incl. enhanced and group-specific)
# According to https://tissueenrich.gdcb.iastate.edu/ -> recommended BH adj p-value cutoff is 0.01.

################################################################################################
### Binary matrix of isoform expression across species to generate COUNT tree (gains/losses) ###
################################################################################################

# Generate the 1/0 matrix for COUNT (only intra-sp consistent isoforms): 
binary_matrix_for_COUNT<-my_threshold %>% dplyr::select(Homo_sapiens_GM19150,Pan_troglodytes_CH170,Gorilla_gorilla_OMOYE,Pongo_abelii_EB185,Macaca_mulatta_R05040) # take any sample (expression is consistent in all LCLs per species)
colnames(binary_matrix_for_COUNT)<-c("Homo sapiens","Pan troglodytes","Gorilla gorilla","Pongo abelii","Macaca mulatta")
binary_matrix_for_COUNT<-1*binary_matrix_for_COUNT # convert to 1/0 (presence/absence)
binary_matrix_for_COUNT$family<-rownames(binary_matrix_for_COUNT)
binary_matrix_for_COUNT<-binary_matrix_for_COUNT %>% dplyr::select("family","Homo sapiens","Pan troglodytes","Gorilla gorilla","Pongo abelii","Macaca mulatta")

# write.table(binary_matrix_for_COUNT,file="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/binary_matrix_consistent_transcripts_COUNT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#####################################################################
### Inspect changes in splice sites (canonical <-> non_canonical) ###
#####################################################################

all_species_SQANTI<-droplevels(all_species_SQANTI)
binary_matrix_for_CANONICITY<-all_species_SQANTI %>% dplyr::select(isoform,species,all_canonical)
binary_matrix_for_CANONICITY$all_canonical<-ifelse(binary_matrix_for_CANONICITY$all_canonical=="canonical",1,0)
binary_matrix_for_CANONICITY<-reshape(binary_matrix_for_CANONICITY, idvar = "isoform", timevar = "species", direction = "wide")
colnames(binary_matrix_for_CANONICITY)<-c("family","Homo sapiens","Pan troglodytes","Gorilla gorilla","Pongo abelii","Macaca mulatta")
binary_matrix_for_COUNT<-binary_matrix_for_COUNT[order(binary_matrix_for_COUNT[,1],decreasing=FALSE),]
binary_matrix_for_CANONICITY<-binary_matrix_for_CANONICITY[order(binary_matrix_for_CANONICITY[,1],decreasing=FALSE),]

# Species-specific transcripts can arise from differences in junction canonicity or not. Separate both and re-run functional enrichment.
binary_matrix_for_CANONICITY_spsp_gains<-binary_matrix_for_CANONICITY[binary_matrix_for_CANONICITY$family %in% spsp_df$isoform,]
# From 1,215 innovations, 833 have conserved canonicity (832 are all canonical, 1 all non_canonical)
spsp_gains_conserved_canonicity<-as.vector(binary_matrix_for_CANONICITY_spsp_gains[(rowSums(binary_matrix_for_CANONICITY_spsp_gains[,2:6])==5) | (rowSums(binary_matrix_for_CANONICITY_spsp_gains[,2:6])==0),"family"])
# Get genes where species-specific transcripts arise from conserved canonicity.
spsp_gains_conserved_canonicity_GENES<-all_species_SQANTI %>% dplyr::filter(species=="Homo_sapiens" & isoform %in% spsp_gains_conserved_canonicity) %>% dplyr::select(associated_gene) %>% unique()
# The rest of transcripts (382) have changed canonicity in any species
spsp_gains_changed_canonicity<-as.vector(binary_matrix_for_CANONICITY_spsp_gains[!(rowSums(binary_matrix_for_CANONICITY_spsp_gains[,2:6])==5) & !(rowSums(binary_matrix_for_CANONICITY_spsp_gains[,2:6])==0),"family"])
# Get genes where species-specific transcripts show any change in splice junction canonicity
spsp_gains_changed_canonicity_GENES<-all_species_SQANTI %>% dplyr::filter(species=="Homo_sapiens" & isoform %in% spsp_gains_changed_canonicity) %>% dplyr::select(associated_gene) %>% unique()
# Among the 382, there are 184 that are all_canonical ONLY in the species where they are expressed
spsp_gains_only_canonical_where_they_expressed<-as.vector(binary_matrix_for_CANONICITY_spsp_gains[rowSums(binary_matrix_for_CANONICITY_spsp_gains[,2:6])==1,"family"])
# Get genes where species-specific transcripts show all canonical junctions only in the species expressing the spsp-transcript
spsp_gains_only_canonical_where_they_expressed_GENES<-all_species_SQANTI %>% dplyr::filter(species=="Homo_sapiens" & isoform %in% spsp_gains_only_canonical_where_they_expressed) %>% dplyr::select(associated_gene) %>% unique()

# Functional enrichments:
# 1) Spsp-gains with conserved canonicity
WebGestaltR(projectName="WebGestalt_spsp_gains_with_conserved_canonicity_GO",interestGene = as.vector(spsp_gains_conserved_canonicity_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")

WebGestaltR(projectName="WebGestalt_spsp_gains_with_conserved_canonicity_Panther",interestGene = as.vector(spsp_gains_conserved_canonicity_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

# 2) Spsp-gains with changed canonicity
WebGestaltR(projectName="WebGestalt_spsp_gains_with_changed_canonicity_GO",interestGene = as.vector(spsp_gains_changed_canonicity_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")

WebGestaltR(projectName="WebGestalt_spsp_gains_with_changed_canonicity_Panther",interestGene = as.vector(spsp_gains_changed_canonicity_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

# 3) Sp-sp gains that are only canonical in the species where they are expressed
WebGestaltR(projectName="WebGestalt_spsp_gains_only_canonical_where_they_expressed_GO",interestGene = as.vector(spsp_gains_only_canonical_where_they_expressed_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")

WebGestaltR(projectName="WebGestalt_spsp_gains_only_canonical_where_they_expressed_Panther",interestGene = as.vector(spsp_gains_only_canonical_where_they_expressed_GENES$associated_gene),referenceGene = as.vector(universe_genes),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

##########################################
### Calculate isoform fraction (usage) ###
##########################################

# IF is calculated for all transcript models, not only intra-sp consistent isoforms. TMM normalization is used.
# Re-calculate TPM
ISA_counts <- TMM_counts %>% dplyr::select(-GM12878.rep1)
ISA_eff_length <- eff_length %>% dplyr::select(-GM12878.rep1)
ISA_TPM<-tpm3(ISA_counts,ISA_eff_length)
ISA_TPM<-as.data.frame(ISA_TPM)
# Experimental design
myDesign <- data.frame(sampleID = colnames(ISA_counts),condition = c("Chimpanzee","Chimpanzee","Chimpanzee","Orangutan","Gorilla","Orangutan","Gorilla","Human","Human","Human","Gorilla","Orangutan","Macaque","Macaque","Macaque"),
batch=c("R","L","R","L","R","L","L","S","R","L","R","R","R","R","L"))
# Load SQANTI gff with all quantified transcripts
annotation<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/Homo_sapiens/QC/Homo_sapiens_5_sp_SQANTI3_corrected.gtf.cds.gff",sep="\t",header=FALSE,quote = "")
annotation$V10<-gsub(".+transcript_id \"(.+?)\";","\\1",annotation$V9)
# Keep quantified transcripts
annotation<-annotation %>% filter(V10 %in% rownames(ISA_TPM))
annotation<-annotation %>% dplyr::select(-V10)
# Write to file for ISA
write.table(annotation,file="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/Homo_sapiens_5_sp_SQANTI3_all_quantified_for_ISA.gtf",sep="\t",row.names = FALSE,quote=FALSE,col.names = FALSE)

# Load SwitchList object
SwitchList <- importRdata(
isoformRepExpression = ISA_TPM,
isoformCountMatrix   = ISA_counts,
designMatrix         = myDesign,
isoformExonAnnoation = "/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/Homo_sapiens_5_sp_SQANTI3_all_quantified_for_ISA.gtf",
isoformNtFasta       = "/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/Homo_sapiens/QC/Homo_sapiens_5_sp_SQANTI3_corrected.fasta",
showProgress = TRUE)

# Matrix with isoform fraction values (IF)
IF<-SwitchList$isoformRepIF
IF[is.na(IF)]<-0

# Plot PCA and HC of IF values
metadata<-data.frame(row.names=colnames(IF[,-1]))
metadata$samples<-gsub("(.+).rep[123]","\\1",rownames(metadata))
metadata$species<-c("Chimpanzee","Chimpanzee","Chimpanzee","Orangutan","Gorilla","Orangutan","Gorilla","Human","Human","Human","Gorilla","Orangutan","Macaque","Macaque","Macaque")
p <- pca(as.matrix(IF[,-1]), metadata = metadata, removeVar = 0.1)
biplot(p,legendPosition='right',colby = 'species',pointSize = 5,title='PCA for isoform usage (derived from TMM norm.TPM)',lab=p$metadata$samples) + theme_bw() + theme(aspect.ratio=0.9) + geom_point(colour="black",pch=21, size=5,alpha=0.5) + scale_fill_npg() + scale_color_npg() + theme(plot.title = element_text(hjust = 0.5,vjust=0.5)) + theme(legend.title = element_blank()) + coord_cartesian(xlim=c(-15,25))

my_spearman<-cor(IF[,-1],method="spearman")
t <- t(my_spearman)
rownames(t)<-gsub("(.+).rep[123]","\\1",rownames(t))
colnames(t)<-gsub("(.+).rep[123]","\\1",colnames(t))
dist.mat<-vegdist(t,method="euclidean")
clust.res<-hclust(dist.mat)
plot(clust.res)

# Inspect average IF in novel vs annotated transcripts
rownames(IF)<-IF$isoform_id
IF<-IF %>% dplyr::select(-isoform_id)
colnames(IF)<-c("Pan_troglodytes_CH114.rep3","Pan_troglodytes_CH170.rep1","Pan_troglodytes_CH391.rep3","Pongo_abelii_CRL-1850.rep1","Gorilla_gorilla_DIAN.rep3",
"Pongo_abelii_EB185.rep1","Gorilla_gorilla_GG05.rep1","Homo_sapiens_GM12878.rep2","Homo_sapiens_GM19150.rep3","Homo_sapiens_GM19238.rep1",
"Gorilla_gorilla_OMOYE.rep3","Pongo_abelii_PPY6_1.rep3","Macaca_mulatta_R02027.rep3","Macaca_mulatta_R05040.rep3","Macaca_mulatta_R94011.rep1")
colnames(IF)<-gsub('(.+_.+?)_.+', '\\1', colnames(IF))
IF<-as.data.frame ( sapply(split.default(IF, names(IF)), rowMeans) )
long_IF<-melt(as.matrix(IF))
colnames(long_IF)<-c("isoform","species","Average_IF")

# Add average IF to the SQANTI classification files
for (species in list_species) {

	# Load SQANTI df computed from each species and its genome/annotation
	SQANTI<-paste("/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/liftOver_merged_models/present_in_5_species/SQANTI_QC/",species,"/QC/",species,"_5_sp_SQANTI3_classification.txt",sep="")
	SQANTI<-read.table(file=SQANTI,sep="\t",header=TRUE)
	# Keep all quantified transcripts	
	SQANTI<-SQANTI %>% filter(isoform %in% rownames(ISA_TPM)) 
	SQANTI$associated_transcript2<-ifelse(SQANTI$associated_transcript=="novel","novel","annotated")
	SQANTI$species<-species
	SQANTI_name<-paste(species,"_SQANTI_all_quantified_transcripts",sep="")
	assign(x=SQANTI_name,value=SQANTI)
	
}

# rbind each species df into a single one
all_species_SQANTI_all_quantified_transcripts<-rbind(Homo_sapiens_SQANTI_all_quantified_transcripts,Pan_troglodytes_SQANTI_all_quantified_transcripts,Gorilla_gorilla_SQANTI_all_quantified_transcripts,Pongo_abelii_SQANTI_all_quantified_transcripts,Macaca_mulatta_SQANTI_all_quantified_transcripts)

# Merge by isoform+species to add Average IF
all_species_SQANTI_all_quantified_transcripts_IF<-merge(all_species_SQANTI_all_quantified_transcripts,long_IF,by=c("isoform","species"))

# Test for differences in average IF between annotated vs unannotated
wilcox.test(all_species_SQANTI_all_quantified_transcripts_IF$Average_IF~all_species_SQANTI_all_quantified_transcripts_IF$associated_transcript2)$p.value

##################################
### Differential isoform usage ### 
##################################

# First, remove genes with internal exon wobble (see Methods)
model_collapsing<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/liftOvers_new_assemblies/TAMA/merged_models_in_human_trans_report.txt",sep="\t",header=TRUE)
model_collapsing$all_wobble_list<-paste(model_collapsing$start_wobble_list,model_collapsing$end_wobble_list,sep=",")
model_collapsing<-model_collapsing %>% dplyr::select(transcript_id,all_wobble_list)
# Transform to long format
s <- strsplit(as.character(model_collapsing$all_wobble_list), split = ",")
model_collapsing_long_format <- data.frame(transcript_id = rep(model_collapsing$transcript_id, sapply(s, length)), wobble = unlist(s))
# Remove mono-exonic (only two wobble values corresponding to extremes)
model_collapsing_long_format<-model_collapsing_long_format %>% group_by(transcript_id) %>% dplyr::filter(n()>2)
# Filter wobble in extremes (start and beginning)
model_collapsing_long_format_only_internal_wobble<-model_collapsing_long_format %>% group_by(transcript_id) %>% dplyr::slice(2:n()) %>% dplyr::slice(1:(n()-1))
# Internal wobble to numeric
model_collapsing_long_format_only_internal_wobble$wobble<-as.numeric(as.character(model_collapsing_long_format_only_internal_wobble$wobble))
# Get isoform model IDs that were collapsed in internal exons
internal_exon_wobble_ids<-model_collapsing_long_format_only_internal_wobble %>% ungroup() %>% dplyr::filter(wobble>0) %>% dplyr::select(transcript_id) %>% unique()
# Get the corresponding genes. They will be excluded in this comparison (2,225 genes to be excluded)
internal_exon_wobble_genes<-complete_human_SQANTI_classification %>% dplyr::filter(isoform %in% internal_exon_wobble_ids$transcript_id) %>% dplyr::select(associated_gene) %>% unique()


### Isoform prefiltering ### 

# Remove 2,225 genes from SwitchList (7,524-2,225 genes) with merged internal exons.
# Author's recommendation for filtering genes (instead of modifying SwitchList$isofomFeatures df)
SwitchList<-subsetSwitchAnalyzeRlist(switchAnalyzeRlist = SwitchList,
subset = !SwitchList$isoformFeatures$gene_id %in% internal_exon_wobble_genes$associated_gene) # 5,299 evaluated genes

# Prefilter based on expression and usage (default)
SwitchListFiltered <- preFilter(switchAnalyzeRlist = SwitchList)

# Run differential isoform usage analysis
SwitchListFilteredAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchListFiltered,
    reduceToSwitchingGenes=TRUE)

##############################################
### Species-specific isoform usage changes ###
##############################################

library("pheatmap")
universe_IF_genes<-SwitchListFiltered$isoformFeatures %>% dplyr::select(gene_id) %>% unique() # universe of genes for DIU

# Define pairwise combinations
SwitchListFilteredAnalyzed$isoformFeatures$combo<-paste(SwitchListFilteredAnalyzed$isoformFeatures$condition_1,SwitchListFilteredAnalyzed$isoformFeatures$condition_2,sep="_")
# Get significant dIU changes by pairwise combination
sig<-SwitchListFilteredAnalyzed$isoformFeatures %>% filter(isoform_switch_q_value<0.05 & abs(dIF)>0.1) %>% dplyr::select(isoform_id,combo)
# List of lists
sig<-split(sig$isoform_id,sig$combo)

# 2 functions to retrieve data: 
# It will tell you the original rownames (default fromList() in UpSetR is not returning the names). 
# Taken from @docmanny https://github.com/hms-dbmi/UpSetR/issues/85
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
# Taken from @docmanny https://github.com/hms-dbmi/UpSetR/issues/85
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


# Retrieve transcripts showing species-specific usage changes
hs<-get_intersect_members(fromList(sig),'Chimpanzee_Human',
'Gorilla_Human',
'Human_Orangutan',
'Human_Macaque') 
cs<-get_intersect_members(fromList(sig),'Chimpanzee_Human',
'Chimpanzee_Gorilla',
'Chimpanzee_Orangutan',
'Chimpanzee_Macaque')
gs<-get_intersect_members(fromList(sig),'Gorilla_Human',
'Chimpanzee_Gorilla',
'Gorilla_Orangutan',
'Gorilla_Macaque')
os<-get_intersect_members(fromList(sig),'Human_Orangutan',
'Chimpanzee_Orangutan',
'Gorilla_Orangutan',
'Macaque_Orangutan')
ms<-get_intersect_members(fromList(sig),'Human_Macaque',
'Chimpanzee_Macaque',
'Gorilla_Macaque',
'Macaque_Orangutan')

# All species-specific isoform usage changes together
spsp_switches<-rbind(hs,cs,gs,os,ms)

# Generate barplot for species-specific DIU (novel vs annotated)
# Upregulated species-specific DIU
hs_up<-IF[rownames(IF) %in% rownames(hs),] %>% dplyr::filter(Homo_sapiens>Pan_troglodytes & Homo_sapiens>Gorilla_gorilla & Homo_sapiens>Pongo_abelii & Homo_sapiens>Macaca_mulatta)
cs_up<-IF[rownames(IF) %in% rownames(cs),] %>% dplyr::filter(Pan_troglodytes>Homo_sapiens & Pan_troglodytes>Gorilla_gorilla & Pan_troglodytes>Pongo_abelii & Pan_troglodytes>Macaca_mulatta)
gs_up<-IF[rownames(IF) %in% rownames(gs),] %>% dplyr::filter(Gorilla_gorilla>Homo_sapiens & Gorilla_gorilla>Pan_troglodytes & Gorilla_gorilla>Pongo_abelii & Gorilla_gorilla>Macaca_mulatta)
os_up<-IF[rownames(IF) %in% rownames(os),] %>% dplyr::filter(Pongo_abelii>Homo_sapiens & Pongo_abelii>Pan_troglodytes & Pongo_abelii>Gorilla_gorilla & Pongo_abelii>Macaca_mulatta)
ms_up<-IF[rownames(IF) %in% rownames(ms),] %>% dplyr::filter(Macaca_mulatta>Homo_sapiens & Macaca_mulatta>Pan_troglodytes & Macaca_mulatta>Gorilla_gorilla & Macaca_mulatta>Pongo_abelii)
# Downregulated species-specific DIU
hs_down<-IF[rownames(IF) %in% rownames(hs),] %>% dplyr::filter(Homo_sapiens<Pan_troglodytes & Homo_sapiens<Gorilla_gorilla & Homo_sapiens<Pongo_abelii & Homo_sapiens<Macaca_mulatta)
cs_down<-IF[rownames(IF) %in% rownames(cs),] %>% dplyr::filter(Pan_troglodytes<Homo_sapiens & Pan_troglodytes<Gorilla_gorilla & Pan_troglodytes<Pongo_abelii & Pan_troglodytes<Macaca_mulatta)
gs_down<-IF[rownames(IF) %in% rownames(gs),] %>% dplyr::filter(Gorilla_gorilla<Homo_sapiens & Gorilla_gorilla<Pan_troglodytes & Gorilla_gorilla<Pongo_abelii & Gorilla_gorilla<Macaca_mulatta)
os_down<-IF[rownames(IF) %in% rownames(os),] %>% dplyr::filter(Pongo_abelii<Homo_sapiens & Pongo_abelii<Pan_troglodytes & Pongo_abelii<Gorilla_gorilla & Pongo_abelii<Macaca_mulatta)
ms_down<-IF[rownames(IF) %in% rownames(ms),] %>% dplyr::filter(Macaca_mulatta<Homo_sapiens & Macaca_mulatta<Pan_troglodytes & Macaca_mulatta<Gorilla_gorilla & Macaca_mulatta<Pongo_abelii)
# As expected, they sum exactly the # of spsp-DIU

# Put all species-specific DIU in df with labels (up/down)
spsp_DIU_up_down<-rbind(data.frame(species="Homo_sapiens", direction="Upregulated",isoform=rownames(hs_up)),
data.frame(species="Homo_sapiens", direction="Downregulated",isoform=rownames(hs_down)),
data.frame(species="Pan_troglodytes", direction="Upregulated",isoform=rownames(cs_up)),
data.frame(species="Pan_troglodytes", direction="Downregulated",isoform=rownames(cs_down)),
data.frame(species="Gorilla_gorilla", direction="Upregulated",isoform=rownames(gs_up)),
data.frame(species="Gorilla_gorilla", direction="Downregulated",isoform=rownames(gs_down)),
data.frame(species="Pongo_abelii", direction="Upregulated",isoform=rownames(os_up)),
data.frame(species="Pongo_abelii", direction="Downregulated",isoform=rownames(os_down)),
data.frame(species="Macaca_mulatta", direction="Upregulated",isoform=rownames(ms_up)),
data.frame(species="Macaca_mulatta", direction="Downregulated",isoform=rownames(ms_down)))
# Add transcript annotation in that species
annotation_transcripts<-all_species_SQANTI_all_quantified_transcripts_IF %>% dplyr::select(isoform,species,associated_transcript2)
spsp_DIU_up_down<-merge(spsp_DIU_up_down,annotation_transcripts,by=c("isoform","species"))
# Formatting
colnames(spsp_DIU_up_down)<-c("isoform","Species","Direction","Annotation")
spsp_DIU_up_down$Annotation<-gsub("annotated","Known",spsp_DIU_up_down$Annotation)
spsp_DIU_up_down$Annotation<-gsub("novel","Novel",spsp_DIU_up_down$Annotation)
spsp_DIU_up_down$Species<-gsub("_"," ", spsp_DIU_up_down$Species)
spsp_DIU_up_down$Species <- factor(spsp_DIU_up_down$Species, levels = c("Macaca mulatta","Pongo abelii","Gorilla gorilla","Pan troglodytes","Homo sapiens"))
spsp_DIU_up_down$Annotation <- factor(spsp_DIU_up_down$Annotation, levels = c("Novel","Known"))
spsp_DIU_up_down$Direction <- factor(spsp_DIU_up_down$Direction, levels = c("Downregulated","Upregulated"))
# Barplot
spsp_DIU_up_down %>% group_by(Species,Direction,Annotation) %>% dplyr::count() %>% mutate(n = ifelse(Direction == "Upregulated",n,-1*n)) %>% ggplot(., aes(x=Species, y=n, alpha=Direction,fill=Annotation)) + geom_bar(stat="identity",width = .9,color="black") + scale_fill_manual(values=c("darkorange1", "dodgerblue4")) + scale_alpha_manual(values=c(0.2, 0.8)) + theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + scale_y_continuous(breaks = c(-300,0,300),labels = c(300,0,300)) + coord_flip(ylim=c(-300,300)) + ylab("Number of isoforms") + theme(aspect.ratio=0.6) +theme(text=element_text(size=21))

# Same barplot as above but now using Ensembl human reference annotation (known vs novel)
annotation_transcripts_human<-annotation_transcripts %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(isoform,associated_transcript2)
colnames(annotation_transcripts_human)<-c("isoform","Annotation_human")
spsp_DIU_up_down<-merge(spsp_DIU_up_down,annotation_transcripts_human)
spsp_DIU_up_down$Annotation_human<-gsub("annotated","Known",spsp_DIU_up_down$Annotation_human)
spsp_DIU_up_down$Annotation_human<-gsub("novel","Novel",spsp_DIU_up_down$Annotation_human)
spsp_DIU_up_down$Annotation_human <- factor(spsp_DIU_up_down$Annotation_human, levels = c("Novel","Known"))

spsp_DIU_up_down %>% group_by(Species,Direction,Annotation_human) %>% dplyr::count() %>% mutate(n = ifelse(Direction == "Upregulated",n,-1*n)) %>% ggplot(., aes(x=Species, y=n, alpha=Direction,fill=Annotation_human)) + geom_bar(stat="identity",width = .9,color="black") + scale_fill_manual(values=c("darkorange1", "dodgerblue4")) + scale_alpha_manual(values=c(0.2, 0.8)) + theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + scale_y_continuous(breaks = c(-300,0,300),labels = c(300,0,300)) + coord_flip(ylim=c(-300,300)) + ylab("Number of isoforms") + theme(aspect.ratio=0.6) +theme(text=element_text(size=21)) + guides(fill=guide_legend(order=1))

#######################################
### Functional enrichments spsp-DIU ###
#######################################

# Get genes showing species-specific isoform usage changes
spsp_switching_genes<-SwitchListFilteredAnalyzed$isoformFeatures[SwitchListFilteredAnalyzed$isoformFeatures$isoform_id %in% rownames(spsp_switches),] %>% dplyr::select(gene_id) %>% unique()

# Get genes showing species-specific isoform usage changes per species
hs_genes<-all_species_SQANTI_all_quantified_transcripts_IF[all_species_SQANTI_all_quantified_transcripts_IF$isoform %in% rownames(hs),] %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(associated_gene) %>% unique()
cs_genes<-all_species_SQANTI_all_quantified_transcripts_IF[all_species_SQANTI_all_quantified_transcripts_IF$isoform %in% rownames(cs),] %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(associated_gene) %>% unique()
gs_genes<-all_species_SQANTI_all_quantified_transcripts_IF[all_species_SQANTI_all_quantified_transcripts_IF$isoform %in% rownames(gs),] %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(associated_gene) %>% unique()
os_genes<-all_species_SQANTI_all_quantified_transcripts_IF[all_species_SQANTI_all_quantified_transcripts_IF$isoform %in% rownames(os),] %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(associated_gene) %>% unique()
ms_genes<-all_species_SQANTI_all_quantified_transcripts_IF[all_species_SQANTI_all_quantified_transcripts_IF$isoform %in% rownames(ms),] %>% dplyr::filter(species=="Homo_sapiens") %>% dplyr::select(associated_gene) %>% unique()

# Functional enrichment
WebGestaltR(projectName="WebGestalt_spsp_DIU_GO",interestGene = as.vector(spsp_switching_genes$gene_id),referenceGene = as.vector(universe_IF_genes$gene_id),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="geneontology_Biological_Process",organism="hsapiens")
WebGestaltR(projectName="WebGestalt_spsp_DIU_Panther",interestGene = as.vector(spsp_switching_genes$gene_id),referenceGene = as.vector(universe_IF_genes$gene_id),interestGeneType="ensembl_gene_id",referenceGeneType="ensembl_gene_id",enrichDatabase="pathway_Panther",organism="hsapiens")

# No significant enrichment for spsp-DIU (all together or separated by species)

###########################################################
### Gene-level classification for evolutionary analyses ###
###########################################################

# Generate isoform and gene classification based on transcript expression patterns (innovations), DEU and DIU

###############
# Innovations #
###############

# Group 1: genes where all isoforms are either expressed or non-expressed in all species (all_conserved)
# Group 2: genes where all isoforms are either human innovation (single species), conserved or other (only_human_sp)
# Group 3: genes where all isoforms are either NHP innovation (single species), conserved or other (only_NHP_sp)
# Group 4: genes with isoforms that are spsp gains in several species, conserved or other (convergence_spsp)
# Group 5: genes with 'other' or 'conserved, other'

# We start from isoforms that are intra-sp consistent
colnames(binary_matrix_for_COUNT)<-c("family","Homo_sapiens","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Macaca_mulatta")
# Classify the behavior of each isoform
classify_innovations<-binary_matrix_for_COUNT %>% dplyr::mutate(
category = case_when(Homo_sapiens=="1" & Pan_troglodytes=="0" & Gorilla_gorilla=="0" & Pongo_abelii=="0" & Macaca_mulatta=="0" ~ "human_sp",
Homo_sapiens=="0" & Pan_troglodytes=="1" & Gorilla_gorilla=="0" & Pongo_abelii=="0" & Macaca_mulatta=="0" ~ "chimpanzee_sp",
Homo_sapiens=="0" & Pan_troglodytes=="0" & Gorilla_gorilla=="1" & Pongo_abelii=="0" & Macaca_mulatta=="0" ~ "gorilla_sp",
Homo_sapiens=="0" & Pan_troglodytes=="0" & Gorilla_gorilla=="0" & Pongo_abelii=="1" & Macaca_mulatta=="0" ~ "orangutan_sp",
Homo_sapiens=="0" & Pan_troglodytes=="0" & Gorilla_gorilla=="0" & Pongo_abelii=="0" & Macaca_mulatta=="1" ~ "macaque_sp",
Homo_sapiens=="1" & Pan_troglodytes=="1" & Gorilla_gorilla=="1" & Pongo_abelii=="1" & Macaca_mulatta=="1" ~ "conserved",
Homo_sapiens=="0" & Pan_troglodytes=="0" & Gorilla_gorilla=="0" & Pongo_abelii=="0" & Macaca_mulatta=="0" ~ "conserved",
TRUE ~ "other"))
# Add gene
isoform_2_gene<-Homo_sapiens_SQANTI %>% dplyr::select(isoform,associated_gene)
classify_innovations<-merge(classify_innovations,isoform_2_gene,by.x="family",by.y="isoform")
# All classes of isoforms per gene (alphabetical sorting)
classify_innovations<-as.data.table(classify_innovations)[, toString(sort(unique(category))), by = "associated_gene"]
colnames(classify_innovations)[2]<-"isoform_class"
# Provide a gene class depending on its isoforms
classify_innovations<-classify_innovations %>% dplyr::mutate(
gene_class = case_when(isoform_class=="conserved" ~ "all_conserved", 
isoform_class=="human_sp" | isoform_class=="conserved, human_sp" | isoform_class=="conserved, human_sp, other" | isoform_class=="human_sp, other"  ~ "only_human_sp",
isoform_class=="chimpanzee_sp" | isoform_class=="chimpanzee_sp, conserved" | isoform_class=="chimpanzee_sp, conserved, other" | isoform_class=="chimpanzee_sp, other" ~ "only_NHP_sp",
isoform_class=="gorilla_sp" | isoform_class=="conserved, gorilla_sp" | isoform_class=="conserved, gorilla_sp, other" | isoform_class=="gorilla_sp, other" ~ "only_NHP_sp",
isoform_class=="orangutan_sp" | isoform_class=="conserved, orangutan_sp" | isoform_class=="conserved, orangutan_sp, other" | isoform_class=="orangutan_sp, other"  ~ "only_NHP_sp",
isoform_class=="macaque_sp" | isoform_class=="conserved, macaque_sp" | isoform_class=="conserved, macaque_sp, other" | isoform_class=="macaque_sp, other" ~ "only_NHP_sp",
isoform_class=="other" | isoform_class=="conserved, other" ~ "other",
TRUE ~ "convergence_spsp"
))


# Plot human population-based dN/dS from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737536/ (Supp. File 17)
dN_dS_human<-read.table(file="/media/luis/Data/PhD/Isoseq_LCLs/dNdS_primates/ExAc_human_dNdS.txt",sep="\t",header=TRUE)
dN_dS_human<-dN_dS_human %>% dplyr::select(Ensembl.Gene.ID,PSG,ExAC.unique_variants.dN_dS)

# dN/dS in our groups
classify_innovations_with_dN_dS<-merge(classify_innovations,dN_dS_human,by.x="associated_gene",by.y="Ensembl.Gene.ID")
classify_innovations_with_dN_dS$gene_class <- factor(classify_innovations_with_dN_dS$gene_class, levels = c("all_conserved", "only_human_sp", "only_NHP_sp","other","convergence_spsp"))
ggboxplot(classify_innovations_with_dN_dS,x="gene_class",y="ExAC.unique_variants.dN_dS",fill="gene_class") + coord_cartesian(ylim=c(0,1)) + geom_violin(alpha=0.1)

#########
## DIU ##
#########

# Group 1: genes where all isoforms showed no change (all_conserved)
# Group 2: genes where all isoforms are either human up (single species) or conserved or other (only_human_up)
# Group 3: genes where all isoforms are either NHP up (single species) or conserved or other (only_NHP_up)
# Group 4: genes with isoforms that are spsp-DIU-up in several species or conserved or other (convergence_spsp)
# Group 5: genes with 'other' patterns of change (3+2 or spsp-down) or 'conserved, other'

# We want to classify isoforms in all evaluated genes (present in SwitchListFiltered)
classify_DIU<-SwitchListFiltered$isoformFeatures %>% dplyr::select(isoform_id) %>% unique()

# Select 3+2 cases to be included as 'Other' and ignore the rest of pattern changes. There are 10 combinations 3+2. Select by hand.
h_c<-get_intersect_members(fromList(sig),'Gorilla_Human','Human_Orangutan','Human_Macaque','Chimpanzee_Gorilla','Chimpanzee_Orangutan','Chimpanzee_Macaque') 
h_g<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Human_Orangutan','Human_Macaque','Chimpanzee_Gorilla','Gorilla_Orangutan','Gorilla_Macaque') 
h_o<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Gorilla_Human','Human_Macaque','Chimpanzee_Orangutan','Gorilla_Orangutan','Macaque_Orangutan') 
h_m<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Gorilla_Human','Human_Orangutan','Chimpanzee_Macaque','Gorilla_Macaque','Macaque_Orangutan') 
c_g<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Chimpanzee_Orangutan','Chimpanzee_Macaque','Gorilla_Human','Gorilla_Orangutan','Gorilla_Macaque') 
c_o<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Chimpanzee_Gorilla','Chimpanzee_Macaque','Human_Orangutan','Gorilla_Orangutan','Macaque_Orangutan') 
c_m<-get_intersect_members(fromList(sig),'Chimpanzee_Human','Chimpanzee_Gorilla','Chimpanzee_Orangutan','Human_Macaque','Gorilla_Macaque','Macaque_Orangutan') 
g_o<-get_intersect_members(fromList(sig),'Gorilla_Human','Chimpanzee_Gorilla','Gorilla_Macaque','Human_Orangutan','Chimpanzee_Orangutan','Macaque_Orangutan')
g_m<-get_intersect_members(fromList(sig),'Gorilla_Human','Chimpanzee_Gorilla','Gorilla_Orangutan','Human_Macaque','Chimpanzee_Macaque','Macaque_Orangutan')
o_m<-get_intersect_members(fromList(sig),'Human_Orangutan','Chimpanzee_Orangutan','Gorilla_Orangutan','Human_Macaque','Chimpanzee_Macaque','Gorilla_Macaque')

# Groups can be checked using upset and pheatmap
#ok<-fromList(sig)
#ok<-ok[rowSums(ok)=="6",]
#upset(ok,nsets=10,nintersects=1000,order.by="freq")
#IF[rownames(IF) %in% rownames(o_m),] %>% pheatmap(.) # see HC grouping

# Group 'Other' will be the union of 3+2 groups and spsp down isoforms. 436 isoforms in total
other<-c(rownames(h_c),rownames(h_g),rownames(h_o),rownames(h_m),rownames(c_g),rownames(c_o),rownames(c_m),rownames(g_o),rownames(g_m),rownames(o_m),rownames(hs_down),rownames(cs_down),rownames(gs_down),rownames(os_down),rownames(ms_down))

# Classify isoforms in DIU
classify_DIU<-classify_DIU %>% dplyr::mutate(
category = case_when(isoform_id %in% rownames(hs_up) ~ "human_up",
isoform_id %in% rownames(cs_up) ~ "chimpanzee_up",
isoform_id %in% rownames(gs_up) ~ "gorilla_up",
isoform_id %in% rownames(os_up) ~ "orangutan_up",
isoform_id %in% rownames(ms_up) ~ "macaque_up",
!isoform_id %in% rownames(fromList(sig)) ~ "conserved",
isoform_id %in% other ~ "other")) 

# This will generate 'NA' corresponding to isoforms that are not conserved, spsp DIU up, or other (spsp DIU down or any other 3+2 change). 
# We will ignore these 'NA' here. They correspond to other changes not consistent with differences between groups of species (see Methods).
classify_DIU<-classify_DIU[!is.na(classify_DIU$category),]

# Add gene
isoform_2_gene<-Homo_sapiens_SQANTI_all_quantified_transcripts %>% dplyr::select(isoform,associated_gene)
classify_DIU<-merge(classify_DIU,isoform_2_gene,by.x="isoform_id",by.y="isoform")

# Classify the type of gene
classify_DIU<-as.data.table(classify_DIU)[, toString(sort(unique(category))), by = "associated_gene"]
colnames(classify_DIU)[2]<-"isoform_class"

# Provide a gene class depending on DIU of its isoforms
classify_DIU<-classify_DIU %>% dplyr::mutate(
gene_class = case_when(isoform_class=="conserved" ~ "all_conserved", 
isoform_class=="human_up" | isoform_class=="conserved, human_up" | isoform_class=="conserved, human_up, other" | isoform_class=="human_up, other" ~ "only_human_up",
isoform_class=="chimpanzee_up" | isoform_class=="chimpanzee_up, conserved" | isoform_class=="chimpanzee_up, conserved, other" | isoform_class=="chimpanzee_up, other" ~ "only_NHP_up",
isoform_class=="gorilla_up" | isoform_class=="conserved, gorilla_up" | isoform_class=="conserved, gorilla_up, other" | isoform_class=="gorilla_up, other" ~ "only_NHP_up",
isoform_class=="orangutan_up" | isoform_class=="conserved, orangutan_up" | isoform_class=="conserved, orangutan_up, other" | isoform_class=="orangutan_up, other" ~ "only_NHP_up",
isoform_class=="macaque_up" | isoform_class=="conserved, macaque_up" | isoform_class=="conserved, macaque_up, other" | isoform_class=="macaque_up, other"  ~ "only_NHP_up",
isoform_class=="other" | isoform_class=="conserved, other" ~ "other",
TRUE ~ "convergence_spsp"
))

# dN/dS in our groups
classify_DIU_with_dN_dS<-merge(classify_DIU,dN_dS_human,by.x="associated_gene",by.y="Ensembl.Gene.ID")
classify_DIU_with_dN_dS$gene_class <- factor(classify_DIU_with_dN_dS$gene_class, levels = c("all_conserved", "only_human_up", "only_NHP_up","other","convergence_spsp"))
ggboxplot(classify_DIU_with_dN_dS,x="gene_class",y="ExAC.unique_variants.dN_dS",fill="gene_class") + coord_cartesian(ylim=c(0,1)) + geom_violin(alpha=0.1)

# The gene classes defined in 'classify_innovations' (based on transcript expression gains/losses) and 'classify_DIU' (based on DIU)
# are combined together with gene classes for differential exon usage (DEU) to generate a combined classification of genes taking into account the 3 analyses.



