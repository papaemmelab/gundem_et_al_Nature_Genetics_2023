library('NMF')
library('GenomicRanges')
library('GenomeInfoDb')
library('plotrix')
library("ggplot2")
library("BSgenome.Hsapiens.UCSC.hg19")
library("reshape2")
library('knitr')
library("plotly")
library("ape")
library("ggdendro")
library("deconstructSigs")
library('MutationalPatterns')
library(ggpubr)
library(cowplot)

WDIR = '~/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source(file="scripts/util_signatures.R")
source('scripts/util_dataset.R')

# Read reference signtures from COSMIC and Kucab et al.
roworder = c(
	"A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
	"A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
	"A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
	"A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
	"A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
	"A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
)
cosmic_file = 'data/cosmic_signatures.csv'
kucab_file = 'data/Mutagen53_sub_signature.txt'
cosmic_signatures = read.csv(file=cosmic_file, header=TRUE)
cosmic_signatures = as.matrix(cosmic_signatures[,3:ncol(cosmic_signatures)])
mutagen_signatures= read.delim(file=kucab_file, header=T)
cancer_signatures=cbind(cosmic_signatures,mutagen_signatures$Temozolomide..200.uM.)
colnames(cancer_signatures)[66]="TMZ"
rownames(cancer_signatures) = roworder
# Read sample-level mutation counts for 96 mutation types
sample_level_mut_mat = read.table('data/sample_level_96_channel_mutations_counts.txt', header=T, sep="\t", stringsAsFactors=F)
colnames(sample_level_mut_mat) = gsub('\\.','-',colnames(sample_level_mut_mat))
sample_channels = t(sample_level_mut_mat)
sample_channels = sample_channels / rowSums(sample_channels)
# Read de novo extracted signatures
set_of_signatures = read.table('data/denovo_signatures.txt', row.names=1, sep="\t", header=T, stringsAsFactors=F)
# Read indel exposures
indel_act_file = 'data/sample_level_signature_exposures_indels.txt'
indel_exposures = read.table(indel_act_file, header=T, sep="\t", stringsAsFactors=F)
indel_exposures = indel_exposures[match(rownames(sample_channels), indel_exposures$Samples),]
indel_exposures$Sum = rowSums(indel_exposures[,colnames(indel_exposures)!="Samples"])
# Only keep the signatures which contribute to >1% of the total mutation burden across the cohort
indel_signs = setdiff(colnames(indel_exposures), c('Samples','Sum'))
selected_indel_signs = c()
for(sign in indel_signs) {
	total_contribution = sum(indel_exposures[,sign]) / sum(indel_exposures$Sum)
	if(total_contribution > 0.01) {selected_indel_signs = c(selected_indel_signs, sign)}
}
indel_exposures = indel_exposures[,c('Samples',selected_indel_signs)]
indel_exposures$Sum = rowSums(indel_exposures[,2:ncol(indel_exposures)])
indel_signs = setdiff(colnames(indel_exposures), c('Samples','Sum'))
for(sign in indel_signs) {
	indel_exposures[,sign] = indel_exposures[,sign] / indel_exposures$Sum
}
# Sample annotations
sample_info = mdata
sample_info$DISEASE_SUBTYPE[sample_info$DISEASE_SUBTYPE=='11q'] = 'SEG_CNV'
sample_info$DISEASE_SUBTYPE[sample_info$DISEASE_SUBTYPE=='NUM_CNV'] = 'NUMERIC_CNV'
# Extended Data Figure 2a
plot_ref = plot_96_profile(cancer_signatures[,c('SBS18','SBS38','SBS40','TMZ','SBS31','SBS35')])
plot_denovo = plot_96_profile(set_of_signatures)
pdf('raw/ext_data_fig_2a.pdf', width=12)
plot_grid(plot_denovo, plot_ref)
dev.off()
# Signature fitting
sieve = perso_fit_to_signatures(as.matrix(sample_level_mut_mat), as.matrix(set_of_signatures), ind_prior_knowledge=c(1,2), n_ite=10, initial_cutoff=0.4, cutoff=0.01, cutoff_add=0.35)
rownames(sieve)=colnames(set_of_signatures)
colnames(sieve)=colnames(sample_level_mut_mat)
# Calculate exposures
sieve[which(sieve>0)] = 1
NBL_dataset_activities = sieve_to_exposure(samples=as.matrix(sample_level_mut_mat), signatures=as.matrix(set_of_signatures), sieve=sieve, threshold=0.85)
sample_exposures = as.data.frame(t(NBL_dataset_activities))
sample_exposures$sample = rownames(sample_exposures)
rownames(sample_exposures) = NULL
sample_exposures = cbind(sample_exposures[,'sample'], sample_exposures[,setdiff(colnames(sample_exposures), 'sample')])
colnames(sample_exposures)[1] = 'sample'
write.table(sample_exposures, file='data/sample_level_signature_exposures.txt', row.names=F, col.names=T, sep="\t", quote=F)
sample_exposures$total_subs = rowSums(sample_exposures[,setdiff(colnames(sample_exposures), 'sample')])
sample_exposures$A_SBS18 = sample_exposures$A_SBS18 / sample_exposures$total_subs
sample_exposures$B_SBS38 = sample_exposures$B_SBS38 / sample_exposures$total_subs
sample_exposures$C_SBS40 = sample_exposures$C_SBS40 / sample_exposures$total_subs
sample_exposures$D_TMZ = sample_exposures$D_TMZ / sample_exposures$total_subs
sample_exposures$E_SBS31 = sample_exposures$E_SBS31 / sample_exposures$total_subs
sample_exposures$F_SBS35 = sample_exposures$F_SBS35 / sample_exposures$total_subs
sample_exposures = merge(x=sample_exposures, y=sample_info[sample_info$SS_ANALYSIS=='yes',c('DNA_SAMPLE','SAMPLE_TYPE','DISEASE_SUBTYPE')], by.x='sample', by.y='DNA_SAMPLE')
sample_exposures$SAMPLE_TYPE = factor(sample_exposures$SAMPLE_TYPE, levels=sample_type_cols$sample_type)
sample_exposures$DISEASE_SUBTYPE = factor(sample_exposures$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
# Extended Data Figure 2b
library(ComplexHeatmap)
library(circlize)
sample_channels = t(sample_level_mut_mat)
rownames(sample_channels) = gsub('\\.','-',rownames(sample_channels))
sample_channels = sample_channels / rowSums(sample_channels)
dplot = sample_channels
colnames(dplot) = NULL
rownames(dplot) = NULL
row_annotation = sample_info[match(rownames(sample_channels), sample_info$DNA_SAMPLE),c('DNA_SAMPLE','DISEASE_SUBTYPE','SAMPLE_TYPE','PLT_THERAPY','TMZ_THERAPY','XRT_TO_TUMOR','NU_SUBS','NU_INDELS')]
exposure_col_fun = colorRamp2(c(0,1), c("white", "brown"))
right_ha = HeatmapAnnotation(
	signature_exposure=cbind(A_SBS18=sample_exposures$A_SBS18, B_SBS38=sample_exposures$B_SBS38, C_SBS40=sample_exposures$C_SBS40, D_TMZ=sample_exposures$D_TMZ, E_SBS31=sample_exposures$E_SBS31, F_SBS35=sample_exposures$F_SBS35),
	indel_exposure=cbind(ID1=indel_exposures$ID1, ID2=indel_exposures$ID2, ID3=indel_exposures$ID3, ID5=indel_exposures$ID5, ID9=indel_exposures$ID9, ID10=indel_exposures$ID10),
	No_subs=anno_barplot(as.matrix(row_annotation$NU_SUBS)), No_indels=anno_barplot(as.matrix(row_annotation$NU_INDELS)), which='row',
	col=list(signature_exposure=exposure_col_fun, indel_exposure=exposure_col_fun), show_legend=c("indel_exposure"=FALSE)
)
top_ha = HeatmapAnnotation(
	mut_type=c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16), rep("T>A", 16), rep("T>C", 16), rep("T>G", 16)),
	col=list(mut_type=c("C>A"="royalblue","C>G"="black","C>T"="red","T>A"="grey","T>C"="green2","T>G"="hotpink"))
)
left_ha = rowAnnotation(
	recurrent_lesion=row_annotation$DISEASE_SUBTYPE, sample_type=row_annotation$SAMPLE_TYPE, CDDP_therapy=row_annotation$PLT_THERAPY,
	TMZ_therapy=row_annotation$TMZ_THERAPY, Radiotherapy=row_annotation$XRT_TO_TUMOR,
	col=list(
		recurrent_lesion=c(
			"MDM2-CDK4"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='MDM2-CDK4','color'],
			"MYCN"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='MYCN','color'],
			"NUMERIC_CNV"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='NUMERIC_CNV','color'],
			"SEG_CNV"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='SEG_CNV','color'],
			"TERT"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='TERT','color'],
			"ATRX"=disease_subtype_cols[disease_subtype_cols$disease_subtype=='ATRX','color']#,
		), sample_type=c(
			"diagnosis"=sample_type_cols[sample_type_cols$sample_type=='diagnosis','color'],
			"reresection"=sample_type_cols[sample_type_cols$sample_type=='reresection','color'],
			"t-resection"=sample_type_cols[sample_type_cols$sample_type=='2nd_look','color'],
			"relapse"=sample_type_cols[sample_type_cols$sample_type=='relapse','color'],
			"frelapse"=sample_type_cols[sample_type_cols$sample_type=='frelapse','color']
		), CDDP_therapy=c("no"="grey80", "yes"="grey30"),
		TMZ_therapy=c("no"="grey80", "yes"="grey30"),
		Radiotherapy=c("no_xrt"="grey80", "xrt_direct"="grey30", "xrt_other"="grey50")
	)
)
hm <- Heatmap(dplot, cluster_columns=FALSE, col=colorRamp2(c(0, 0.04), c("white", "darkgrey")), top_annotation=top_ha, left_annotation=left_ha, right_annotation=right_ha)
pdf('raw/ext_data_fig_2b.pdf', width=12)
draw(hm)
dev.off()
