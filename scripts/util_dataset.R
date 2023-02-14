library(ggpubr)
library(cowplot)
library(plyr)
library(reshape2)
library(gtable)
library(grid)
library(gridExtra)

WDIR = "/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final"
setwd(WDIR)
# Colors for disease subtypes
disease_subtype_cols = as.data.frame(rbind(
	c("MYCN","#009E73"),c("TERT","#9ad0f3"),c('ATRX',"grey45"),c("MDM2-CDK4","#CC79A7"), c("SEG_CNV","#0072B2"), c("NUMERIC_CNV","#e79f00")
))
colnames(disease_subtype_cols) = c('disease_subtype','color')
disease_subtype_cols$disease_subtype = as.character(disease_subtype_cols$disease_subtype)
disease_subtype_cols$color = as.character(disease_subtype_cols$color)
# Colors for sample types
sample_type_cols = as.data.frame(cbind(
	c('diagnosis','reresection','2nd_look','relapse','frelapse'),
	c("#4A708B","#87CEEB","#6E8B3D","#CD8500","#CD5C5C")
))
colnames(sample_type_cols) = c('sample_type','color')
sample_type_cols$sample_type = as.character(sample_type_cols$sample_type)
sample_type_cols$color = as.character(sample_type_cols$color)
# Colors for tumor site
tumor_site_colors = as.data.frame(cbind(c('primary_disease','local_disease','distant_metastasis'), c('#4A708B','#AA4400','#800000')))
colnames(tumor_site_colors) = c('tumor_site','color')
tumor_site_colors$tumor_site = as.character(tumor_site_colors$tumor_site)
tumor_site_colors$color = as.character(tumor_site_colors$color)
# Exposures for substitution signatures
abs_exposures = read.table('data/sample_level_signature_exposures.txt', header=T, sep="\t", stringsAsFactors=F)
subs_exposures = abs_exposures[,colnames(abs_exposures)!='sample'] / rowSums(abs_exposures[,colnames(abs_exposures)!='sample'])
subs_exposures$sample = abs_exposures$sample
colnames(abs_exposures)[grep('sample',colnames(abs_exposures),invert=T)] = paste(grep('sample',colnames(abs_exposures),invert=T,value=T), '_n', sep="")
# Read the cohort information
mdata = read.delim('data/final_cohort.txt', sep="\t", stringsAsFactors=F)
mdata$PATIENT = sub('E-H-','H',sub('I-H-','H',mdata$PATIENT))
mdata = mdata[order(mdata$PATIENT),]
mdata$DISEASE_SUBTYPE[mdata$patient %in% c(112909, 134817, 134822)] = 'ATRX'
mdata$DISEASE_SUBTYPE[mdata$patient %in% c(112910, 132376, 132379)] = 'TERT'
mdata$AGE_GROUP = 1
mdata$AGE_GROUP[which(mdata$AADx>=18 & mdata$AADx<(12*12))] = 2
mdata$AGE_GROUP[which(mdata$AADx>=(12*12))] = 3
mdata$AGE_GROUP = paste('A', mdata$AGE_GROUP, sep='')
mdata$STAGE[mdata$STAGE %in% c('2A','2B')] = 2
mdata = mdata[which(mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','2ndlook','relapse','frelapse','autopsy')),]
mdata$no_tumor_samples = NULL
patients_with_multiple_tumors = as.data.frame(table(mdata$PATIENT[which(mdata$MS_ANALYSIS=='yes')]))
colnames(patients_with_multiple_tumors) = c('patient','no_tumor_samples')
mdata = merge(x=mdata, y=patients_with_multiple_tumors, by.x='PATIENT', by.y='patient', all.x=T)
mdata$no_tumor_samples[is.na(mdata$no_tumor_samples)] = 1
mdata = merge(x=mdata, y=subs_exposures, by.x='DNA_SAMPLE', by.y='sample', all.x=T)
mdata = merge(x=mdata, y=abs_exposures, by.x='DNA_SAMPLE', by.y='sample', all.x=T)
mdata$SAMPLE_TYPE[mdata$SAMPLE_TYPE=='2ndlook'] = '2nd_look'
mdata$SAMPLE_TYPE = factor(mdata$SAMPLE_TYPE, levels=sample_type_cols$sample_type)
mdata$DISEASE_SUBTYPE = factor(mdata$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
mdata$AGE_GROUP2 = mdata$AGE_GROUP
mdata$AGE_GROUP2[mdata$AGE_GROUP2=='A1'] = '<18m'
mdata$AGE_GROUP2[mdata$AGE_GROUP2=='A2'] = '>=18m & <12y'
mdata$AGE_GROUP2[mdata$AGE_GROUP2=='A3'] = '>=12y'
mdata$AGE_GROUP2 = factor(mdata$AGE_GROUP2, levels=c('<18m','>=18m & <12y','>=12y'))
mdata$STAGE = factor(mdata$STAGE, levels=c('1','2','3','4N','4'))
# Read all mutations for WGS samples
muts_wgs = read.table('data/wgs_mutations_annotated.txt', header=T, sep="\t", stringsAsFactors=F)
muts_wgs$CHANGE = NA
muts_wgs$CHANGE[muts_wgs$VAG_VT=='Sub'] = apply(muts_wgs[muts_wgs$VAG_VT=='Sub',c('REF','ALT')], 1, paste, collapse=">")
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='A>T')] = 'T>A'
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='A>G')] = 'T>C'
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='G>T')] = 'C>A'
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='G>C')] = 'C>G'
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='G>A')] = 'C>T'
muts_wgs$CHANGE[which(muts_wgs$CHANGE=='A>C')] = 'T>G'
muts_wgs$depth = muts_wgs$wt.count + muts_wgs$mut.count
muts_wgs$vaf = muts_wgs$mut.count / muts_wgs$depth
wgs_selected_columns = c(
	'TARGET_NAME','POS_KEY','VAG_VT','VAG_EFFECT','VAG_GENE','VAG_PROTEIN_CHANGE','CHANGE','PATHWAY','CHR','START','REF','ALT','vaf','depth',
	'CHR1','START1','CHR2','START2','STRAND1','STRAND2','ccf','cluster','cluster_descr','subclone_type'
)
oncogenic_wgs = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')],
	y=muts_wgs[muts_wgs$DRIVER=='yes',wgs_selected_columns],
	by.x='DNA_SAMPLE', by.y='TARGET_NAME'
)
all_wgs = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')], y=muts_wgs[,c(wgs_selected_columns,'DRIVER')], by.x='DNA_SAMPLE', by.y='TARGET_NAME'
)
# Read all mutations for IMPACT samples
muts_impact = read.delim('data/mskimpact_coding_mutations_annotated.txt', header=T, sep="\t", stringsAsFactors=F)
muts_impact$POS_KEY = paste(muts_impact$Chromosome, muts_impact$Start_Position, sep="_")
muts_impact$change = NA
sel = muts_impact$Variant_Type=='SNP' & !is.na(muts_impact$Tumor_Seq_Allele1) & !is.na(muts_impact$Tumor_Seq_Allele2)
muts_impact$change[sel] = apply(muts_impact[sel,c('Tumor_Seq_Allele1','Tumor_Seq_Allele2')], 1, paste, collapse=">")
muts_impact$change[which(muts_impact$change=='A>T')] = 'T>A'
muts_impact$change[which(muts_impact$change=='A>G')] = 'T>C'
muts_impact$change[which(muts_impact$change=='G>T')] = 'C>A'
muts_impact$change[which(muts_impact$change=='G>C')] = 'C>G'
muts_impact$change[which(muts_impact$change=='G>A')] = 'C>T'
muts_impact$change[which(muts_impact$change=='A>C')] = 'T>G'
muts_impact$t_depth = muts_impact$t_wt_count + muts_impact$t_alt_count
#Oncogenic mutations
sel = which(is.na(muts_impact$isabl_patient) & muts_impact$Oncogenic==1 & muts_impact$Mutation_Status=='SOMATIC')
impact_selected_columns = c(
	'Tumor_Sample_Barcode','POS_KEY','Variant_Type','Consequence','Hugo_Symbol','HGVSp_Short','change','Pathway','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2',
	't_vaf','t_depth','Chr1','Start1','Chr2','Start2','Strand1','Strand2','CCF','clone_id','clone_descr','subclone_type'
)
oncogenic_impact = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')],
	y=muts_impact[sel,impact_selected_columns],
	by.x='DNA_SAMPLE', by.y='Tumor_Sample_Barcode'
)
all_impact = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')],
	y=muts_impact[which(is.na(muts_impact$isabl_patient) & muts_impact$Mutation_Status=='SOMATIC'),c(impact_selected_columns, 'Oncogenic')],
	by.x='DNA_SAMPLE', by.y='Tumor_Sample_Barcode'
)
colnames(all_impact) = colnames(all_wgs)
all_mutations = rbind(all_wgs, all_impact)
colnames(oncogenic_impact) = colnames(oncogenic_wgs)
oncogenic = rbind(oncogenic_impact, oncogenic_wgs)
oncogenic$VAG_VT[oncogenic$VAG_VT=='DEL'] = 'Del'
oncogenic$VAG_VT[oncogenic$VAG_VT=='SNP'] = 'Sub'
oncogenic$VAG_VT[oncogenic$VAG_VT=='DNP'] = 'Sub'
oncogenic$VAG_VT[oncogenic$VAG_VT=='INS'] = 'Ins'
oncogenic$VAG_VT[oncogenic$VAG_VT=='ONP'] = 'Sub'
oncogenic$VAG_EFFECT[oncogenic$VAG_EFFECT=="frameshit_variant"] = "frameshift_variant"
oncogenic$VAG_EFFECT[oncogenic$VAG_EFFECT=="Nonsense"] = "stop_gained"
oncogenic$VAG_EFFECT[oncogenic$VAG_EFFECT=="frameshift_variant,splice_region_variant"] = "frameshift_variant"
oncogenic$VAG_EFFECT[oncogenic$VAG_EFFECT %in% c("missense_variant",'missense','missense_variant,splice_region_variant')] = "non_synonymous_codon"
oncogenic$VAG_EFFECT[oncogenic$VAG_EFFECT %in% c('splice_acceptor_variant','splice_donor_variant','splice_acceptor_variant,intron_variant')] = "splice_site_variant"
oncogenic$VAG_GENE[oncogenic$VAG_GENE=='IGF2BP3-MDM2-CDK4'] = 'IGF2BP3'
oncogenic = unique(oncogenic)
oncogenic$cohort = 'MSK-WGS'
oncogenic$cohort[grepl('^P',oncogenic$PATIENT)] = 'MSK-IMPACT'
# All mutations included in clonality analysis
clonality_wgs = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')],
	y=muts_wgs[muts_wgs$DRIVER %in% c('yes','clonality'),wgs_selected_columns],
	by.x='DNA_SAMPLE', by.y='TARGET_NAME'
)
clonality_impact = merge(
	x=mdata[,c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE','AADx','SAMPLE_TYPE')],
	y=muts_impact[,impact_selected_columns],
	by.x='DNA_SAMPLE', by.y='Tumor_Sample_Barcode'
)
colnames(clonality_impact) = colnames(clonality_wgs)
clonality = rbind(clonality_impact, clonality_wgs)
clonality$VAG_VT[clonality$VAG_VT=='DEL'] = 'Del'
clonality$VAG_VT[clonality$VAG_VT=='SNP'] = 'Sub'
clonality$VAG_VT[clonality$VAG_VT=='DNP'] = 'Sub'
clonality$VAG_VT[clonality$VAG_VT=='INS'] = 'Ins'
clonality$VAG_VT[clonality$VAG_VT=='ONP'] = 'Sub'
clonality$VAG_EFFECT[clonality$VAG_EFFECT=="frameshit_variant"] = "frameshift_variant"
clonality$VAG_EFFECT[clonality$VAG_EFFECT=="Nonsense"] = "stop_gained"
clonality$VAG_EFFECT[clonality$VAG_EFFECT=="frameshift_variant,splice_region_variant"] = "frameshift_variant"
clonality$VAG_EFFECT[clonality$VAG_EFFECT %in% c("missense_variant",'missense','missense_variant,splice_region_variant')] = "non_synonymous_codon"
clonality$VAG_EFFECT[clonality$VAG_EFFECT %in% c('splice_acceptor_variant','splice_donor_variant','splice_acceptor_variant,intron_variant')] = "splice_site_variant"
clonality$VAG_GENE[clonality$VAG_GENE=='IGF2BP3-MDM2-CDK4'] = 'IGF2BP3'
clonality = unique(clonality)
clonality$cohort = 'MSK-WGS'
clonality$cohort[grepl('^P',clonality$PATIENT)] = 'MSK-IMPACT'

mdata$SAMPLE_TYPE = as.character(mdata$SAMPLE_TYPE)
mdata$SAMPLE_TYPE[mdata$SAMPLE_TYPE=='2nd_look'] = 't-resection'
mdata$SAMPLE_TYPE = factor(mdata$SAMPLE_TYPE, levels=c('diagnosis','reresection','t-resection','relapse','frelapse'))
mdata$DISEASE_SUBTYPE = factor(mdata$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
mdata$TUMOR_SITE3 = factor(mdata$TUMOR_SITE3, levels=c('primary_disease','local_disease','distant_metastasis'))
print('Cohort:')
print(paste(nrow(mdata), 'tumors from ', length(unique((mdata$PATIENT))), ' patients'), sep='')
print('Sample types:')
print(table(mdata$SAMPLE_TYPE[sel]))
print('Disease subtypes:')
print(table(mdata$DISEASE_SUBTYPE[sel]))
print('Tumor locations:')
print(table(mdata$TUMOR_SITE3[sel]))

