WDIR = '/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/util_dataset.R')
library(gridExtra)
library(ghibli)

######################################################### Extended Data Figure 1a_1
data_summary = mdata[,c('PATIENT','DNA_SAMPLE','SEQ','RNA_SAMPLE','SAMPLE_TYPE','DISEASE_SUBTYPE','TUMOR_SITE3')]
data_summary$RNA_SAMPLE[!is.na(data_summary$RNA_SAMPLE)] = 'yes'
data_summary$RNA_SAMPLE[is.na(data_summary$RNA_SAMPLE)] = 'no'
data_summary$SEQ[data_summary$SEQ=='IMPACT'] = 'MSK-IMPACT'
data_summary = merge(x=data_summary, y=as.data.frame(table(data_summary$PATIENT[data_summary$SEQ=='WGS'])), by.x='PATIENT', by.y='Var1', all.x=T)
colnames(data_summary)[ncol(data_summary)] = 'no_WGS_tumors'
data_summary$no_WGS_tumors[is.na(data_summary$no_WGS_tumors)] = 0
data_summary = merge(x=data_summary, y=as.data.frame(table(data_summary$PATIENT[data_summary$SEQ=='MSK-IMPACT'])), by.x='PATIENT', by.y='Var1', all.x=T)
colnames(data_summary)[ncol(data_summary)] = 'no_MSKIMPACT_tumors'
data_summary$no_MSKIMPACT_tumors[is.na(data_summary$no_MSKIMPACT_tumors)] = 0
data_summary = merge(x=data_summary, y=as.data.frame(table(data_summary$PATIENT)), by.x='PATIENT', by.y='Var1', all.x=T)
colnames(data_summary)[ncol(data_summary)] = 'no_tumors'
data_summary$SAMPLE_TYPE = as.character(data_summary$SAMPLE_TYPE)
data_summary$SAMPLE_TYPE[data_summary$SAMPLE_TYPE=='2nd_look'] = 'therapy_resection'
data_summary$SAMPLE_TYPE[data_summary$SAMPLE_TYPE=='frelapse'] = 'further_relapse'
data_summary$SAMPLE_TYPE = factor(data_summary$SAMPLE_TYPE, levels=c('diagnosis','reresection','therapy_resection','relapse','further_relapse'))
disease_site_cols = c(rgb(74/255, 112/255, 139/255), rgb(200/255, 113/255, 55/255), rgb(128/255, 0/255, 0/255))
data_summary$SEQ = factor(data_summary$SEQ, levels=c('WGS','MSK-IMPACT'))
seqtype_cols = c(rgb(199/255, 177/255, 156/255), rgb(253/255, 221/255, 160/255))
########### number of patients vs number of tumors
data_summary$clonality = 'single_tumor'
data_summary$clonality[data_summary$no_WGS_tumors>1] = "multi WGS"
data_summary$clonality[data_summary$no_WGS_tumors==1 & data_summary$no_MSKIMPACT_tumors>=1] = "WGS & MSK-IMPACT"
data_summary$clonality[data_summary$no_WGS_tumors==0 & data_summary$no_MSKIMPACT_tumors>=2] = "multi MSK-IMPACT"

plot_edf1a_1 = ggbarplot(
	as.data.frame(table(data_summary$SEQ, data_summary$RNA_SAMPLE)), x='Var1', y='Freq', fill='Var2',
	xlab="", ylab="No. tumors", palette=seqtype_cols, legend="right"
) + rotate_x_text(angle=45)
plot_edf1a_1 = ggpar(plot_edf1a_1, legend.title="RNA-seq available")

plot_edf1a_2 = ggbarplot(
	as.data.frame(table(data_summary$SEQ, data_summary$SAMPLE_TYPE)), x='Var1', y='Freq', fill='Var2',
	xlab="", ylab="No. tumors", palette=sample_type_cols$color, legend='right'
) + rotate_x_text(angle=45)
plot_edf1a_2 = ggpar(plot_edf1a_2, legend.title="Sample type")

plot_edf1a_3 = ggbarplot(
	as.data.frame(table(data_summary$SEQ, data_summary$DISEASE_SUBTYPE)), x='Var1', y='Freq', fill='Var2',
	xlab="", ylab="No. tumors", palette=disease_subtype_cols$color, legend='right'
) + rotate_x_text(angle=45)
plot_edf1a_3 = ggpar(plot_edf1a_3, legend.title="Disease subtype")

plot_edf1a_4 = ggbarplot(
	as.data.frame(table(data_summary$SEQ, data_summary$TUMOR_SITE3)), x='Var1', y='Freq', fill='Var2',
	xlab="Tumor site", ylab="No. tumors", palette=disease_site_cols, legend='right'
) + rotate_x_text(angle=45)
plot_edf1a_4 = ggpar(plot_edf1a_4, legend.title="Disease site")

pdf('raw/ext_data_fig_1a_1.pdf', height=7, width=7)
plot_grid(plot_edf1a_1, plot_edf1a_2, plot_edf1a_3, plot_edf1a_4, nrow=2)
dev.off()

pdf('raw/ext_data_fig_1a_2.pdf', height=3.5, width=3.5)
ggbarplot(as.data.frame(table(unique(data_summary[data_summary$clonality!="single_tumor",c('PATIENT','clonality')])$clonality)), x='Var1', y='Freq', fill='grey', ylab="No. patients", xlab="Type of sequencing") + rotate_x_text(angle=45)
dev.off()

######################################################### Extended Data Figure 1a_2
evol = read.delim('data/evolutionary_patterns.txt', header=T, sep="\t")
plot_grid(plot_edf1a_1, plot_edf1a_2, plot_edf1a_3, plot_edf1a_4, plot_edf1a_5, nrow=2)

######################################################### Extended Data Figure 1b
# Sample-level mutation counts
# Therapy status
mdata_wgs = mdata[mdata$SS_ANALYSIS=='yes',]
mdata_wgs$DNA_SAMPLE = factor(mdata_wgs$DNA_SAMPLE, levels=mdata_wgs$DNA_SAMPLE[order(mdata_wgs$NU_SUBS)])
mdata_wgs$THERAPY_VALUE = 1
p_therapy <- ggbarplot(mdata_wgs, x='DNA_SAMPLE', y='THERAPY_VALUE', color='THERAPY', fill='THERAPY', palette=c('grey','brown'))
p_therapy <- p_therapy + rremove('x.axis') + rremove('y.axis') + rremove('xylab') + rremove('xy.title') + rremove('xy.text') + rremove('x.ticks') + rremove('y.ticks') + rremove('legend')
p_therapy <- facet(p_therapy, facet.by="DISEASE_SUBTYPE", scales='free_x', nrow=1)
gp_therapy <- ggplotGrob(p_therapy)
facet.columns <- gp_therapy$layout$l[grepl("panel", gp_therapy$layout$name)]
x.var <- sapply(ggplot_build(p_therapy)$layout$panel_scales_x, function(l) length(l$range$range))
gp_therapy$widths[facet.columns] <- gp_therapy$widths[facet.columns] * x.var
p_subs <- ggbarplot(mdata_wgs, x='DNA_SAMPLE', y='NU_SUBS', color='grey45', fill='SAMPLE_TYPE', palette=sample_type_cols$color, ylab="No. substitutions")
p_subs <- facet(p_subs, facet.by = "DISEASE_SUBTYPE", scale="free_x", nrow=1)
p_subs <- p_subs + rremove('x.axis') + rremove('xlab') + rremove('x.title') + rremove('x.text') + rremove('x.ticks') + rremove('legend')
gp_subs <- ggplotGrob(p_subs)
facet.columns <- gp_subs$layout$l[grepl("panel", gp_subs$layout$name)]
x.var <- sapply(ggplot_build(p_subs)$layout$panel_scales_x, function(l) length(l$range$range))
gp_subs$widths[facet.columns] <- gp_subs$widths[facet.columns] * x.var
p_indel <- ggbarplot(mdata_wgs, x='DNA_SAMPLE', y='NU_INDELS', color='grey45', fill='SAMPLE_TYPE', palette=sample_type_cols$color, ylab="No. indels")
p_indel <- facet(p_indel, facet.by = "DISEASE_SUBTYPE", scale="free_x", nrow=1)
p_indel <- p_indel + rremove('x.axis') + rremove('xlab') + rremove('x.title') + rremove('x.text') + rremove('x.ticks') + rremove('legend')
gp_indel <- ggplotGrob(p_indel)
facet.columns <- gp_indel$layout$l[grepl("panel", gp_indel$layout$name)]
x.var <- sapply(ggplot_build(p_indel)$layout$panel_scales_x, function(l) length(l$range$range))
gp_indel$widths[facet.columns] <- gp_indel$widths[facet.columns] * x.var
p_sv <- ggbarplot(mdata_wgs, x='DNA_SAMPLE', y='NU_SV', color='grey45', fill='SAMPLE_TYPE', palette=sample_type_cols$color, ylab="No. SVs")
p_sv <- facet(p_sv, facet.by = "DISEASE_SUBTYPE", scale="free_x", nrow=1)
p_sv <- p_sv + rremove('x.axis') + rremove('xlab') + rremove('x.title') + rremove('x.text') + rremove('x.ticks') + rremove('legend')
gp_sv <- ggplotGrob(p_sv)
facet.columns <- gp_sv$layout$l[grepl("panel", gp_sv$layout$name)]
x.var <- sapply(ggplot_build(p_sv)$layout$panel_scales_x, function(l) length(l$range$range))
gp_sv$widths[facet.columns] <- gp_sv$widths[facet.columns] * x.var
pdf('raw/ext_data_fig_1b.pdf', height=10, width=20)
plot_grid(gp_therapy, gp_subs, gp_indel, gp_sv, align='v', ncol=1, rel_heights=c(1,4,4,4))
dev.off()
######################################################### Extended Data Figure 1c
mut_counts = as.data.frame(sort(table(unique(oncogenic[,c('PATIENT','VAG_GENE')])$VAG_GENE), decreasing=T))
colnames(mut_counts) = c('Gene','N_patients')
mut_counts$Prop_patients = mut_counts$N_patients / length(unique(mdata$PATIENT))
mut_counts = mut_counts[mut_counts$N_patients>1,]
selected_genes = as.character(mut_counts$Gene)
gene_mut_type = unique(oncogenic[,c('PATIENT','VAG_GENE','VAG_VT')])
gene_mut_type = as.data.frame(table(gene_mut_type$VAG_GENE, gene_mut_type$VAG_VT))
colnames(gene_mut_type) = c('Gene','Variant_type','N_patients')
gene_mut_type = gene_mut_type[gene_mut_type$Gene %in% selected_genes,]
gene_mut_type$Gene = factor(gene_mut_type$Gene, levels=selected_genes)
for(gene in unique(gene_mut_type$Gene)) {
	gene_mut_type[gene_mut_type$Gene==gene,'Prop_mutations'] = gene_mut_type[gene_mut_type$Gene==gene,'N_patients'] / sum(gene_mut_type[gene_mut_type$Gene==gene,'N_patients'])
}
plot_genes = ggbarplot(mut_counts, x='Gene', y="Prop_patients", fill='darkgrey') + rotate_x_text(90) + geom_hline(yintercept=0.1, linetype="dashed", color='red')
plot_mut_type = ggbarplot(gene_mut_type, x='Gene', y="Prop_mutations", fill="Variant_type", palette=ghibli_palettes$LaputaMedium[2:7])
plot_mut_type = plot_mut_type + scale_y_reverse() + rremove("x.axis") + rremove("x.text") + rremove('xlab') + rremove('x.ticks')
pdf('raw/ext_data_fig_1c.pdf', width=20)
plot_grid(plot_genes, plot_mut_type, align='v', nrow=2)
dev.off()
##################### Extended Data Figure 1d
library(survival)
library(survminer)
patients = unique(mdata[,c('PATIENT','DISEASE_SUBTYPE','AADx','os_status','os_time_months','AGE_GROUP','STAGE')])
patients_sub = patients[!is.na(patients$os_time_months),]
oncogenic_uniq = unique(oncogenic[,c('PATIENT','VAG_GENE')])
muts_table = table(oncogenic_uniq$PATIENT, oncogenic_uniq$VAG_GENE)
muts_table = as.data.frame(muts_table[,which(colSums(muts_table)>5)])
muts_table$Freq[muts_table$Freq>0] = 1
muts_table = dcast(muts_table, Var1 ~ Var2, value.var="Freq")
colnames(muts_table)[1] = 'patient'
colnames(muts_table)[grep('gain',colnames(muts_table))] = paste('cna',grep('gain',colnames(muts_table),value=T),sep="_")
colnames(muts_table)[grep('loss',colnames(muts_table))] = paste('cna',grep('loss',colnames(muts_table),value=T),sep="_")
colnames(muts_table)[grep('cnloh',colnames(muts_table))] = paste('cna',grep('cnloh',colnames(muts_table),value=T),sep="_")
colnames(muts_table)[!grepl('cna',colnames(muts_table))] = paste('mut_',colnames(muts_table)[!grepl('cna',colnames(muts_table))], sep="")
patients_sub = merge(x=patients_sub, y=muts_table, by.x='PATIENT', by.y='mut_patient', all.x=T)
################### Check TERT promoter mutations
patients_mycn = patients_sub[patients_sub$DISEASE_SUBTYPE=='MYCN',]
patients_mycn$mut_TERTp = 0
patients_mycn$mut_TERTp[patients_mycn$PATIENT %in% unique(oncogenic[oncogenic$VAG_GENE=='TERT' & oncogenic$VAG_VT=='Sub' & oncogenic$DISEASE_SUBTYPE=='MYCN','PATIENT'])] = 1
pval_gene = summary(coxph(Surv(os_time_months, os_status) ~ AGE_GROUP + mut_TERTp, data=patients_mycn))$coefficients['mut_TERTp','Pr(>|z|)']
patients_mycn$mut_TERTp[patients_mycn$PATIENT %in% unique(oncogenic[oncogenic$VAG_GENE=='TERT' & oncogenic$VAG_VT!='Sub' & oncogenic$DISEASE_SUBTYPE=='MYCN','PATIENT'])] = 2
pt_surv_fit = survfit(Surv(time=patients_mycn$os_time_months, event=patients_mycn$os_status) ~ patients_mycn$mut_TERTp, data=patients_mycn)
splot_tertp = ggsurvplot(
	fit=pt_surv_fit, data=patients_mycn, xlab='Months', ylab='Survival probability', conf.int = TRUE, conf.int.fill="purple", surv.median.line="hv",, title=paste('MYCN-A patients N=',nrow(patients_mycn),sep=''),
	risk.table=TRUE, risk.table.height=0.25, risk.table.fontsize=4.5, pval=round(pval_gene, digits=4), pval.size=5, pval.coord=c(0,0), palette="Dark2", legend.labs=c('TERT WT','TERTp','TERTsv')
)
pdf('raw/ext_data_fig_1d.pdf')
splot_tertp
dev.off()
##################### Extended Data Figure 1e Comutation patterns
source('scripts/util_plot_comutation.R')
oncogenic_uniq = unique(oncogenic[,c('PATIENT','VAG_GENE','POS_KEY')])
comutation_data = table(oncogenic_uniq$PATIENT, oncogenic_uniq$VAG_GENE)
comutation_data[comutation_data>1] = 1
plot_comutation(comutation_data, frequent=T, cutoff.frequent=2, cache=F, out.dir='raw')
system('mv raw/comutation_pattern.pdf raw/ext_data_fig_1e.pdf')
system(paste('rm -r raw/cache', sep=''))
