WDIR = '~/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/util_dataset.R')
# Extended Data Figure 3a
plot_sbs18 <- ggscatter(
	data=mdata[mdata$SEQ=='WGS' & mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','t-resection'),], x='A_SBS18_n', y='AADx', add="reg.line", conf.int=TRUE,
	add.params=list(color="blue", fill="lightgray"), xlab='#Mutations A_SBS18', ylab='Age at diagnosis (m)', title='All diagnostic/t. resection tumors'
)
plot_sbs18 <- plot_sbs18 + stat_cor(method = "pearson", label.x=2000, label.y=200)
plot_sbs38 <- ggscatter(
	data=mdata[mdata$SEQ=='WGS' & mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','t-resection'),], x='B_SBS38_n', y='AADx', add="reg.line", conf.int=TRUE,
	add.params=list(color="blue", fill="lightgray"), xlab='#Mutations B_SBS38', ylab='Age at diagnosis (m)', title='All diagnostic/t. resection tumors'
) 
plot_sbs38 <- plot_sbs38 + stat_cor(method = "pearson", label.x=500, label.y=100)
plot_sbs40 <- ggscatter(
	data=mdata[mdata$SEQ=='WGS' & mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','t-resection'),], x='C_SBS40_n', y='AADx', add="reg.line", conf.int=TRUE,
	add.params=list(color="blue", fill="lightgray"), xlab='#Mutations C_SBS40', ylab='Age at diagnosis (m)', title='All diagnostic/t. resection tumors'
) 
plot_sbs40 <- plot_sbs40 + stat_cor(method = "pearson", label.x=500, label.y=100)
pdf('raw/ext_data_fig_3a.pdf', width=14)
plot_grid(plot_sbs18, plot_sbs40, nrow=1)
dev.off()
################ Check correlation of A_SBS18 with glutaminolysis signature
glutaminolysis_genes = c(
	'ME2','SLC7A11','SLC7A2','SLC6A6','GLS2','SLC6A15','ME3','SLC16A10','SLC7A6','SLC7A10','SLC7A5',
	'SLC3A2','SLC38A5','SLC1A5','SLC43A1','GPT','PFAS','CAD','SLC7A1','PPAT','ALDH18A1','GPT2','GOT2'
)
# Read expression data
mdata$RNA_SAMPLE2 = gsub('-','.',mdata$RNA_SAMPLE)
exp_rpkm = read.table('data/all_expression_coding_rpkm.txt.gz', row.names=1, header=T, sep='\t')
selected_genes = which(rowSums(exp_rpkm) > summary(rowSums(exp_rpkm))['Median'])
exp_filtered = exp_rpkm[selected_genes,]
exp_filtered = log2(exp_filtered+1)
exp_filtered = exp_filtered - apply(exp_filtered, 1, median)

mean_exp = as.data.frame(colSums(exp_rpkm[glutaminolysis_genes,]) / length(glutaminolysis_genes))
colnames(mean_exp) = 'glutaminolysis_signature'
mean_exp$sample = rownames(mean_exp)
mean_exp = merge(y=mdata[mdata$PURITY>0.2,c('RNA_SAMPLE2','SAMPLE_TYPE','DISEASE_SUBTYPE','A_SBS18','A_SBS18_n')], x=mean_exp, by.y='RNA_SAMPLE2', by.x='sample')

compare_dsubtype <- list( c("MYCN", "TERT"), c("MYCN","ATRX"), c('MYCN','SEG_CNV'), c("TERT","ATRX"), c("TERT", "SEG_CNV"), c("ATRX","SEG_CNV"))
symnum_args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
plot_exp_dx <- ggboxplot(mean_exp[mean_exp$SAMPLE_TYPE %in% c('diagnosis','reresection'),], x='DISEASE_SUBTYPE', y='glutaminolysis_signature', fill='DISEASE_SUBTYPE', palette=disease_subtype_cols$color, add='jitter')
plot_exp_dx <- plot_exp_dx + stat_compare_means(comparisons = compare_dsubtype, symnum.args=symnum_args) + rremove("legend")

compare_stype <- list( c("diagnosis", "2nd_look"), c("diagnosis","relapse"), c('diagnosis','frelapse'), c("relapse","frelapse"), c("TERT", "SEG_CNV"), c("ATRX","SEG_CNV"))
plot_exp_mycn <- ggboxplot(mean_exp[mean_exp$DISEASE_SUBTYPE=='MYCN',], x='SAMPLE_TYPE', y='glutaminolysis_signature', fill='SAMPLE_TYPE', palette=sample_type_cols$color, add='jitter')
plot_exp_mycn <- plot_exp_mycn + stat_compare_means(comparisons = compare_stype, symnum.args=symnum_args) + rremove("legend")

pdf('raw/ext_data_fig_3b.pdf', width=14)
plot_grid(plot_exp_dx, plot_exp_mycn, nrow=1)
dev.off()
# Extended Data Figure 3c
xcell = read.table('data/xCell_predictions.txt', header=T, sep=",")
colnames(xcell)[1] = 'sample'
xcell$sample = gsub('\\.','-',xcell$sample)
xcell = merge(x=xcell, y=mdata[,c('DNA_SAMPLE','RNA_SAMPLE')], by.x='sample', by.y='RNA_SAMPLE')
neoantigen = read.table('data/neoantigen_predictions.csv', header=T, sep=",")
colnames(neoantigen)[1] = 'sample'
xcell = merge(x=xcell, y=neoantigen, by.x='DNA_SAMPLE', by.y='sample')

plot_neo = ggscatter(neoantigen, x='Mutations', y='Neoantigens', add="reg.line", conf.int=T, add.params=list(color='blue', fill='lightgray'), ylab="Number of predicted neoantigens", xlab='Number of SNVs') + stat_cor(method = "pearson", label.x=1000, label.y=60)
plot_xcell = ggscatter(xcell, x='MicroenvironmentScore', y='Neoantigens', add="reg.line", conf.int=T, add.params=list(color='blue', fill='lightgray'), ylab='Number of predicted neoantigens', xlab='xCell ImmuneScore') + stat_cor(method = "pearson", label.x=0.3, label.y=60)

pdf('raw/ext_data_fig_3c.pdf', width=14)
plot_grid(plot_neo, plot_xcell, nrow=1)
dev.off()

# Extended Data Figure 3d

source('scripts/util_compare_clones.R')
mdata$DNA_SAMPLE2  = gsub('-','.',mdata$DNA_SAMPLE)
samples_rx_naive = mdata[which(mdata$SEQ=='WGS' & mdata$THERAPY=='no'),c('DNA_SAMPLE2')]
samples_post_rx = mdata[which(mdata$SEQ=='WGS' & mdata$THERAPY=='yes'),c('DNA_SAMPLE2')]
subclones_post_rx = unique(subs_clones[subs_clones$sample %in% samples_post_rx & subs_clones$descr!="trunk",c('patient','descr','no.muts.in.cluster')])
subclones_post_rx$clone_id = apply(subclones_post_rx[,c('patient','descr')], 1, paste, collapse="__")

clone_signatures = read.table('data/patient_level_signature_exposures.txt', header=T, sep="\t", stringsAsFactors=F)

clone_signatures$relapse_specific = 'no'
clone_signatures$relapse_specific[clone_signatures$clone %in% subclones_post_rx$clone_id] = 'yes'
total_clone_signature_counts = melt(ddply(clone_signatures, c('relapse_specific'), summarise, A_SBS18=sum(A_SBS18), B_SBS38=sum(B_SBS38), C_SBS40=sum(C_SBS40), D_TMZ=sum(D_TMZ), E_SBS31=sum(E_SBS31), F_SBS35=sum(F_SBS35)), id.vars=c('relapse_specific'))
total_clone_signature_counts$proportion = NA
values = total_clone_signature_counts[total_clone_signature_counts$relapse_specific=='yes','value']
total_clone_signature_counts[total_clone_signature_counts$relapse_specific=='yes','proportion'] = values / sum(values)
values = total_clone_signature_counts[total_clone_signature_counts$relapse_specific=='no','value']
total_clone_signature_counts[total_clone_signature_counts$relapse_specific=='no','proportion'] = values / sum(values)
plot_all = ggbarplot(total_clone_signature_counts, x='relapse_specific', y='proportion', fill='variable', palette=c('#1A83B8','#015580','#899DA4','#178520','#AD381F','#801701'), xlab="is_relapse_specific", legend='left', title="All SNVs")

drivers = read.table('data/driver_substitutions_for_signature_analysis_output.txt', header=T, sep="\t", stringsAsFactors=F)
drivers$class[drivers$class!="relapse_only"] = 'no'
drivers$class[drivers$class=="relapse_only"] = 'yes'
total_driver_signature_counts = as.data.frame(table(drivers$assigned_sig, drivers$class))
total_driver_signature_counts$proportion = NA
values = total_driver_signature_counts[total_driver_signature_counts$Var2=='no','Freq']
total_driver_signature_counts[total_driver_signature_counts$Var2=='no','proportion'] = values / sum(values)
values = total_driver_signature_counts[total_driver_signature_counts$Var2=='yes','Freq']
total_driver_signature_counts[total_driver_signature_counts$Var2=='yes','proportion'] = values / sum(values)
plot_drivers = ggbarplot(total_driver_signature_counts, x='Var2', y='proportion', fill='Var1', palette=c('#1A83B8','#015580','#899DA4','#AD381F','#801701'), xlab="is_relapse_specific", ylab="No. SNVs", legend='none', title='Driver SNVs')

pdf('raw/ext_data_fig_3d.pdf')
plot_grid(plot_all, plot_drivers, rel_widths=c(1,0.65))
dev.off()
