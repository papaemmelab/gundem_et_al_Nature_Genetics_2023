WDIR = '/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/util_dataset.R')

library(circlize)
library('Ternary')
library(wesanderson)
library(cowplot)
library(ComplexHeatmap)
comparison_colors = c('grey', wes_palette('Darjeeling1', 5, type="discrete")[5])
##################### MSK patients and samples broken down by sequencing
mdata$SAMPLE_TYPE[mdata$SAMPLE_TYPE=='autopsy'] = 'frelapse'
mdata = mdata[mdata$SS_ANALYSIS!='no' | mdata$MS_ANALYSIS=='yes',]
patients_impact = unique(mdata[mdata$COHORT=='IMPACT',c('PATIENT','AADx','DISEASE_SUBTYPE','os_status','os_time_months','no_tumor_samples','COHORT')])
patients_wgs = unique(mdata[mdata$COHORT=='WGS',c('PATIENT','AADx','DISEASE_SUBTYPE','os_status','os_time_months','no_tumor_samples','COHORT')])
patients_all = merge(x=rbind(patients_impact, patients_wgs), y=as.data.frame(table(unique(mdata[,c('PATIENT','SAMPLING_DATE')])$PATIENT)), by.x='PATIENT', by.y='Var1')
colnames(patients_all) = c('PATIENT','AADx','DISEASE_SUBTYPE','OS_STATUS','OS_TIME_MONTHS','NO_SAMPLES_PT','ASSAY','NO_TIME_POINTS_PT')
samples_wgs = mdata[mdata$COHORT=='WGS',c('PATIENT','COHORT','SAMPLING_DATE','DNA_SAMPLE','SAMPLE_TYPE','DODx')]
samples_impact = mdata[mdata$COHORT=='IMPACT',c('PATIENT','COHORT','SAMPLING_DATE','DNA_SAMPLE','SAMPLE_TYPE','DODx')]
##################### Ternary plot with genes
samples = rbind(samples_impact, samples_wgs)
patients_with_multiple_tumors = names(which(table(samples$PATIENT)>1))
oncogenic_sub = unique(oncogenic[oncogenic$PATIENT %in% patients_with_multiple_tumors & !is.na(oncogenic$subclone_type),c('PATIENT','VAG_GENE','POS_KEY','subclone_type','VAG_EFFECT')])
counts = table(oncogenic_sub$VAG_GENE, oncogenic_sub$subclone_type)
shared = as.data.frame(counts[,'trunk'])
colnames(shared)[1] = 'N'
shared$gene = rownames(shared)
rownames(shared) = NULL
primary = as.data.frame(rowSums(counts[,grepl('primary',colnames(counts)) & !grepl('met',colnames(counts))]))
colnames(primary)[1] = 'N'
primary$gene = rownames(primary)
rownames(primary) = NULL
metastasis = as.data.frame(rowSums(counts[,grepl('met',colnames(counts))]))
colnames(metastasis)[1] = 'N'
metastasis$gene = rownames(metastasis)
rownames(metastasis) = NULL
counts = merge(x=merge(x=shared, y=primary, by.x='gene', by.y='gene', all=T), y=metastasis, by.x='gene', by.y='gene', all=T)
colnames(counts) = c('gene','shared','primary','metastasis')
counts = counts[rowSums(data.matrix(counts[,2:4]))>=5,]
counts_matrix = t(data.matrix(counts[,c('shared','metastasis','primary')]))
data_points = split(counts_matrix, rep(1:ncol(counts_matrix), each = nrow(counts_matrix)))
names(data_points) = counts$gene
pdf('raw/figure_3a_genes.pdf')
point_cols = vapply(data_points, function (x) rgb(x[1]/sum(x)*255, x[2]/sum(x)*255, x[3]/sum(x)*255, maxColorValue=255), character(1))
TernaryPlot(atip='Shared', btip='Metastasis/Relapse', ctip='Primary', tip.col=c('red','green','blue'), grid.lines=5, grid.lty='dotted',grid.minor.lines=1, grid.minor.lty='dotted')
AddToTernary(points, data_points, pch=21, cex=colSums(counts_matrix)*0.1, bg=point_cols)
AddToTernary(text, data_points, names(data_points), cex=0.7, font=2)
legend('topright', title='Count', legend=c(10,20,30), cex=0.7, bty='n', pch=21, pt.cex=c(10*0.2,20*0.2,30*0.2))
legend('topleft', 
	legend = c("Shared", "Primary", "Metastasis/relapse"),
	cex = 0.8, bty = "n", pch = 21, pt.cex = 1.8,
	pt.bg = c(rgb(255/255, 0, 0), 
		rgb(0, 0, 255/255),
		rgb(0, 255/255, 0))
)
dev.off()
##################### Ternary plot with pathways
samples = rbind(samples_impact, samples_wgs)
patients_with_multiple_tumors = names(which(table(samples$PATIENT)>1))
oncogenic_sub = unique(oncogenic[oncogenic$PATIENT %in% patients_with_multiple_tumors & !is.na(oncogenic$subclone_type),c('PATIENT','VAG_GENE','POS_KEY','subclone_type')])
counts = table(oncogenic_sub$VAG_GENE, oncogenic_sub$subclone_type)
shared = as.data.frame(counts[,'trunk'])
colnames(shared)[1] = 'N'
shared$gene = rownames(shared)
rownames(shared) = NULL
primary = as.data.frame(rowSums(counts[,grepl('primary',colnames(counts)) & !grepl('met',colnames(counts))]))
colnames(primary)[1] = 'N'
primary$gene = rownames(primary)
rownames(primary) = NULL
metastasis = as.data.frame(rowSums(counts[,grepl('met',colnames(counts))]))
colnames(metastasis)[1] = 'N'
metastasis$gene = rownames(metastasis)
rownames(metastasis) = NULL
counts = merge(x=merge(x=shared, y=primary, by.x='gene', by.y='gene', all=T), y=metastasis, by.x='gene', by.y='gene', all=T)
colnames(counts) = c('gene','shared','primary','metastasis')
#counts = counts[!grepl('_loss', counts$gene) & !grepl('_gain', counts$gene) & !grepl('_cnloh', counts$gene),]
genes_MAPK_pathway = c('ALK','BRAF','CIC','EGFR','ERBB4','FGFR1','FLT3','HRAS','KRAS','NF1','NRAS','PTPN11','SOS1')
genes_differentiation = c('RARA','RARB','PHOX2B','WNT5A','IGF2BP3','SPRY2')
genes_PI3K_mTOR_pathway = c('PIK3CA','FBXW7','MTOR','PIK3R1','TSC1','TSC2')
genes_cell_cycle = c('CCND1','CDK4','CDKN1B','CDKN1C','CDKN2A')
genes_chromatin_remodeling = c('ARID1A','ARID1B','CREBBP','EP300','SMARCA4','BCOR')
genes_neurodevelopment = c('PTPRD','DLG2','AUTS2','SHANK2','TIAM1','DLC1')
genes_TP53_pathway = c('TP53','MDM2','MDM4')
counts_pathways = as.data.frame(rbind(
        c('MAPK_pathway',colSums(counts[counts$gene %in% genes_MAPK_pathway,2:4])),
        c('PI3K_mTOR_pathway',colSums(counts[counts$gene %in% genes_PI3K_mTOR_pathway,2:4])),
        c('Chromatin_remodeling',colSums(counts[counts$gene %in% genes_chromatin_remodeling,2:4])),
        c('Cell_cycle',colSums(counts[counts$gene %in% genes_cell_cycle,2:4])),
        c('Differentiation',colSums(counts[counts$gene %in% genes_differentiation,2:4])),
        c('Neurodevelopment',colSums(counts[counts$gene %in% genes_neurodevelopment,2:4])),
        c('TP53_pathway',colSums(counts[counts$gene %in% genes_TP53_pathway,2:4]))
), stringsAsFactors=F)
colnames(counts_pathways)[1] = 'gene'
counts = counts[!(counts$gene %in% c(genes_TP53_pathway, genes_MAPK_pathway, genes_PI3K_mTOR_pathway, genes_cell_cycle, genes_chromatin_remodeling, genes_differentiation, genes_neurodevelopment)),]
counts = counts[rowSums(data.matrix(counts[,2:4]))>=5,]
counts = rbind(counts_pathways, counts)
counts$shared = as.integer(counts$shared)
counts$primary = as.integer(counts$primary)
counts$metastasis = as.integer(counts$metastasis)
counts_matrix = t(data.matrix(counts[,c('shared','metastasis','primary')]))
data_points = split(counts_matrix, rep(1:ncol(counts_matrix), each = nrow(counts_matrix)))
names(data_points) = counts$gene
pdf('raw/figure_3a_pathways.pdf')
point_cols = vapply(data_points, function (x) rgb(x[1]/sum(x)*255, x[2]/sum(x)*255, x[3]/sum(x)*255, maxColorValue=255), character(1))
par('family'='sans')
TernaryPlot(atip='Shared', btip='Metastasis/Relapse', ctip='Primary', tip.col=c('red','green','blue'), grid.lines=5, grid.lty='dotted',grid.minor.lines=1, grid.minor.lty='dotted')
AddToTernary(points, data_points, pch=21, cex=colSums(counts_matrix)*0.1, bg=point_cols)
AddToTernary(text, data_points, names(data_points), cex=0.7, font=2)
legend('topright', title='Count', legend=c(10,20,30), cex=0.7, bty='n', pch=21, pt.cex=c(10*0.2,20*0.2,30*0.2))
legend('topleft',
	legend = c("Shared", "Primary", "Metastasis/relapse"),
	cex = 0.8, bty = "n", pch = 21, pt.cex = 1.8,
	pt.bg = c(rgb(255/255, 0, 0),
		rgb(0, 0, 255/255),
		rgb(0, 255/255, 0))
)

dev.off()
##################### Prepare data for heatmap
oncogenic_sub = unique(oncogenic[oncogenic$PATIENT %in% patients_with_multiple_tumors,c('PATIENT','VAG_GENE','POS_KEY','subclone_type')])
oncogenic_sub$trunk = 1
oncogenic_sub$trunk[oncogenic_sub$subclone_type!="trunk"] = 0
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_MAPK_pathway] = 'MAPK_pathway'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_PI3K_mTOR_pathway] = 'PI3K_mTOR_pathway'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_chromatin_remodeling] = 'Chromatin_remodeling'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_neurodevelopment] = 'Neurodevelopment'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_differentiation] = 'Differentiation'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_cell_cycle] = 'Cell_cycle'
oncogenic_sub$VAG_GENE[oncogenic_sub$VAG_GENE %in% genes_TP53_pathway] = 'TP53_pathway'
oncogenic_sub = oncogenic_sub[!(oncogenic_sub$VAG_GENE %in% c('BCOR','ARID4B','BRCA2')),]
oncogenic_sub$gene = paste(oncogenic_sub$VAG_GENE, oncogenic_sub$trunk, sep=" ")
# Plotting order of for gene
dplot = table(oncogenic_sub$PATIENT, oncogenic_sub$gene)
patients_wo_mutation = setdiff(patients_with_multiple_tumors, rownames(dplot))
patients_wo_mutation_data = mat.or.vec(length(patients_wo_mutation), ncol(dplot))
rownames(patients_wo_mutation_data) = patients_wo_mutation
dplot = rbind(dplot, patients_wo_mutation_data)
freq_genes = names(which(colSums(dplot)>1))
#freq_genes = colnames(dplot)
offtrunk_genes_order = intersect(names(sort(table(oncogenic_sub$gene[grep('0$',oncogenic_sub$gene)]), decreasing=T)), freq_genes)
offtrunk_cna_order = offtrunk_genes_order[grepl('gain',offtrunk_genes_order) | grepl('loss',offtrunk_genes_order) | grepl('cnloh',offtrunk_genes_order)]
offtrunk_genes_order = setdiff(offtrunk_genes_order, offtrunk_cna_order)
offtrunk_cna_order = setdiff(offtrunk_cna_order, c('6q_loss 0','9p_loss 0','1p_loss 0','17p_cnloh 0'))
trunk_genes_order = intersect(names(sort(table(oncogenic_sub$gene[grep('1$',oncogenic_sub$gene)]), decreasing=T)), freq_genes)
trunk_cna_order = trunk_genes_order[grepl('gain',trunk_genes_order) | grepl('loss',trunk_genes_order) | grepl('cnloh',trunk_genes_order)]
trunk_genes_order = setdiff(trunk_genes_order, trunk_cna_order)
trunk_cna_order = setdiff(trunk_cna_order,c('17p_cnloh 1','4p_loss 1','9p_loss 1','6q_loss 1'))
gorder = c(trunk_cna_order, trunk_genes_order, offtrunk_cna_order, offtrunk_genes_order)
# Plotting order for patients
patients_mycn = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='MYCN'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='MYCN']), rownames(dplot))
patients_tert = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='TERT'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='TERT']), rownames(dplot))
patients_mdm2_ckd4 = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='MDM2-CDK4'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='MDM2-CDK4']), rownames(dplot))
patients_seg_cnv = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='SEG_CNV'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='SEG_CNV']), rownames(dplot))
patients_num_cnv = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='NUMERIC_CNV'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='NUMERIC_CNV']), rownames(dplot))
patients_atrx = intersect(c(patients_wgs$PATIENT[patients_wgs$DISEASE_SUBTYPE=='ATRX'], patients_impact$PATIENT[patients_impact$DISEASE_SUBTYPE=='ATRX']), rownames(dplot))
porder = c(patients_mycn, patients_tert, patients_atrx, patients_mdm2_ckd4, patients_seg_cnv, patients_num_cnv)
porder = intersect(porder, rownames(dplot))
# Re-order data matrix
dplot = dplot[porder,gorder]
#Patient info for heatmap row annotation
patients_all = patients_all[patients_all$PATIENT %in% rownames(dplot),]
patients_all = patients_all[match(porder,patients_all$PATIENT),]
patients_with_multiwgs = names(which(table(mdata$PATIENT[mdata$SEQ=='WGS'])>1))
patients_all$multi_WGS = 'no'
patients_all$multi_WGS[patients_all$PATIENT %in% patients_with_multiwgs] = 'yes'
col_fun_nsamples = colorRamp2(c(1, max(patients_all$NO_SAMPLES_PT)), c("white", "darkgreen"))
col_fun_ntimes = colorRamp2(c(1, max(patients_all$NO_TIME_POINTS_PT)), c("white", "purple"))
bottom_ha_patients = HeatmapAnnotation(
	multi_WGS=patients_all$multi_WGS,
	no_samples=anno_text(
		patients_all$NO_SAMPLES_PT, location=0.5,just='center', rot=0,
		gp=gpar(col="black", fill=col_fun_nsamples(patients_all$NO_SAMPLES_PT), fontfamily='sans', fontsize=10)#, width=max_text_width(patients_all$NO_SAMPLES_PT)*2)
        ),
	no_time_points=anno_text(
	patients_all$NO_TIME_POINTS_PT, location=0.5,just='center', rot=0,
	gp=gpar(col="black", fill=col_fun_nsamples(patients_all$NO_TIME_POINTS_PT), fontfamily='sans', fontsize=10)#, width=max_text_width(patients_all$NO_SAMPLES_PT)*2)
	),
	col = list(multi_WGS=c('no'='white', 'yes'='grey'))
)
# Fix the patients with parallel evolution at the same gene/pathway
dplot[dplot>3] = 3
dplot_table = t(dplot)
dplot_table[dplot_table==0] = ""
dplot_table[dplot_table=="1"] = "1 event"
dplot_table[dplot_table=="2"] = "2 events"
dplot_table[dplot_table=="3"] = ">=3 events"
sel_patients = c('H103207','H116991','H133121','H136083','P-0036611','H135467')
dplot_table['MAPK_pathway 0',sel_patients] = paste(dplot_table['MAPK_pathway 0',sel_patients], '|parallel', sep='')
sel_patients = c('H118733','H134817')
dplot_table['Neurodevelopment 0',sel_patients] = paste(dplot_table['Neurodevelopment 0',sel_patients], '|parallel', sep='')
sel_patients = c('P-0003149')
dplot_table['Chromatin_remodeling 0',sel_patients] = paste(dplot_table['Chromatin_remodeling 0',sel_patients], '|parallel', sep='')
sel_patients = c('H116987','H134722')
dplot_table['TP53_pathway 0',sel_patients] = paste(dplot_table['TP53_pathway 0',sel_patients], '|parallel', sep='')
sel_patients = c('H134817','H135089')
dplot_table['ATRX 0',sel_patients] = paste(dplot_table['ATRX 0',sel_patients], '|parallel', sep='')
sel_patients = c('H134819')
dplot_table['Differentiation 0',sel_patients] = paste(dplot_table['Differentiation 0',sel_patients], '|parallel', sep='')
sel_patients = c('H134821','P-0006215')
dplot_table['PI3K_mTOR_pathway 0',sel_patients] = paste(dplot_table['PI3K_mTOR_pathway 0',sel_patients], '|parallel', sep='')
col = c('1 event'="#D69C4E", '2 events'="#74A089", '>=3 events'="#046C9A")
# Fix the patients with parallel/continuous SV evolution at the same gene
sel_patients = c('H103207','H116984','H116986','H116987','H116989','H116991','H118706','H132374','H132381')
dplot_table['MYCN 0',sel_patients] = paste(dplot_table['MYCN 0',sel_patients], '|continuous', sep='')
sel_patients = c('H112908','H132379','H132380','H132386','H132384')
dplot_table['TERT 0',sel_patients] = paste(dplot_table['TERT 0',sel_patients], '|continuous', sep='')
sel_patients = c('H132384','H134819')
dplot_table['TP53_pathway 0',sel_patients] = paste(dplot_table['TP53_pathway 0',sel_patients], '|continuous', sep='')
dplot_table['Cell_cycle 0',sel_patients] = paste(dplot_table['Cell_cycle 0',sel_patients], '|continuous', sep='')
# Fix the patients with parallel acquisition of the same CN changes
sel_patients = c('H132383')
dplot_table['1q_gain 0',sel_patients] = paste(dplot_table['1q_gain 0',sel_patients], '|parallel', sep='')
sel_patients = c('H108313')
dplot_table['7q_gain 0',sel_patients] = paste(dplot_table['7q_gain 0',sel_patients], '|parallel', sep='')
sel_patients = c('H118706','H132392')
dplot_table['4p_loss 0',sel_patients] = paste(dplot_table['4p_loss 0',sel_patients], '|parallel', sep='')
sel_patients = c('H103207')
dplot_table['17q_gain 0',sel_patients] = paste(dplot_table['17q_gain 0',sel_patients], '|parallel', sep='')
sel_patients = c('H112909','H134822','H134819','H135088')
dplot_table['3p_loss 0',sel_patients] = paste(dplot_table['3p_loss 0',sel_patients], '|parallel', sep='')
sel_patients = c('H134819')
dplot_table['11q_loss 0',sel_patients] = paste(dplot_table['11q_loss 0',sel_patients], '|parallel', sep='')
sel_patients = c('H132396')
dplot_table['12q_gain 0',sel_patients] = paste(dplot_table['12q_gain 0',sel_patients], '|parallel', sep='')
# Heatmap alteration function
alter_fun = list(
	background = function(x, y, w, h) 
		grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#BEBEBE", col = NA)),
	# n=1 mutation
	'1 event' = function(x, y, w, h) 
		grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col['1 event'], col = NA)),
	# m=2 mutations
	'2 events' = function(x, y, w, h) 
		grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col['2 events'], col = NA)),
	# n>=3 mutations
	'>=3 events' = function(x, y, w, h)
		grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col['>=3 events'], col = NA)),
	# parallel events
	'parallel' = function(x, y, w, h) 
		grid.points(x, y, pch = 16, size = unit(5,"points")),
	# continuous events
	continuous = function(x, y, w, h) {
		grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(lwd = 1))
		grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(lwd = 1))
	}
)
rownames(dplot_table) = colsplit(rownames(dplot_table), pattern=" ", names=1:2)[,1]
rownames(dplot_table) = gsub('_',' ',sub('_cnloh','=',sub('_loss','-',sub('_gain','+',rownames(dplot_table)))))
rnames = rownames(dplot_table)
dplot_table = matrix(dplot_table,nr=nrow(dplot_table))
colnames(dplot_table) = porder
rownames(dplot_table) = rnames
pdf('raw/figure_3b_pathways.pdf', width=15, height=10)
oncoPrint(
	dplot_table, alter_fun = alter_fun, col = col, show_pct=FALSE, top_annotation=NULL,
	show_row_names=TRUE, show_column_names=TRUE, column_order=colnames(dplot_table),
	column_names_gp=gpar(fontsize=10, fontfamily='sans'), row_names_gp=gpar(fontsize=10, fontfamily='sans'),
        row_split=factor(c(
                rep("CNA on trunk",length(trunk_cna_order)), rep("mutation on trunk",length(trunk_genes_order)),
                rep("CNA off-trunk",length(offtrunk_cna_order)), rep("mutation off-trunk",length(offtrunk_genes_order))
                ),
                levels=c('CNA on trunk','mutation on trunk','CNA off-trunk','mutation off-trunk')
        ),
        column_split=factor(
                c(rep("MYCN",length(patients_mycn)), rep("TERT",length(patients_tert)), rep("ATRX",length(patients_atrx)),
                rep('MDM2\nCDK4',length(patients_mdm2_ckd4)), rep('Seg_CNV',length(patients_seg_cnv)), rep('Num_CNV',length(patients_num_cnv))),
                levels=c('MYCN','ATRX','TERT','MDM2\nCDK4','Seg_CNV','Num_CNV')
        ),
	bottom_annotation=bottom_ha_patients
)
dev.off()
