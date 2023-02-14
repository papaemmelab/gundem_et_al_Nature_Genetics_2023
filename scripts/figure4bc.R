library(plyr)
library(ggpubr)
library(cowplot)
library(plyr)
library(reshape2)
library(pals)

evolution  = read.delim('analysis/revision1/data/evolutionary_patterns.txt', header=T, sep="\t", stringsAsFactors=F)
relapses = read.table('~/Documents/projects/neuroblastoma/triple_callers/analysis/revision1/data/evolutionary_patterns_relapses.txt', header=T, sep="\t", stringsAsFactors=F)
relapses = merge(x=relapses, y=evolution[,c('Patient','main_driver')], by.x='Patient', by.y='Patient')
relapses$id = rownames(relapses)
# Plot the number of clonal transition between dx-relapse and relapses
barplot_transitions_reltype <- ggbarplot(as.data.frame(table(relapses$relapse_type, relapses$relapse_mode)), x='Var1', y='Freq', fill='Var2', palette=c('#875692','#F2F3F4','#222222','#F3C300'), legend='none') + rotate_x_text(angle=45)
# Plot the number clonal transition per subtype
barplot_transitions_subtypes <- ggbarplot(as.data.frame(table(relapses$main_driver, relapses$relapse_mode)), x='Var1', y='Freq', fill='Var2', palette=c('#875692','#F2F3F4','#222222','#F3C300'), legend='right') + rotate_x_text(angle=45)
cnas = c('1p=','1p-','1q+','2p+','2q-','3p-','3p+','4-','4p-','4q-','5-','5q-','6p-','6p+','6q-','7-','7p-','7q+','7qter+','8p-','8q-','8q+','9=','9p-','9p=','9p+','9q+','11q-','12q+','14q-','14=','16p-','16p=','16q-','17p-','17q+','19+','19p-','19q-','21q-','21q+','chromothripsis_20','chromothripsis_16','diploid','tetraploid','Xp-','10qter-','APC','del_7p','11pter+','del_5q')
################################
temporal_info = NULL
for(i in 1:nrow(relapses)) {
	temporal_info = rbind(temporal_info, cbind(
		relapses$id[i], relapses$Patient[i], relapses$main_driver[i], relapses$relapse_mode[i], relapses$relapse_type[i],
		relapses$driver[i], unlist(strsplit(relapses$relapse_driver_to[i], split=","))
	))
}
temporal_info = as.data.frame(temporal_info)
colnames(temporal_info) = c('id','patient','disease_subtype','relapse_mode','relapse_type','driver','relapse_gene')
temporal_info$relapse_pathway = NA
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ALK','HRAS','NF1','BRAF','SOS1','CIC','FGFR1','NRAS','KRAS','MAP3K1')] = 'RTK_RAS_MAPK'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('PTPRD','DLG2','TIAM1','AUTS2','SHANK2')] = 'NEURODEVELOPMENT'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ARID1B','ARID1A')] = 'ARID1A/B'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ATRX')] = 'ATRX'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('TERT')] = 'TERT'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('PIK3CA','PIK3R1','MTOR','FBXW7')] = 'PI3K_MTOR'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('CDKN2A','CDKN1C')] = 'CELL_CYCLE'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('CREBBP','SMARCA4','EP300','BCOR')] = 'CHROMATIN'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('MYCN','MYCL1')] = 'MYCN'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c(c('RARA','RARB','PHOX2B','WNT5A','IGF2BP3','SPRY2','KLF4'))] = 'DIFFERENTIATION'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('17p-','TP53')] = 'TP53/17p-'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% cnas] = 'CNA'
temporal_info[is.na(temporal_info$relapse_pathway) & !is.na(temporal_info$relapse_gene),'relapse_pathway'] = 'other'
temporal_info_to = temporal_info
colnames(temporal_info_to)[7:8] = c('relapse_gene_to','relapse_pathway_to')
################################
temporal_info = NULL
for(i in 1:nrow(relapses)) {
	temporal_info = rbind(temporal_info, cbind(
		relapses$id[i], relapses$Patient[i], relapses$main_driver[i], relapses$relapse_mode[i], relapses$relapse_type[i],
		relapses$driver[i], unlist(strsplit(relapses$relapse_driver_from[i], split=","))
	))
}
temporal_info = as.data.frame(temporal_info)
colnames(temporal_info) = c('id','patient','disease_subtype','relapse_mode','relapse_type','driver','relapse_gene')
temporal_info$relapse_pathway = NA
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ALK','HRAS','NF1','BRAF','SOS1','CIC','FGFR1','NRAS','KRAS','MAP3K1')] = 'RTK_RAS_MAPK'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('PTPRD','DLG2','TIAM1','AUTS2','SHANK2')] = 'NEURODEVELOPMENT'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ARID1B','ARID1A')] = 'ARID1A/B'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('ATRX')] = 'ATRX'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('TERT')] = 'TERT'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('PIK3CA','PIK3R1','MTOR','FBXW7')] = 'PI3K_MTOR'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('CDKN2A','CDKN1C')] = 'CELL_CYCLE'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('CREBBP','SMARCA4','EP300','BCOR')] = 'CHROMATIN'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('MYCN','MYCL1')] = 'MYCN'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c(c('RARA','RARB','PHOX2B','WNT5A','IGF2BP3','SPRY2','KLF4'))] = 'DIFFERENTIATION'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% c('17p-','TP53')] = 'TP53/17p-'
temporal_info$relapse_pathway[temporal_info$relapse_gene %in% cnas] = 'CNA'
temporal_info[is.na(temporal_info$relapse_pathway) & !is.na(temporal_info$relapse_gene),'relapse_pathway'] = 'other'
temporal_info_from = temporal_info
colnames(temporal_info_from)[7:8] = c('relapse_gene_from','relapse_pathway_from')
######################################
clonal_transitions = merge(x=temporal_info_from, y=temporal_info_to[,c('id','relapse_gene_to','relapse_pathway_to')], by.x='id', by.y='id')
clonal_transitions_emergent = unique(clonal_transitions[clonal_transitions$relapse_mode=='emergent',c('id','patient','relapse_pathway_from','relapse_pathway_to')])
clonal_transitions_emergent[is.na(clonal_transitions_emergent$relapse_pathway_to),'relapse_pathway_to'] = 'CNA'
clonal_transitions_switch = unique(clonal_transitions[clonal_transitions$relapse_mode=='switch',c('id','patient','relapse_pathway_from','relapse_pathway_to')])
clonal_transitions_switch[is.na(clonal_transitions_switch$relapse_pathway_to),'relapse_pathway_to'] = 'CNA'

tr_switch = cbind('switch',unique(clonal_transitions_switch[,c('id','patient','relapse_pathway_to')]))
tr_emergent = cbind('emergent', unique(clonal_transitions_emergent[,c('id','relapse_pathway_to')]))
colnames(tr_switch)[1] = 'type'
colnames(tr_emergent)[1] = 'type'
tr_summary = ddply(rbind(tr_switch[,c('type','id','relapse_pathway_to')], tr_emergent[,c('type','id','relapse_pathway_to')]), c('relapse_pathway_to','type'), summarise, n=length(id))
tr_total = ddply(rbind(tr_switch[,c('type','id','relapse_pathway_to')], tr_emergent[,c('type','id','relapse_pathway_to')]), c('relapse_pathway_to'), summarise, n=length(id))
tr_summary$relapse_pathway_to = factor(tr_summary$relapse_pathway_to, levels=tr_total$relapse_pathway_to[order(tr_total$n, decreasing=T)])
# Plot the number of clonal transitions affecting each pathway
barplot_transitions_pathways <- ggbarplot(tr_summary, x='relapse_pathway_to', y='n', fill='type', palette=c('#F3C300','#875692'), legend='none') + rotate_x_text(angle=45)

pdf('~/Documents/projects/neuroblastoma/triple_callers/analysis/revision1/raw/figure4bc.pdf', height=4, width=12)
plot_grid(barplot_transitions_reltype, barplot_transitions_subtypes, barplot_transitions_pathways, nrow=1, align='h', rel_widths=c(0.7,2,2.5))
dev.off()

#################################################
clonal_transitions_switch$is_same = 'yes'
clonal_transitions_switch[which(clonal_transitions_switch$relapse_pathway_from!=clonal_transitions_switch$relapse_pathway_to),'is_same'] = 'no'
clonal_transitions_switch$summary = apply(clonal_transitions_switch[,c('relapse_pathway_from','relapse_pathway_to')], 1, paste, collapse=" > ")
barplot_clonal_switches <- ggbarplot(as.data.frame(sort(table(clonal_transitions_switch$summary[!is.na(clonal_transitions_switch$relapse_pathway_from) & clonal_transitions_switch$is_same=='yes']), decreasing=T)), x='Var1', y='Freq', fill='grey', ylab='n', xlab='clonal switch') + rotate_x_text(angle=45)


plot_grid(barplot_transitions_reltype, barplot_transitions_subtypes, barplot_transitions_pathways, barplot_clonal_switches, nrow=1, align='h', rel_widths=c(0.7,2,2.5,2))

#plot_grid(plot_relapse_time, plot_subclone, plot_relapse_modes, plot_relapse_drivers, nrow=1, rel_widths=c(1.2,1.2,1,3.5))

pdf('~/Documents/projects/neuroblastoma/triple_callers/analysis/revision1/raw/figure4c_inside.pdf', height=4, width=4)
print(barplot_clonal_switches)
dev.off()
##########################################################
btw_relapses = evolution[,c('Patient','main_driver','relapse_analysis','relapse_mode_btw_relapses','relapse_driver_btw_relapses_gene_to')]
colnames(btw_relapses)[1] = 'relapse_type'
colnames(btw_relapses)[4] = 'relapse_mode'
colnames(btw_relapses)[5] = 'relapse_gene_to'
btw_relapse_dx = evolution[,c('Patient','main_driver','relapse_analysis','relapse_mode_btw_dx_and_relapse','relapse_driver_btw_dx_relapse_gene_to')]
colnames(btw_relapse_dx)[1] = 'relapse_type'
colnames(btw_relapse_dx)[4] = 'relapse_mode'
colnames(btw_relapse_dx)[5] = 'relapse_gene_to'
relapse_info = rbind(btw_relapses, btw_relapse_dx)
relapse_info = relapse_info[!is.na(relapse_info$relapse_mode) & relapse_info$relapse_mode!="",]

ggbarplot(as.data.frame(table(relapse_info$relapse_analysis, relapse_info$relapse_mode) / rowSums(table(relapse_info$relapse_analysis, relapse_info$relapse_mode))), x='Var1', y='Freq', fill='Var2') + rotate_x_text(angle=45)

temporal_info = NULL
for(i in 1:nrow(temporal)) {
	temporal_info = rbind(temporal_info, cbind(
		temporal$Patient[i], temporal$main_driver[i], temporal$relapse_analysis[i], temporal$relapse_mode[i],
		unlist(strsplit(temporal[temporal$Patient==temporal$Patient[i],'relapse_driver_gene_to'], split=","))
	))
}
temporal_info = as.data.frame(temporal_info)
colnames(temporal_info) = c('Patient','disease_subtype','relapse_analysis','relapse_mode','relapse_driver_gene_to')
temporal_info$relapse_analysis = as.character(temporal_info$relapse_analysis)
temporal_info$relapse_analysis[grep('relapses',temporal_info$relapse_analysis)] = 'btw_relapses'

ggbarplot(as.data.frame(table(temporal_info$relapse_analysis, temporal_info$relapse_mode) / rowSums(table(temporal_info$relapse_analysis, temporal_info$relapse_mode))), x='Var1', y='Freq', fill='Var2') + rotate_x_text(angle=45)
