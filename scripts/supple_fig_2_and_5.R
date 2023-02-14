# Read all subclone data
WDIR = '~/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/compare_clones.R')
no_wgs_tumors = as.data.frame(table(mdata$PATIENT[mdata$SEQ=='WGS']))
colnames(no_wgs_tumors) = c('patient','no_wgs_tumor_samples')
subs_clones = merge(x=subs_clones, y=no_wgs_tumors, by.x='patient', by.y='patient')
# Compare substitution clusters to indel clusters

# Compare CCFs of matched clusters
merged_ccfs = unique(subs_clones[,c('patient','clone','sample','descr','ccf','matching_indel_clone','matching_indel')])
colnames(merged_ccfs)[5] = 'subs_ccf'
merged_ccfs$key = paste(merged_ccfs$sample, merged_ccfs$matching_indel_clone)
indel_clones$key = paste(indel_clones$sample, indel_clones$clone)
merged_ccfs = merge(x=merged_ccfs, y=indel_clones[,c('key','ccf')], by.x='key', by.y='key')
colnames(merged_ccfs)[ncol(merged_ccfs)] = 'indel_ccf'
plot_ccf_scatter = ggscatter(merged_ccfs, x='subs_ccf', y='indel_ccf', add="reg.line", conf.int=TRUE, add.params=list(color="blue", fill="lightgray"), col='darkgrey') + stat_cor(method="pearson", label.x=0.1, label.y=0.8)

trunk_ccfs = melt(merged_ccfs[merged_ccfs$descr=='trunk',c('key','subs_ccf','indel_ccf')], id.vars=c('key'))
trunk_ccfs$variable = sub('_ccf','',trunk_ccfs$variable)
plot_boxplot_trunkal_ccf = ggboxplot(trunk_ccfs, x='variable', y='value', ylab='CCF of trunk', col='darkgrey', add='jitter') + stat_compare_means(comparisons=list(c('subs','indel')))

plot_boxplot_nomuts_matched_subs_clusters = ggboxplot(unique(subs_clones[,c('clone','patient','no.muts.in.cluster','matching_indel')]), x='matching_indel', y='no.muts.in.cluster', col='matching_indel', add='jitter', palette=c("#E7B800","#00AFBB"))
plot_boxplot_nomuts_matched_subs_clusters = plot_boxplot_nomuts_matched_subs_clusters + stat_compare_means(comparisons=list(c('yes','no'))) + rremove('legend')
plot_boxplot_nomuts_matched_indel_clusters = ggboxplot(unique(indel_clones[,c('clone','patient','no.muts.in.cluster','matching_subs')]), x='matching_subs', y='no.muts.in.cluster', col='matching_subs', add='jitter', palette=c("#E7B800","#00AFBB"))
plot_boxplot_nomuts_matched_indel_clusters = plot_boxplot_nomuts_matched_indel_clusters + stat_compare_means(comparisons=list(c('yes','no'))) + rremove('legend')

plot_density_subs_ccfs_unmatched_clusters =  ggdensity(unique(subs_clones[,c('clone','patient','ccf','matching_indel')]), x='ccf', color='matching_indel', fill='matching_indel', palette=c("#00AFBB", "#E7B800"), xlab="CCF in SNV cluster")
plot_density_indel_ccfs_unmatched_clusters =  ggdensity(unique(indel_clones[,c('clone','patient','ccf','matching_subs')]), x='ccf', color='matching_subs', fill='matching_subs', palette=c("#00AFBB", "#E7B800"), xlab='CCF in Indel cluster')

# Determine the number of samples substitution clones are present in
subs_clones_present_in_gt_1_sample = names(which(table(unique(subs_clones[subs_clones$ccf>0.05,c('clone','sample')])$clone)>1))
sclones = as.data.frame(unique(subs_clones[,c('clone','matching_indel')]))
sclones$type = 'sample_specific'
sclones$type[sclones$clone %in% subs_clones_present_in_gt_1_sample] = 'shared'
bdata1 = as.data.frame(table(sclones$matching_indel,sclones$type))
colnames(bdata1) = c('matching_indel_clone','subs_clone_type','number')
bdata1$proportion = NA
bdata1[bdata1$subs_clone_type=='sample_specific','proportion'] = bdata1[bdata1$subs_clone_type=='sample_specific','number'] / sum(bdata1[bdata1$subs_clone_type=='sample_specific','number'])
bdata1[bdata1$subs_clone_type=='shared','proportion'] = bdata1[bdata1$subs_clone_type=='shared','number'] / sum(bdata1[bdata1$subs_clone_type=='shared','number'])
# Determine the number samples indel clones are present in
indel_clones_present_in_gt_1_sample = names(which(table(unique(indel_clones[indel_clones$ccf>0.05,c('clone','sample')])$clone)>1))
inclones = as.data.frame(unique(indel_clones[,c('clone','matching_subs')]))
inclones$type = 'sample_specific'
inclones$type[inclones$clone %in% indel_clones_present_in_gt_1_sample] = 'shared'
bdata2 = as.data.frame(table(inclones$matching_subs,inclones$type))
colnames(bdata2) = c('matching_subs_clone','indel_clone_type','number')
bdata2[bdata2$indel_clone_type=='sample_specific','proportion'] = bdata2[bdata2$indel_clone_type=='sample_specific','number'] / sum(bdata2[bdata2$indel_clone_type=='sample_specific','number'])
bdata2[bdata2$indel_clone_type=='shared','proportion'] = bdata2[bdata2$indel_clone_type=='shared','number'] / sum(bdata2[bdata2$indel_clone_type=='shared','number'])

plot_number_of_samples_subs = ggbarplot(bdata1, x='subs_clone_type', y='proportion', fill='matching_indel_clone', palette=c("#00AFBB", "#E7B800"))
plot_number_of_samples_indels = ggbarplot(bdata2, x='indel_clone_type', y='proportion', fill='matching_subs_clone', palette=c("#00AFBB", "#E7B800"))

pdf('raw/supple_fig_5.pdf', width=8, height=16)
plot_grid(
	plot_density_subs_ccfs_unmatched_clusters, plot_density_indel_ccfs_unmatched_clusters,
	plot_boxplot_nomuts_matched_subs_clusters, plot_boxplot_nomuts_matched_indel_clusters,
	plot_number_of_samples_subs, plot_number_of_samples_indels, plot_ccf_scatter, plot_boxplot_trunkal_ccf, nrow=4
)
dev.off()

### Read data on evolutionary patterns
evolution = read.delim('data/evolutionary_patterns.txt', header=T, sep="\t", stringsAsFactors=F)
evolution$spread[is.na(evolution$spread)] = 'only_primary'
evolution$spread[evolution$spread=='D'] = 'metastasis'
evolution$spread[evolution$spread=='L'] = 'locoregional'
evolution$spread = factor(evolution$spread, levels=c('only_primary','locoregional','metastasis'))
evolution$stage[evolution$stage!=4] = 'low-stage'
evolution$stage[evolution$stage==4] = 'stage-4'
evolution$stage = factor(evolution$stage, levels=c('low-stage','stage-4'))
### Number and proportion of mutations on the trunk
prop_subs_on_trunk = unique(subs_clones[subs_clones$descr=='trunk',c('patient','prop.muts.in.cluster','no.muts.in.cluster','no_wgs_tumor_samples')])
colnames(prop_subs_on_trunk)[2:4] = c('prop_subs_on_trunk','no_subs_on_trunk','no_sample')
prop_subs_on_trunk = merge(x=prop_subs_on_trunk, y=unique(mdata[,c('PATIENT','DISEASE_SUBTYPE','AADx','AGE_GROUP')]), by.x='patient', by.y='PATIENT')
prop_indels_on_trunk = unique(indel_clones[indel_clones$descr=='trunk',c('patient','prop.muts.in.cluster')])
colnames(prop_indels_on_trunk)[2] = 'prop_indels_on_trunk'
prop_muts_on_trunk = merge(x=prop_subs_on_trunk, y=prop_indels_on_trunk, by.x='patient', by.y='patient')
prop_muts_on_trunk = merge(x=prop_muts_on_trunk, y=as.data.frame(table(unique(subs_clones[,c('patient','clone')])$patient)), by.x='patient', by.y='Var1')
colnames(prop_muts_on_trunk)[ncol(prop_muts_on_trunk)] = 'no_subclones'
prop_muts_on_trunk$no_subclones_per_sample = prop_muts_on_trunk$no_subclones / prop_muts_on_trunk$no_sample
prop_muts_on_trunk$no_sample[prop_muts_on_trunk$no_sample>=4] = "4+"
prop_muts_on_trunk$no_sample = factor(prop_muts_on_trunk$no_sample, levels=c('2','3','4+'))

# Check the relationship between subclonal structure and the number of WGS samples
symnum_args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
my_comparisons1 <- list(c("2","3"), c("2","4+"), c("3", "4+"))
plot_nosubclone_vs_nosample <- ggboxplot(prop_muts_on_trunk, x='no_sample', y='no_subclones', add='jitter', fill='no_sample', palette=get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 4))
plot_nosubclone_vs_nosample <- plot_nosubclone_vs_nosample + stat_compare_means(comparisons=my_comparisons1, symnum.args=symnum_args) + rremove("legend")

# Number of mutations on trunk vs no of tumors
plot_nosubs_vs_nosample <- ggboxplot(prop_muts_on_trunk, x='no_sample', y='no_subs_on_trunk', add='jitter', fill='no_sample', palette=get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 4))
plot_nosubs_vs_nosample <- plot_nosubs_vs_nosample + stat_compare_means(comparisons = my_comparisons1, symnum.args=symnum_args) + rremove("legend")

# Check the relationship between subclonal structure and disease subtype
my_comparisons2 <- list( c("MYCN", "TERT"), c("MYCN","ATRX"), c('MYCN','SEG_CNV'), c("TERT","ATRX"), c("TERT", "SEG_CNV"), c("ATRX","SEG_CNV"))
plot_nosubclone_vs_dtype <- ggboxplot(
        prop_muts_on_trunk[prop_muts_on_trunk$no_sample!='1',], x='DISEASE_SUBTYPE', y='no_subclones_per_sample', add='jitter', fill='DISEASE_SUBTYPE',
        xlab="", legend="none", palette=disease_subtype_cols$color
)
plot_nosubclone_vs_dtype <- plot_nosubclone_vs_dtype + stat_compare_means(comparisons=my_comparisons2, symnum.args=symnum_args) + rotate_x_text(angle=45)
plot_nosubs_vs_dtype <- ggboxplot(
        prop_muts_on_trunk[prop_muts_on_trunk$no_sample!='1',], x='DISEASE_SUBTYPE', y='no_subs_on_trunk', add='jitter', fill='DISEASE_SUBTYPE',
        xlab="", legend="none", palette=disease_subtype_cols$color
)
plot_nosubs_vs_dtype <- plot_nosubs_vs_dtype + stat_compare_means(comparisons=my_comparisons2, symnum.args=symnum_args) + rotate_x_text(angle=45)
# Check the relationship between subclonal structure and type of spread
prop_muts_on_trunk_evolution = merge(x=prop_muts_on_trunk, y=evolution[,c('Patient','spread','stage','X.Pri')], by.x='patient', by.y='Patient')
my_comparisons3 <- list( c("only_primary", "locoregional"), c("only_primary","metastasis"), c('locoregional','metastasis'))
plot_nosubs_vs_spread <- ggboxplot(prop_muts_on_trunk_evolution, x='spread', y='no_subs_on_trunk', fill='spread', add='jitter', xlab='', palette=c('#4A708B','#AA4400','#800000'))
plot_nosubs_vs_spread <- plot_nosubs_vs_spread + stat_compare_means(comparisons=my_comparisons3, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(angle=45)
plot_nosubclone_vs_spread <- ggboxplot(prop_muts_on_trunk_evolution, x='spread', y='no_subclones_per_sample', fill='spread', add='jitter', xlab='', palette=c('#4A708B','#AA4400','#800000'))
plot_nosubclone_vs_spread <- plot_nosubclone_vs_spread + stat_compare_means(comparisons=my_comparisons3, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(angle=45)

# Check the relationship between subclonal structure and stage
my_comparisons4 <- list( c("low-stage", "stage-4"))
plot_nosubclone_vs_stage <- ggboxplot(prop_muts_on_trunk_evolution[prop_muts_on_trunk_evolution$no_sample!='1',], x='stage', y='no_subclones_per_sample', fill='stage', add='jitter', xlab='', palette=c('#C3CCCC','#6A8181'))
plot_nosubclone_vs_stage <- plot_nosubclone_vs_stage + stat_compare_means(comparisons=my_comparisons4, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(angle=45)
plot_nosubs_vs_stage <- ggboxplot(prop_muts_on_trunk_evolution[prop_muts_on_trunk_evolution$no_sample!='1',], x='stage', y='no_subs_on_trunk', fill='stage', add='jitter', xlab='', palette=c('#C3CCCC','#6A8181'))
plot_nosubs_vs_stage <- plot_nosubs_vs_stage + stat_compare_means(comparisons=my_comparisons4, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(angle=45)


pdf('raw/supple_fig_2.pdf', width=8, height=16)
plot_grid(plot_nosubclone_vs_nosample, plot_nosubs_vs_nosample, plot_nosubclone_vs_dtype, plot_nosubs_vs_dtype, plot_nosubclone_vs_spread, plot_nosubs_vs_spread, plot_nosubclone_vs_stage, plot_nosubs_vs_stage, nrow=4)
dev.off()
