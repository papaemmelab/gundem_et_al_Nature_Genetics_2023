library(ggpubr)
library(cowplot)

muts = read.table('data/mutations_for_matched_wgs_and_impact_samples.txt', header=T, sep="\t", stringsAsFactors=F)

plot1 <- ggscatter(muts, x='wgs_vaf', y='impact_vaf', col='call_status', xlab='WGS VAF', ylab='MSK-IMPACT VAF')
plot2 <- ggpaired(muts, cond1='impact_depth', cond2='wgs_depth', col='call_status', ylab='Coverage', xlab='', legend='none', line.color='grey')
plot2 <- facet(plot2, facet.by='call_status')
pdf('raw/supple_fig_3.pdf', width=12)
plot_grid(plot1, plot2, nrow=1)
dev.off()


