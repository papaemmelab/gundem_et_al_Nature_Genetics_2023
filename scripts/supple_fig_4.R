plot_genomewide_cn_landscape <- function(bbg_segments, title="", cn_ymin_high=(-1), cn_ymax_high=1, cn_ymin=(-1), cn_ymax=1) {
	total_no_samples = length(unique(bbg_segments$sample))
	cGain.segs = NULL
	cLoss.segs = NULL
	for(chr in unique(bbg_segments$chr)) {
		datasub = bbg_segments[bbg_segments$chr==chr,]
		all.positions = sort(unique(c(datasub$startpos, datasub$endpos)))
		all.positions = cbind(all.positions[1:length(all.positions)-1], all.positions[2:length(all.positions)])
		colnames(all.positions) = c('start','end')
		for(i in 1:nrow(all.positions)) {
			sel = datasub$startpos<=all.positions[i,'start'] & datasub$endpos>=all.positions[i,'end']
			clonal.cn = datasub[datasub$frac1_A==1 & sel,c('sample','nMaj1_A','nMin1_A','ntot')]
			n.gain = 0
			n.loss = 0
			if (nrow(clonal.cn)>0) {
				n.gain = length(which(rowSums(clonal.cn[,c('nMaj1_A','nMin1_A')]) > sample_info[clonal.cn$sample,'ploidy']))
				n.loss = length(which(clonal.cn$nMin1_A==0))
			}
			cGain.segs = rbind(cGain.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.gain))
			cLoss.segs = rbind(cLoss.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.loss))
		}
	}
	colnames(cGain.segs) = c('chr','sp','ep','val')
	cGain.segs = as.data.frame(cGain.segs)
	cGain.segs$ep = as.integer(as.character(cGain.segs$ep))
	cGain.segs$sp = as.integer(as.character(cGain.segs$sp))
	cGain.segs$val = as.integer(as.character(cGain.segs$val))
	colnames(cLoss.segs) = c('chr','sp','ep','val')
	cLoss.segs = as.data.frame(cLoss.segs)
	cLoss.segs$ep = as.integer(as.character(cLoss.segs$ep))
	cLoss.segs$sp = as.integer(as.character(cLoss.segs$sp))
	cLoss.segs$val = as.integer(as.character(cLoss.segs$val))
	# Chromosomes
	chrs = read.table(file='data/gr37.fasta.fai', header=F, stringsAsFactors=F)
	maxchr = chrs$V2
	names(maxchr) = chrs$V1
	loci = rep(0,23)
	names(loci) = names(maxchr[1:23])
	for (i in 1:(length(loci)-1)) {
		chr = names(loci)[i]
		loci[i+1] = loci[i] + maxchr[chr]
	}
	chrg = NULL
	for (i in 2:length(loci)) {
		chrg = c(chrg, (loci[i]+loci[i-1])/2)
	}
	# Plot
	clonal.data = cGain.segs[!(cGain.segs$chr %in% c('23','X')),]
	clonal.data$chrlen = rep(0,nrow(clonal.data))
	for(chr in unique(clonal.data$chr)) {clonal.data[clonal.data$chr==chr,'chrlen'] = loci[chr]}
	clonal.data[,'sp'] = clonal.data[,'sp'] + clonal.data[,'chrlen']
	clonal.data[,'ep'] = clonal.data[,'ep'] + clonal.data[,'chrlen']
	plot(c(), ylim=c(cn_ymin,cn_ymax), xlim=c(1,max(clonal.data$ep)), bty="n", xlab="", ylab="", xaxt="n", main=title)
	segments(x0=clonal.data[,'sp'], x1=clonal.data[,'ep'], y0=clonal.data[,'val']/total_no_samples, lwd=2, col='orangered3')

	clonal.data = cLoss.segs[!(cLoss.segs$chr %in% c('23','X')),]
	clonal.data$chrlen = rep(0,nrow(clonal.data))
	for(chr in unique(clonal.data$chr)) {clonal.data[clonal.data$chr==chr,'chrlen'] = loci[chr]}
	clonal.data[,'sp'] = clonal.data[,'sp'] + clonal.data[,'chrlen']
	clonal.data[,'ep'] = clonal.data[,'ep'] + clonal.data[,'chrlen']
	segments(x0=clonal.data[,'sp'], x1=clonal.data[,'ep'], y0=(clonal.data[,'val']/(-total_no_samples)), lwd=2, col='navyblue')
	abline(v=loci, lty = 2)
	text(x = chrg, y = rep(1, 22), labels=1:22)
	legend(0,0.9, horiz=T, lty=c(1,1), lwd=c(2.5,2.5), c('gain','loss'), col = c('orangered3','navyblue'), bg="white")
}
### Read FACET outputs
sample_info = read.table('data/impact_samples.txt', header=F, stringsAsFactors=F)
sample_info = sample_info[!(sample_info$V1 %in% c('P-0016160-T07-IH3','P-0039609-T01-IM6','P-0041009-T02-IM6','P-0058028-T01-IM6')),]
colnames(sample_info) = c('sample','fitname','purity','ploidy')
all_facet = NULL
for(i in 1:nrow(sample_info)) {
	facet_dir = grep(sample_info$fitname[i], list.files(list.files('~/Documents/projects/neuroblastoma/triple_callers/impact_data', pattern=sample_info$sample[i], full.names=T), full.names=T), value=T)
	facet_file = list.files(facet_dir, pattern="_purity.cncf.txt", full.names=T)
	facet = read.table(facet_file, header=T, sep="\t", stringsAsFactors=F)
	facet$sample = sample_info$sample[i]
	all_facet = rbind(all_facet, facet)
}
colnames(sample_info) = c('sample','fitname','purity','ploidy')
rownames(sample_info) = sample_info$sample
colnames(all_facet)[colnames(all_facet)=="chrom"] = 'chr'
colnames(all_facet)[colnames(all_facet)=="loc.start"] = 'startpos'
colnames(all_facet)[colnames(all_facet)=="loc.end"] = 'endpos'
colnames(all_facet)[colnames(all_facet)=="tcn.em"] = 'ntot'
colnames(all_facet)[colnames(all_facet)=="lcn.em"] = 'nMin1_A'
all_facet$nMaj1_A = all_facet$ntot - all_facet$nMin1_A
all_facet$frac1_A = 1 # Assume all CNAs are clonal

pdf('raw/supple_fig_4.pdf', width=18, height=12)
par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
plot_genomewide_cn_landscape(all_facet, title=paste("CNAs from MSK-IMPACT n=", length(unique(all_facet$sample)), ' tumors', sep=""))

### Read Battenberg output
sample_info = read.table('data/wgs_samples.txt', header=T, stringsAsFactors=F)
all_bbg = NULL
for(i in 1:nrow(sample_info)) {
	bbg_file = list.files(paste('~/Documents/projects/neuroblastoma/triple_callers/', sample_info$sample[i], sep=''), pattern="*subclones_raw.txt", full.names=T, recursive=T)
	bbg = read.table(bbg_file, header=T, sep="\t", stringsAsFactors=F)
	bbg$sample = sample_info$sample[i]
	all_bbg = rbind(all_bbg, bbg)
}
rownames(sample_info) = sample_info$sample
all_bbg$frac1_A = 1

plot_genomewide_cn_landscape(all_bbg, title=paste("CNAs from WGS n=", length(unique(all_bbg$sample)), ' tumors', sep=""))
dev.off()
