library(ggpubr)

cosine.similarity <- function(A,B){
	return(sum(A*B)/sqrt(sum(A^2)*sum(B^2)))
}

find_matching_clones <- function(clones1, clones2, patient, min_cos_sim=0.7, max_dif=0.15, min_ccf_to_call_clone=0.05) {
	# Find the mtaching clones from clones1 to clones2
	# Calculate cosine similarity amongst all substitution and indel clones in the patient
	ccfs1 = dcast(unique(clones1[clones1$patient==patient,c('sample','clone','ccf')]), clone ~ sample, value.var="ccf")
	ccfs2 = dcast(unique(clones2[clones2$patient==patient,c('sample','clone','ccf')]), clone ~ sample, value.var="ccf")
	sel_samples = intersect(colnames(ccfs1), colnames(ccfs2))
	ccfs1 = ccfs1[,sel_samples]
	ccfs2 = ccfs2[,sel_samples]
	rownames(ccfs2) = ccfs2$clone
	rownames(ccfs1) = ccfs1$clone
	ccfs2$clone = NULL
	ccfs1$clone = NULL
	ccfs1 = as.matrix(ccfs1)
	ccfs2 = as.matrix(ccfs2)
	ccfs1_rownames = rownames(ccfs1)
	ccfs1 = apply(ccfs1, 2, as.numeric)
	if(length(ccfs1_rownames)==1) {
		ccfs1 = as.matrix(t(ccfs1))
		rownames(ccfs1) = ccfs1_rownames
	} else {
		rownames(ccfs1) = ccfs1_rownames
	}
	rownames(ccfs1) = ccfs1_rownames
	ccfs2_rownames = rownames(ccfs2)
	ccfs2 = apply(ccfs2, 2, as.numeric)
	if(length(ccfs2_rownames)==1) {
		ccfs2 = as.matrix(t(ccfs2))
		rownames(ccfs2) = ccfs2_rownames
	} else {
		rownames(ccfs2) = ccfs2_rownames
	}
	cos_sims = mat.or.vec(nrow(ccfs1), nrow(ccfs2))
	if(!('matrix' %in% class(cos_sims))) {cos_sims = as.matrix(cos_sims)}
	for(k in 1:nrow(cos_sims)) {
		for(m in 1:ncol(cos_sims)) {
			cos_sims[k,m] = cosine.similarity(ccfs1[k,], ccfs2[m,])
		}
	}
	colnames(cos_sims) = rownames(ccfs2)
	rownames(cos_sims) = rownames(ccfs1)
	cos_sims[is.na(cos_sims)] = 0
	cos_sims[cos_sims<=min_cos_sim] = 0
	if(nrow(cos_sims)>=2 & ncol(cos_sims)>=2) {
		heatmap.plus(
			cos_sims, scale="none", labCol=sub(paste(patient,'__',sep=""),'',rownames(ccfs2)),
			labRow=sub(paste(patient,'__',sep=""),'',rownames(ccfs1)), margin=c(15,15),
			ylab="Clones1", xlab="Clones2"
		)
	}
	# Compare only those clones with cosine similarity > 0.9
	matching_clones = NULL
	if(ncol(cos_sims)>1) { # If there are >1 clones in clones2 to compare to
		for(sclone in rownames(cos_sims)) {
			match = c()
			mclones = names(which(cos_sims[sclone,]>0))
			if(length(mclones)>0) {
				for(iclone in mclones) {
					match_samples = rownames(ccfs2)[which(abs(ccfs2[iclone,] - ccfs1[sclone,])<max_dif & ccfs2[iclone,] >= min_ccf_to_call_clone)]
					ssamples = rownames(ccfs1)[which(ccfs1[sclone,]>=min_ccf_to_call_clone)]
					if(is.null(ssamples)) {
						match = c(match, FALSE)
					} else if(length(match_samples)==1) {
						match = c(match, TRUE)
					} else if(length(intersect(ssamples, match_samples))==length(ssamples)) {
						match = c(match, TRUE)
					} else if(length(match_samples)>1) {
						match = c(match, TRUE)
					} else {
						match = c(match, FALSE)
					}
				}
				if(length(which(match))>0) {
					matching_clones = rbind(matching_clones, cbind(sclone,names(which(cos_sims[sclone,]>0))[match]))
				}
			}
		}
	} else if (nrow(cos_sims)==1) {
		matching_clones = rbind(matching_clones, c(rownames(cos_sims), colnames(cos_sims)))
	} else { # If there is only a single clone in clones2 to compare to
		if(length(which(cos_sims[,1]>0))>1) { # If more than more clone in clones1 matches with the only clone in clones2,
						      # Choose the one with the smallest summed difference in CCF
			all_cossim = abs(ccfs1[names(which(cos_sims[,1]>0)),] - ccfs2[colnames(cos_sims),])
			if(class(all_cossim)=='matrix') {
				all_cossim = rowSums(all_cossim)
			} 
			clone1_with_min_cossim = names(which(all_cossim==min(all_cossim)))
			matching_clones = rbind(matching_clones, cbind(clone1_with_min_cossim, colnames(cos_sims)))
		} else {                              # If there is only one clone in clones1 matching with the only clone in clones2,
						      # Take that pair
			matching_clones = rbind(matching_clones, cbind(names(which(cos_sims[,1]>0)), colnames(cos_sims)))
		}
	}
	matching_clones = as.data.frame(matching_clones)
	colnames(matching_clones) = c('clone','matching_clone')
	# Append those that don't match
	rest = setdiff(rownames(ccfs1), matching_clones$clone)
	if(length(rest)>0) {
		rest = cbind(rest,NA)
		colnames(rest) = colnames(matching_clones)
		matching_clones = rbind(matching_clones, rest)
	}
	# If a clone in clones1 matches more than one clone in clones2, pick the one with the smallest summed difference in CCF
	multi_match = as.character(matching_clones$clone[duplicated(matching_clones$clone)])
	if(length(multi_match)>0) {
		for(clone in multi_match) {
			mclones = as.character(matching_clones$matching_clone[matching_clones$clone==clone])
			mccfs = ccfs2[mclones,]
			if('matrix' %in% class(mccfs)) {
				for(p in 1:nrow(mccfs)) {mccfs[p,] = mccfs[p,] - ccfs1[clone,]}
				fmclone = names(which(rowSums(abs(mccfs))==min(rowSums(abs(mccfs)))))
			} else {
				fmclone = names(which(abs(mccfs - ccfs1[clone,])==min(abs(mccfs - ccfs1[clone,]))))
			}
			matching_clones$matching_clone[matching_clones$clone==clone] = fmclone
		}
	}
	matching_clones = unique(matching_clones)
	return(matching_clones)
}

##################################################################################################
## Read all clone information
library(reshape2)
library(heatmap.plus)
library(plyr)

WDIR = '/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/dataset.R')
# Read subclones from SNVs, indels and SVs
subs_clones = read.delim('data/patient_level_subclones_substitutions.txt', header=T, sep="\t", stringsAsFactors=F)
indel_clones = read.delim('data/patient_level_subclones_indels.txt', header=T, sep="\t", stringsAsFactors=F)
sv_clones = read.delim('data/patient_level_subclones_svs.txt', header=T, sep="\t", stringsAsFactors=F)

## Find the substitution clones with matching SV clusters
for(sclone in unique(subs_clones$clone)) {
	presence_summary = paste(sort(unique(subs_clones[subs_clones$ccf>0.05 & subs_clones$clone==sclone,c('sample_short')])), collapse="/")
	subs_clones[subs_clones$clone==sclone,'presence_summary'] = presence_summary
}
subs_clones$matching_sv = 'no'
subs_clones$matching_sv_clone = NA
all_matches = NULL
for(patient in unique(subs_clones$patient)) {
	clones1 = unique(subs_clones[subs_clones$patient==patient,c('clone','descr','presence_summary')])
	clones2 = unique(sv_clones[sv_clones$patient==patient,c('sv_clone_id','presence.summary')])
	for(i in 1:nrow(clones2)) {
		presence_summary = paste(sort(unlist(strsplit(clones2$presence.summary[i], split="/"))), collapse="/")
		clones2$sv_presence_summary[i] = presence_summary
	}
	matches = merge(x=clones1[,c('clone','presence_summary')], y=clones2[,c('sv_clone_id','sv_presence_summary')], by.x='presence_summary', by.y='sv_presence_summary')
	all_matches = rbind(all_matches, matches)
	subs_clones$matching_sv[subs_clones$clone %in% matches$clone] = 'yes'
	for(clone in matches$clone) {
		subs_clones$matching_sv_clone[subs_clones$clone==clone] = length(matches$sv_clone_id[matches$clone==clone])
	}
}
## Find the substitution clones with matching indel clones
all_matches = NULL
for(patient in unique(subs_clones$patient)) {
	matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
	all_matches = rbind(all_matches, matches)
}
all_matches = unique(all_matches)
subs_clones = merge(x=subs_clones, y=all_matches, by.x='clone', by.y='clone', all.x=T)
subs_clones$matching_indel_clone = subs_clones$matching_clone
subs_clones$matching_clone = NULL
subs_clones$matching_indel = NA
subs_clones$matching_indel[which(is.na(subs_clones$matching_indel_clone) & subs_clones$patient!='118706')] = "no"
subs_clones$matching_indel[which(!is.na(subs_clones$matching_indel_clone) & subs_clones$patient!='118706')] = "yes"
subs_clones$matching_indel_clone = as.character(subs_clones$matching_indel_clone)
indel_clones$matching_subs = 'no'
indel_clones[indel_clones$clone %in% subs_clones$matching_indel_clone,'matching_subs'] = 'yes'
#
merged_ccfs = unique(subs_clones[,c('patient','clone','sample','descr','ccf','matching_indel_clone','matching_indel')])
colnames(merged_ccfs)[5] = 'subs_ccf'
merged_ccfs$key = paste(merged_ccfs$sample, merged_ccfs$matching_indel_clone)
indel_clones$key = paste(indel_clones$sample, indel_clones$clone)
merged_ccfs = merge(x=merged_ccfs, y=indel_clones[,c('key','ccf')], by.x='key', by.y='key')
colnames(merged_ccfs)[ncol(merged_ccfs)] = 'indel_ccf'

erronously_matched = c(
	merged_ccfs[which(merged_ccfs$indel_ccf>0.8 & merged_ccfs$subs_ccf<0.5),c('clone')],
	merged_ccfs[which(merged_ccfs$indel_ccf<0.4 & merged_ccfs$subs_ccf>0.5),c('clone')]
)
subs_clones[subs_clones$clone %in% erronously_matched,'matching_indel'] = 'no'
subs_clones[subs_clones$clone %in% erronously_matched,'matching_indel_clone'] = NA

