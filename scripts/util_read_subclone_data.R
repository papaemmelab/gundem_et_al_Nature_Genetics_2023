library(ggpubr)

cosine.similarity <- function(A,B){
	return(sum(A*B)/sqrt(sum(A^2)*sum(B^2)))
}

get_purity_color_palette <- function(x,n=10){
	colorRampPalette(c("orange", "brown"))(n)[cut(x,n)]
}

## Correct errors in BBG output
correct_clonality <- function(bbg, sf_cutoff=0.1) {
	# Do not call anything with subclonal fraction less than 10%
	sel = bbg[,'frac1_A'] <= sf_cutoff
	bbg[sel,'nMaj1_A'] = bbg[sel,'nMaj2_A']
	bbg[sel,'nMin1_A'] = bbg[sel,'nMin2_A']
	bbg[sel,'frac1_A'] = 1
	bbg[sel,c('nMaj2_A','nMin2_A','frac2_A')] = NA

	sel = !is.na(bbg[,'frac2_A']) & bbg[,'frac2_A'] < sf_cutoff
	bbg[sel,'frac1_A'] = 1
	bbg[sel,c('nMaj2_A','nMin2_A','frac2_A')] = NA

	return(bbg)
}

get_bbg_summary <- function(bbg_file, wga=FALSE){
	bbg = read.delim(bbg_file, header=T, sep="\t", stringsAsFactors=F)
	bbg = correct_clonality(bbg, sf_cutoff=0.11)
	bbg$id = paste(bbg$sample, bbg$seg_order, sep="__")
	# Count the number of segments off-trunk
	cna_summary = NULL
	if(length(unique(bbg$clone))>1) {
		cna_summary = merge(
			x=as.data.frame(table(bbg[!grepl('_T0_',bbg$sample) & !grepl('-T0-',bbg$sample) & !grepl('trunk',bbg$clone),c('clone')])),
			y=as.data.frame(table(unique(bbg[!grepl('_T0_',bbg$sample) & !grepl('-T0-',bbg$sample) & !grepl('trunk',bbg$clone),c('sample','clone')])$clone)),
			by.x='Var1', by.y='Var1'
		)
		colnames(cna_summary) = c('subclone','total_segs','no_samples')
		cna_summary$total_segs = as.integer(cna_summary$total_segs / cna_summary$no_samples)
		cna_summary$subclone = as.character(cna_summary$subclone)
	}
	cna_summary = rbind(cna_summary, c('trunk', length(bbg[grepl('_T0_',bbg$sample) | grepl('-T0-',bbg$sample) & bbg$color!="darkgrey",'sample']), length(unique(bbg$sample))))
	colnames(cna_summary) = c('subclone','total_segs','no_samples')
	cna_summary = as.data.frame(cna_summary)
	cna_summary$total_segs = as.integer(as.character(cna_summary$total_segs))
	cna_summary$prop_segs = cna_summary$total_segs / sum(cna_summary$total_segs)
	return(cna_summary)
}

read_sv_clones <- function(sv_file) {
        svs = read.table(sv_file, header=T, sep="\t", stringsAsFactors=F, comment.char="")
	if("evidence" %in% colnames(svs)) {svs$confidence = svs$evidence}
        svs = svs[svs$confidence %in% c("high","high_conf"),]
	svs$sv_clone_id = paste(svs$presence.summary, svs$subs_clone, sep="__")
	svs$rearr_type = "translocation"
	svs$rearr_type[svs$chr1==svs$chr2 & paste(svs$strand1, svs$strand2)=="- -"] = "tandem_duplication"
	svs$rearr_type[svs$chr1==svs$chr2 & paste(svs$strand1, svs$strand2)=="+ +"] = "deletion"
	svs$rearr_type[svs$chr1==svs$chr2 & paste(svs$strand1, svs$strand2)=="+ -"] = "inversion"
	svs$rearr_type[svs$chr1==svs$chr2 & paste(svs$strand1, svs$strand2)=="- +"] = "inversion"
	svs$csv_locus = tolower(svs$csv_locus)
	for(id in unique(svs$event_id)) {
		if(length(unique(svs[svs$event_id==id,'csv_locus']))>1) {
			svs[svs$event_id==id,'csv_locus'] = paste(unique(svs[svs$event_id==id,'csv_locus']), collapse="|")
		}
	}
	sv_events = unique(svs[,c('sv_clone_id','presence.summary','subs_clone','event_id','event_type','csv_locus')])
	event_counts = ddply(sv_events, c('sv_clone_id','presence.summary','subs_clone','event_type','csv_locus'), summarise, N=length(event_id))
	event_counts$driver_locus = event_counts$csv_locus
	event_counts$csv_locus = NULL
	return(event_counts)
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
	if(class(cos_sims)!='matrix') {cos_sims = as.matrix(cos_sims)}
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
			if(class(mccfs)=='matrix') {
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
source('/scripts/util_dataset.R'
clone_files = read.table('data/patient_level_clone_files.txt', header=T, sep="\t", stringsAsFactors=F)
subs_clones = NULL
indel_clones = NULL
cn_clones = NULL
sv_clones = NULL
sel_columns1 = c("cluster.no","descr","descr2","no.muts.in.cluster","prop.muts.in.cluster","colour","order")
sel_columns2 = c("cluster.no","descr","no.muts.in.cluster","prop.muts.in.cluster","colour","subs_clone")
sel_columns3 = c("clone","subs_clone","color","prop_segs","total_segs")
for(i in 1:nrow(clone_files)) {
	## Read substitution clones
	cur_clone = read.table(clone_files$subs_file[i], header=T, sep="\t", comment.char="", stringsAsFactors=F)
	if(!('no.muts.in.cluster' %in% colnames(cur_clone))) {cur_clone$no.muts.in.cluster = cur_clone$no.of.mutations}
	cur_clone$prop.muts.in.cluster = cur_clone$no.muts.in.cluster / sum(cur_clone$no.muts.in.cluster)
	subs_clones = rbind(subs_clones, cbind(clone_files$patient[i], melt(cur_clone[,c(sel_columns1, grep('T',colnames(cur_clone),value=T))], id.vars=sel_columns1)))
	## Read indel clones
	cur_clone = read.table(clone_files$indel_file[i], header=T, sep="\t", comment.char="", stringsAsFactors=F)
	if(!('no.muts.in.cluster' %in% colnames(cur_clone))) {cur_clone$no.muts.in.cluster = cur_clone$no.of.mutations}
	cur_clone$prop.muts.in.cluster = cur_clone$no.muts.in.cluster / sum(cur_clone$no.muts.in.cluster)
	indel_clones = rbind(indel_clones, cbind(clone_files$patient[i], melt(cur_clone[,c(sel_columns2, grep('T',colnames(cur_clone),value=T))], id.vars=sel_columns2)))
	## Read CN clones
	if(!is.na(clone_files$cn_file[i])) {
		cur_clone = read.table(clone_files$cn_file[i], header=T, sep="\t", comment.char="", stringsAsFactors=F)
		cur_clone$clone = gsub('-','_',gsub(' ','_',gsub(' / ','__',gsub('%','',cur_clone$clone))))
		bbg_file = sub('cn_clones.txt','all_subclones.txt',clone_files$cn_file[i])
		bbg_summary = get_bbg_summary(bbg_file, wga=FALSE)
		bbg_summary$subclone = gsub('-','_',gsub(' ','_',gsub(' / ','__',gsub('%','',bbg_summary$subclone))))
		cur_clone = merge(x=cur_clone, y=bbg_summary, by.x='clone', by.y='subclone', all.x=T)
		cn_clones = rbind(cn_clones, cbind(clone_files$patient[i], melt(cur_clone[,c(sel_columns3, grep('T',colnames(cur_clone),value=T))], id.vars=sel_columns3)))
	}
	## SV clones
	if(!is.na(clone_files$sv_file[i])) {
		cur_clone = read_sv_clones(clone_files$sv_file[i])
		cur_clone = cbind(clone_files$patient[i], cur_clone)
		colnames(cur_clone)[1] = 'patient'
		sv_clones = rbind(sv_clones, cur_clone)
	}
}
subs_clones = as.data.frame(subs_clones)
colnames(subs_clones)[1] = 'patient'
subs_clones$patient = sub('I-H-','H',subs_clones$patient)
subs_clones$clone = apply(subs_clones[,c('patient','descr')], 1, paste, collapse="__")
subs_clones$sample = subs_clones$variable
subs_clones$ccf = subs_clones$value
subs_clones$variable = NULL
subs_clones$value = NULL
subs_clones = merge(x=subs_clones, y=mdata[,c('DNA_SAMPLE','PURITY','no_tumor_samples','AADx')], by.x='sample', by.y='DNA_SAMPLE', all.x=T)
subs_clones$purity_color = get_purity_color_palette(subs_clones$PURITY, 100)
subs_clones = merge(x=subs_clones, y=ddply(unique(subs_clones[subs_clones$ccf>0.01,c('clone','sample')]), c("clone"), summarise, n=length(sample)), by.x='clone', by.y='clone')
colnames(subs_clones)[ncol(subs_clones)] = "n_samples"
indel_clones = as.data.frame(indel_clones)
colnames(indel_clones)[1] = 'patient'
indel_clones$patient = sub('I-H-','H',indel_clones$patient)
indel_clones$clone = apply(indel_clones[,c('patient','descr')], 1, paste, collapse="__")
indel_clones$sample = indel_clones$variable
indel_clones$ccf = indel_clones$value
indel_clones$variable = NULL
indel_clones$value = NULL
indel_clones$subs_clone[!is.na(indel_clones$subs_clone)] = apply(indel_clones[!is.na(indel_clones$subs_clone),c('patient','subs_clone')], 1, paste, collapse="__")
cn_clones = as.data.frame(cn_clones)
colnames(cn_clones)[1] = 'patient'
colnames(cn_clones)[2] = 'descr'
cn_clones$patient = sub('I-H-','H',cn_clones$patient)
cn_clones$sample = cn_clones$variable
cn_clones$ccf = cn_clones$value
cn_clones$variable = NULL
cn_clones$value = NULL
cn_clones$clone = apply(cn_clones[,c('patient','descr')], 1, paste, collapse="__")
cn_clones$subs_clone[!is.na(cn_clones$subs_clone)] = apply(cn_clones[!is.na(cn_clones$subs_clone),c('patient','subs_clone')], 1, paste, collapse="__")
sv_clones = as.data.frame(sv_clones)
colnames(sv_clones)[1] = 'patient'
sv_clones$patient = sub('I-H-', 'H', sv_clones$patient)
sv_clones$N = as.numeric(sv_clones$N)
## Find the substitution clones with matchin CN clones
sel_patients = setdiff(unique(subs_clones$patient), c(
	'H112910','H116990','H118706','H118733','H131236','H131239','H132380','H132383','H132384','H132387','H132388','H134819','H135089'))
all_matches = NULL
for(patient in sel_patients) {
	matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
	all_matches = rbind(all_matches, matches)
}
all_matches$matching_clone[all_matches$clone %in% c(
	'H103207__T6_22','H103207__T10_100__T5_4__T6_6','H116989__T2_56__T3_100','H132374__T4_6__T5_T6_100','H134722__T21_23__T22_5__T23_6','H134722__T22_8__T23_93',
	'H134821__T1_4__T17_T18_100','H134821__T17_15__T18_100','H134822__T10_100__T7_12','H135421__s65','H135467__s63','H158184__s20','H131238__s68'
)] = NA
all_matches = all_matches[!(all_matches$clone=='134817__T12_2__T13_T14_100' & all_matches$matching_clone=='134817__T14_65__T13_100'),]
patient = "H118706"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.08);
matches$matching_clone = as.character(matches$matching_clone)
matches[matches$clone %in% c('118706__T20_mets_100'),'matching_clone'] = NA
matches[matches$clone=='118706__mets_100','matching_clone'] = 'mets_100'
matches[matches$clone=='118706__T191_100','matching_clone'] = '118706__T191_100'
matches[matches$clone=='118706__T192_89','matching_clone'] = '118706__T192_67'
matches[matches$clone=='118706__T193_83','matching_clone'] = '118706__T193_100'
all_matches = rbind(all_matches, matches)
patient = "H118733"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05);
matches$matching_clone = as.character(matches$matching_clone)
matches[matches$clone=='118733__T16_47__T17_T18_T19_T20_100','matching_clone'] = '118733__T17_55__T18_65__T19_34__T20_70'
all_matches = rbind(all_matches, matches)
patient = "H132383"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05);
matches$matching_clone = as.character(matches$matching_clone)
matches[matches$clone=='132383__T10_6__T12_80__T7_84','matching_clone'] = '132383__T12_77__T7_100'
matches[matches$clone=='H132383__T10_88__T11_59__T8_100','matching_clone'] = NA
all_matches = rbind(all_matches, matches)
patient = "H132384"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.1);
all_matches = rbind(all_matches, matches)
patient = "H132387"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.09)
matches[matches$clone=='H132387__T5_53__T6_T7_T1_100','matching_clone'] = NA
all_matches = rbind(all_matches, matches)
patient = "H132388"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
matches$matching_clone = as.character(matches$matching_clone)
all_matches = rbind(all_matches, matches)
patient = "H134819"; matches = find_matching_clones(clones1=subs_clones, clones2=cn_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.07)
all_matches = rbind(all_matches, matches)
all_matches = unique(all_matches)
subs_clones = merge(x=subs_clones, y=all_matches, by.x='clone', by.y='clone', all.x=T)
subs_clones$matching_cn_clone = subs_clones$matching_clone
subs_clones$matching_clone = NULL
subs_clones$matching_cn = NA
subs_clones$matching_cn[which(is.na(subs_clones$matching_cn_clone) & !(subs_clones$patient %in% c('112910','116990','132380','135089')))] = "no"
subs_clones$matching_cn[which(!is.na(subs_clones$matching_cn_clone) & !(subs_clones$patient %in% c('112910','116990','132380','135089')))] = "yes"
## Find the substitution clones with matching indel clones
all_matches = NULL
for(patient in setdiff(unique(subs_clones$patient),c('H116989','H116990','H131213','H118733','H132384','H133120','H134722','H134818','H134819'))) {
	matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
	all_matches = rbind(all_matches, matches)
}
patient = "H116989"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.20, min_cos_sim=0.7)
all_matches = rbind(all_matches, matches)
patient = "H116990"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9)
all_matches = rbind(all_matches, matches)
patient = "H118733"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.04)
all_matches = rbind(all_matches, matches)
patient = "H131213"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.04)
all_matches = rbind(all_matches, matches)
patient = "H132384"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
all_matches = rbind(all_matches, matches)
patient = "H133120"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.07)
all_matches = rbind(all_matches, matches)
patient = "H134722"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.15, min_cos_sim=0.9, min_ccf_to_call_clone=0.07)
all_matches = rbind(all_matches, matches)
patient = "H134818"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.20, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
all_matches = rbind(all_matches, matches)
patient = "H134819"; matches = find_matching_clones(clones1=subs_clones, clones2=indel_clones, patient=patient, max_dif=0.05, min_cos_sim=0.9, min_ccf_to_call_clone=0.05)
all_matches = rbind(all_matches, matches)
all_matches$matching_clone[all_matches$clone %in% c(
	'H103207__T11_15','H116984__T1_100_T2_27','H131208__s51','H131236__s64','H132372__T3_50__T4_100','H132383__T10_88__T11_59__T8_100',
	'H132387__T5_37__T1_8',
	'H132388__T8_36__T1_23','H132388__T7_100__T1_61__T9_74','H132397__T3_65__T4_100','H134821__T14_24','H135421__s28','H118733__T20_100','H131213__s55'
)] = NA
all_matches$matching_clone[all_matches$clone=='H131236__trunk'] = 'H131236__trunk'
all_matches = unique(all_matches)
subs_clones = merge(x=subs_clones, y=all_matches, by.x='clone', by.y='clone', all.x=T)
subs_clones$matching_indel_clone = subs_clones$matching_clone
subs_clones$matching_clone = NULL
subs_clones$matching_indel = NA
subs_clones$matching_indel[which(is.na(subs_clones$matching_indel_clone) & subs_clones$patient!='118706')] = "no"
subs_clones$matching_indel[which(!is.na(subs_clones$matching_indel_clone) & subs_clones$patient!='118706')] = "yes"
subs_clones$matching_indel_clone = as.character(subs_clones$matching_indel_clone)
subs_clones$matching_cn_clone = as.character(subs_clones$matching_cn_clone)
## Read mutational signature results for substitution clones
patient_exposures = read.table('data/patients_signature_exposures.txt', header=T, sep="\t")
patient_exposures$clone_id = patient_exposures$clone
sign_columns = c('A_SBS18','B_SBS38','C_SBS40','D_TMZ','E_SBS31','F_SBS35')
patient_exposures[,sign_columns] = patient_exposures[,sign_columns] / rowSums(patient_exposures[,sign_columns])
patient_exposures$patient = NULL
patient_exposures$clone = NULL
subs_clones = merge(x=subs_clones, y=patient_exposures, by.x='clone', by.y='clone_id', all.x=T)
## Merge mutational signature with substitution clone info and find the substitution cluster exposed to Plt and TMZ
sel_clones = which(!is.na(subs_clones$E_SBS31) & subs_clones$descr=='trunk' & subs_clones$no.muts.in.cluster>200 & !(subs_clones$patient %in% c('H118733','H132379','H116988')))
pt_cutoff = max(unique(subs_clones[sel_clones,c('patient','E_SBS31')])$E_SBS31)
tmz_cutoff = max(subs_clones[sel_clones,c('D_TMZ')])
paste('Pt exposure cutoff = ', pt_cutoff, sep="")
paste('TMZ exposure cutoff = ', tmz_cutoff, sep="")
subs_clones$plt_exposed = rep(NA, nrow(subs_clones))
subs_clones$tmz_exposed = rep(NA, nrow(subs_clones))
subs_clones[which(subs_clones$E_SBS31>pt_cutoff & subs_clones$no.muts.in.cluster>=200),'plt_exposed'] = "yes"
subs_clones[which(subs_clones$E_SBS31<=pt_cutoff & subs_clones$no.muts.in.cluster>=200),'plt_exposed'] = "no"
subs_clones[which(subs_clones$no.muts.in.cluster<200),'plt_exposed'] = NA
subs_clones[which(subs_clones$D_TMZ>tmz_cutoff & subs_clones$no.muts.in.cluster>=200),'tmz_exposed'] = "yes"
subs_clones[which(subs_clones$D_TMZ<=tmz_cutoff & subs_clones$no.muts.in.cluster>=200),'tmz_exposed'] = "no"
subs_clones[subs_clones$no.muts.in.cluster<200,'tmz_exposed'] = NA
## For each clone, determine the number of mutations attributed to therapy and
## Calculate the proportion of mutations not attributed to therapy compared all the mutations in the patient
subs_clones$no.muts.in.cluster.plt = subs_clones[,'no.muts.in.cluster'] * subs_clones[,'E_SBS31']
subs_clones$no.muts.in.cluster.tmz = subs_clones[,'no.muts.in.cluster'] * subs_clones[,'D_TMZ']
subs_clones$no.muts.in.cluster.notherapy = subs_clones$no.muts.in.cluster - (subs_clones$no.muts.in.cluster.plt + subs_clones$no.muts.in.cluster.tmz)
subs_clones$prop.muts.in.cluster.notherapy = -1
for(patient in unique(subs_clones$patient)) {
	mut_sum = sum(unique(subs_clones[subs_clones$patient==patient,c('clone','no.muts.in.cluster.notherapy')])$no.muts.in.cluster.notherapy)
	mut_prop = subs_clones[subs_clones$patient==patient,'no.muts.in.cluster.notherapy'] / mut_sum
	subs_clones[subs_clones$patient==patient,'prop.muts.in.cluster.notherapy'] = mut_prop
}
############################################################################# Compare clones
