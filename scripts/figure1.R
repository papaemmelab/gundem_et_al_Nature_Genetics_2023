WDIR = '/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final/'
setwd(WDIR)
source('scripts/util_dataset.R')
library(gridExtra)

######################################################### Figure 1a
# Panels for Figure 1a. Tile plots for age and stage at diagnosis, disease subtype and sequencing data for patients with >=2 tumors
mmdata = mdata[mdata$MS_ANALYSIS=='yes',]
mmdata = mmdata[order(mmdata$PATIENT),]
multi_mycn1 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='MYCN' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_mycn2 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='MYCN' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
multi_atrx1 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='ATRX' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_atrx2 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='ATRX' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
multi_tert1 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='TERT' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_tert2 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='TERT' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
multi_mdm2_cdk41 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='MDM2-CDK4' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_mdm2_cdk42 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='MDM2-CDK4' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
multi_seg_cnv1 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='SEG_CNV' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_seg_cnv2 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='SEG_CNV' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
multi_num_cnv1 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='NUMERIC_CNV' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='WGS')]
multi_num_cnv2 = mmdata$DNA_SAMPLE[which(mmdata$DISEASE_SUBTYPE=='NUMERIC_CNV' & mmdata$no_tumor_samples>1 & mmdata$COHORT=='IMPACT')]
sample_order = c(multi_mycn1, multi_mycn2, multi_atrx1, multi_atrx2, multi_tert1, multi_tert2, multi_mdm2_cdk41, multi_mdm2_cdk42, multi_seg_cnv1, multi_seg_cnv2, multi_num_cnv1, multi_num_cnv2)
fname = 'data/dummy.tsv'

write.table(mdata[match(sample_order, mdata$DNA_SAMPLE), c('DNA_SAMPLE','PATIENT','AGE_GROUP')], file=fname, col.names=F, row.names=F, sep=" ", quote=F)
system(paste('python', 'scripts/figure1a_cohort_summary_age.py', fname))

write.table(mdata[match(sample_order, mdata$DNA_SAMPLE), c('DNA_SAMPLE','PATIENT','DISEASE_SUBTYPE')], file=fname, col.names=F, row.names=F, sep=" ", quote=F)
system(paste('python', 'scripts/figure1a_cohort_summary_mutation.py', fname))

write.table(mdata[match(sample_order, mdata$DNA_SAMPLE), c('DNA_SAMPLE','PATIENT','STAGE')], file=fname, col.names=F, row.names=F, sep=" ", quote=F)
system(paste('python', 'scripts/figure1a_cohort_summary_stage.py', fname))

write.table(mdata[match(sample_order, mdata$DNA_SAMPLE), c('DNA_SAMPLE','PATIENT','SAMPLE_TYPE')], file=fname, col.names=F, row.names=F, sep=" ", quote=F)
system(paste('python', 'scripts/figure1a_cohort_summary_sample_type.py', fname))

write.table(mdata[match(sample_order, mdata$DNA_SAMPLE), c('DNA_SAMPLE','PATIENT','SEQ')], file=fname, col.names=F, row.names=F, sep=" ", quote=F)
system(paste('python', scripts/figure1a_cohort_summary_seq.py, fname))
#########################################################
# Calculate mutation burden
mdata_wgs = mdata[mdata$SS_ANALYSIS=='yes',]
# Patient-level mutation counts. Take mean whenever multiple diagnostic, 2nd-look, relapse and further relapse samples are available
patient_level = ddply(
	mdata_wgs, c('SAMPLE_TYPE','PATIENT','DISEASE_SUBTYPE','STAGE'), summarise,
	NU_SUBS=mean(NU_SUBS), NU_INDELS=mean(NU_INDELS), NU_SV=mean(NU_SV), no.of.subs.clones=mean(no.of.subs.clones),
	rel_tel_length=mean(rel_tel_length),
	A_SBS18=mean(A_SBS18), B_SBS38=mean(B_SBS38), C_SBS40=mean(C_SBS40), D_TMZ=mean(D_TMZ), E_SBS31=mean(E_SBS31), F_SBS35=mean(F_SBS35),
	A_SBS18_n=mean(A_SBS18_n), B_SBS38_n=mean(B_SBS38_n), C_SBS40_n=mean(C_SBS40_n), D_TMZ_n=mean(D_TMZ_n), E_SBS31_n=mean(E_SBS31_n), F_SBS35_n=mean(F_SBS35_n),
	cna_prop_aberrant=mean(cna_prop_aberrant), cna_prop_clonal=mean(cna_prop_clonal), cna_prop_subclonal=mean(cna_prop_subclonal),
	AADx=mean(AADx), ploidy=mean(PLOIDY)
#	sv_mycn=mean(sv_mycn), sv_nonmycn=mean(sv_nonmycn), cn_mycn=mean(cn_mycn),
)
patient_level = patient_level[patient_level$SAMPLE_TYPE!='reresection',]
patient_level$SAMPLE_TYPE = factor(patient_level$SAMPLE_TYPE, levels=c("diagnosis","t-resection","relapse","frelapse"))
sample_type_cols = sample_type_cols[sample_type_cols$sample_type!="reresection",]
sample_type_cols$sample_type = factor(sample_type_cols$sample_type, levels=c("diagnosis","t-resection","relapse","frelapse"))

### Calculate mutation burden
summary(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='diagnosis']) / 3234.83
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05688 0.32973 0.51641 0.67420 0.94279 2.58252 
summary(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='t-resection']) / 3234.83
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1203  0.3286  0.5468  0.6590  0.7896  2.8845 
summary(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='relapse']) / 3234.83
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1122  0.8658  1.3452  1.5700  1.9942  5.8835 
summary(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='frelapse']) / 3234.83
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2068  1.2761  2.2083  2.6626  3.8631  9.0289 
t.test(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='relapse'], patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='frelapse'])
#	Welch Two Sample t-test
#data:  patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "relapse"] and patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "frelapse"]
#t = -2.6305, df = 38.959, p-value = 0.01215
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -6252.1821  -816.5278
#sample estimates:
#mean of x mean of y 
# 5078.739  8613.094 
t.test(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='diagnosis'], patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='relapse'])
#	Welch Two Sample t-test
#data:  patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "diagnosis"] and patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "relapse"]
#t = -4.7956, df = 56.628, p-value = 1.217e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -4107.989 -1687.615
#sample estimates:
#mean of x mean of y 
# 2180.937  5078.739 
t.test(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='diagnosis'], patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='t-resection'])
#	Welch Two Sample t-test
#data:  patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "diagnosis"] and patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "2nd_look"]
#t = 0.15018, df = 85.855, p-value = 0.881
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -600.5406  698.6928
#sample estimates:
#mean of x mean of y 
# 2180.937  2131.861 
t.test(patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='diagnosis'], patient_level$NU_SUBS[patient_level$SAMPLE_TYPE=='frelapse'])
#	Welch Two Sample t-test
#data:  patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "diagnosis"] and patient_level$NU_SUBS[patient_level$SAMPLE_TYPE == "frelapse"]
#t = -5.2111, df = 28.53, p-value = 1.477e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -8958.443 -3905.871
#sample estimates:
#mean of x mean of y 
# 2180.937  8613.094 
######################################################### Figure 1b-c
# Disease at diagnosis: Number of mutations
library(ggpubr)
colnames(mdata_wgs)[colnames(mdata_wgs)=="NU_SUBS"] = "no_substitutions"
colnames(mdata_wgs)[colnames(mdata_wgs)=="NU_INDELS"] = "no_indels"
colnames(mdata_wgs)[colnames(mdata_wgs)=="NU_SV"] = "no_SV"
colnames(mdata_wgs)[colnames(mdata_wgs)=="DISEASE_SUBTYPE"] = "recurrent_lesion"
colnames(mdata_wgs)[colnames(mdata_wgs)=="SAMPLE_TYPE"] = "sample_type"
colnames(mdata_wgs)[colnames(mdata_wgs)=="PURITY"] = "tumor_purity"
colnames(mdata_wgs)[colnames(mdata_wgs)=="rel_tel_length"] = "relative_telomere_length"
colnames(mdata_wgs)[colnames(mdata_wgs)=="AADx"] = "age_at_Dx_months"
colnames(mdata_wgs)[colnames(mdata_wgs)=="cna_prop_aberrant"] = "prop_aberrant_genome"
colnames(mdata_wgs)[colnames(mdata_wgs)=="PLOIDY"] = "tumor_ploidy"
mdata_wgs$tumor_ploidy = as.numeric(mdata_wgs$tumor_ploidy)
mdata_wgs$tmb = mdata_wgs$no_substitutions / 3234.83
mdata_wgs$plt_nrounds = mdata_wgs$PLT_NROUNDS
mdata_wgs$tmz_nrounds = mdata_wgs$TMZ_NROUNDS

symnum_args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
my_comparisons <- list( c("MYCN", "TERT"), c("MYCN","SEG_CNV"), c("MYCN","ATRX"), c("SEG_CNV", "TERT"), c("SEG_CNV", "ATRX"))
plot_dx_tmb <- ggboxplot(
	mdata_wgs[mdata_wgs$sample_type=="diagnosis",], x="recurrent_lesion", y="tmb", color="recurrent_lesion",
	palette=as.character(disease_subtype_cols$color), add="jitter", ylab='TMB', xlab=""
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Diagnosis")

# Disease at diagnosis: Mutational signatures
my_comparisons <- list( c("MYCN","ATRX"), c("MYCN","TERT"), c("MYCN","SEG_CNV"), c("ATRX","TERT"), c("ATRX","SEG_CNV"), c("TERT","SEG_CNV"))
plot_dx_sbs18 <- ggboxplot(
	mdata_wgs[mdata_wgs$sample_type=="diagnosis",], x="recurrent_lesion", y="A_SBS18", color="recurrent_lesion",
	palette=as.character(disease_subtype_cols$color), add="jitter", ylab='Exposure to A_SBS18', xlab="", ylim=c(0,1)
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Diagnosis")
plot_dx_sbs40 <- ggboxplot(
	mdata_wgs[mdata_wgs$sample_type=="diagnosis",], x="recurrent_lesion", y="C_SBS40", color="recurrent_lesion",
	palette=as.character(disease_subtype_cols$color), add="jitter", ylab='Exposure to C_SBS40', xlab="", ylim=c(0,1)
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Diagnosis")
# Disease at progression: Number of mutations
colnames(patient_level)[colnames(patient_level)=="NU_SUBS"] = "no_substitutions"
colnames(patient_level)[colnames(patient_level)=="NU_INDELS"] = "no_indels"
colnames(patient_level)[colnames(patient_level)=="NU_SV"] = "no_SV"
colnames(patient_level)[colnames(patient_level)=="sv_mycn"] = "no_SV_MYCN"
colnames(patient_level)[colnames(patient_level)=="sv_nonmycn"] = "no_SV_nonMYCN"
colnames(patient_level)[colnames(patient_level)=="DISEASE_SUBTYPE"] = "recurrent_lesion"
colnames(patient_level)[colnames(patient_level)=="SAMPLE_TYPE"] = "sample_type"
colnames(patient_level)[colnames(patient_level)=="rel_tel_length"] = "relative_telomere_length"
colnames(patient_level)[colnames(patient_level)=="no.of.subs.clones"] = "no_subclones"
colnames(patient_level)[colnames(patient_level)=="cn_mycn"] = "cn_MYCN"
patient_level$tmb = patient_level$no_substitutions / 3234.83

my_comparisons <- list(c("diagnosis","t-resection"), c("diagnosis", "relapse"), c("diagnosis", "frelapse"),c("relapse","frelapse"))
plot_prog_tmb <- ggboxplot(
	patient_level, x="sample_type", y="tmb", color="sample_type", palette=as.character(sample_type_cols$color), add="jitter", xlab="", ylab="TMB"
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Progression")
plot_prog_tmb <- ggpar(plot_prog_tmb)

# Disease at progression: Mutational signatures
plot_prog_tmz <- ggboxplot(
	patient_level, x="sample_type", y="D_TMZ", color="sample_type", palette=as.character(sample_type_cols$color), add="jitter", xlab="", ylab="Exposure to D_TMZ"
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Progression")
plot_prog_sbs31 <- ggboxplot(
	patient_level, x="sample_type", y="E_SBS31", color="sample_type", palette=as.character(sample_type_cols$color), add="jitter", xlab="", ylab="Exposure to E_SBS31"
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Progression")
plot_prog_sbs35 <- ggboxplot(
	patient_level, x="sample_type", y="F_SBS35", color="sample_type", palette=as.character(sample_type_cols$color), add="jitter", xlab="", ylab="Exposure to F_SBS35"
) + stat_compare_means(comparisons = my_comparisons, symnum.args=symnum_args) + rremove("legend") + rotate_x_text(45) + labs(title="Progression")

#########################################################
# Plot the relationship btw #rounds of chemotherapy received and the exposure to the associated mutational signature
long_exposures = melt(mdata_wgs[mdata_wgs$STAGE==4,c('DNA_SAMPLE','A_SBS18','B_SBS38','C_SBS40','D_TMZ','E_SBS31','F_SBS35','plt_nrounds','tmz_nrounds')], id.vars=c("DNA_SAMPLE","plt_nrounds","tmz_nrounds"))
long_exposures$tmz_nrounds[long_exposures$tmz_nrounds %in% c(3,4,5,7,16)] = "3+"
long_exposures$tmz_nrounds[long_exposures$tmz_nrounds %in% c(1,2)] = "1-2"
long_exposures$plt_nrounds[long_exposures$plt_nrounds %in% c(4,5,6,7,8)] = "4+"
colnames(long_exposures)[4:5] = c('signature','exposure')
plt_exposures = long_exposures[!is.na(long_exposures$plt_nrounds),]
plt_exposures$plt_nrounds = factor(plt_exposures$plt_nrounds, levels=c('0','1','2','3','4+'))
tmz_exposures = long_exposures[!is.na(long_exposures$tmz_nrounds),]
tmz_exposures$tmz_nrounds = factor(tmz_exposures$tmz_nrounds, levels=c('0','1-2','3+'))
# Color scheme for signatures: A_SBS18, B_SBS38, C_SBS40, D_TMZ, E_SBS31 and F_SBS35
signature_palette = c('#1A83B8','#015580','#899DA4','#178520','#AD381F','#801701')

symnum_args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
my_comparisons <- list( c("0", "1"), c("1", "2"), c("1","3"), c("1","4+"), c("2", "3"), c("2", "4+"), c("3", "4+"))
plot_plt_sbs31 <- ggboxplot(
	plt_exposures[plt_exposures$signature=='E_SBS31',], x="plt_nrounds", y="exposure", color="plt_nrounds",
	palette=rep('#AD381F', times=length(unique(plt_exposures$plt_nrounds))), add="jitter", ylab="Exposure to E_SBS31", xlab="#rounds Rx platinum"
)
plot_plt_sbs31 <- plot_plt_sbs31 + stat_compare_means(comparisons=my_comparisons, symnum.args=symnum_args)
plot_plt_sbs31 <- plot_plt_sbs31 + theme(legend.position="none") + labs(title="Rx platinum vs E_SBS31")

plot_plt_sbs35 <- ggboxplot(
	plt_exposures[plt_exposures$signature=='F_SBS35',], x="plt_nrounds", y="exposure", color="plt_nrounds",
	palette=rep('#801701', times=length(unique(long_exposures$plt_nrounds))), add="jitter", ylab='Exposure to F_SBS35', xlab="#rounds Rx platinum"
)
plot_plt_sbs35 <- plot_plt_sbs35 + stat_compare_means(comparisons=my_comparisons, symnum.args=symnum_args)
plot_plt_sbs35 <- plot_plt_sbs35 + theme(legend.position="none") + labs(title="Rx platinum vs F_SBS35")

my_comparisons <- list( c("0", "1-2"), c("1-2","3+"))
plot_tmz <- ggboxplot(
	tmz_exposures[tmz_exposures$signature=='D_TMZ',], x="tmz_nrounds", y="exposure", color="tmz_nrounds",
	palette=rep('#178520', times=length(unique(tmz_exposures$tmz_nrounds))), add="jitter", ylab="Exposure to D_TMZ", xlab="#rounds Rx TMZ"
)
plot_tmz <- plot_tmz + stat_compare_means(comparisons=my_comparisons, symnum.args=symnum_args)
plot_tmz <- plot_tmz + theme(legend.position="none") + labs(title="Rx TMZ vs D_TMZ")
######################################################### Plot increase in TMB in patients with both diagnostic/2nd_look and relapse/frelapse tumors
mdata_wgs$sample_type[mdata_wgs$sample_type %in% c('reresection','t-resection')] = 'diagnosis'
patient_level = ddply(mdata_wgs, c('PATIENT','sample_type'), summarise, tmb=mean(tmb), no_SV=mean(no_SV), no_indel=mean(no_indels))
those_with_dx_and_relapse = intersect(patient_level[patient_level$sample_type=='diagnosis','PATIENT'], patient_level[patient_level$sample_type=='relapse','PATIENT'])
those_with_dx_and_frelapse = intersect(patient_level[patient_level$sample_type=='diagnosis','PATIENT'], patient_level[patient_level$sample_type=='frelapse','PATIENT'])
compare_dx_and_relapse = merge(
	x=patient_level[patient_level$PATIENT %in% those_with_dx_and_relapse & patient_level$sample_type=='diagnosis',c('PATIENT','tmb','no_SV','no_indel')],
	y=patient_level[patient_level$PATIENT %in% those_with_dx_and_relapse & patient_level$sample_type=='relapse',c('PATIENT','tmb','no_SV','no_indel')],
	by.x='PATIENT', by.y='PATIENT'
)
compare_dx_and_frelapse = merge(
	x=patient_level[patient_level$PATIENT %in% those_with_dx_and_relapse & patient_level$sample_type=='diagnosis',c('PATIENT','tmb','no_SV','no_indel')],
	y=patient_level[patient_level$PATIENT %in% those_with_dx_and_relapse & patient_level$sample_type=='frelapse',c('PATIENT','tmb','no_SV','no_indel')],
	by.x='PATIENT', by.y='PATIENT'
)
compare_dx_and_relapse$tmb.increase = compare_dx_and_relapse$tmb.y / compare_dx_and_relapse$tmb.x
compare_dx_and_frelapse$tmb.increase = compare_dx_and_frelapse$tmb.y / compare_dx_and_frelapse$tmb.x
compare_dx_and_relapse$sample = 'relapse'
compare_dx_and_frelapse$sample = 'frelapse'
my_comparisons <- list( c("relapse", "frelapse"))
plot_tmb_change = ggboxplot(
	rbind(compare_dx_and_relapse, compare_dx_and_frelapse), x='sample', y='tmb.increase', color='black', fill='sample', add="jitter",
	palette=sample_type_cols$color[sample_type_cols$sample_type %in% c('relapse','frelapse')], legend='none',
	title='TMB in patients with Dx and relapse tumors', ylab='Fold change in TMB'
)
plot_tmb_change = plot_tmb_change + stat_compare_means(comparisons=my_comparisons, symnum.args=symnum_args)

library(gtable)
g1 <- ggplotGrob(plot_dx_tmb); g2 <- ggplotGrob(plot_dx_sbs18); g3 <- ggplotGrob(plot_dx_sbs40); g4 <- ggplotGrob(plot_prog_tmb); g5 <- ggplotGrob(plot_tmb_change);
g6 <- ggplotGrob(plot_plt_sbs31); g7 <- ggplotGrob(plot_plt_sbs35); g8 <- ggplotGrob(plot_tmz); g9 <- ggplotGrob(plot_prog_sbs31); g10 <- ggplotGrob(plot_prog_sbs35); g11 <- ggplotGrob(plot_prog_tmz);
g_final = rbind(cbind(g1, g2, g3, g4, g5, g5, size="last"), cbind(g6, g7, g8, g9, g10, g11, size="last"), size="first")
grid.newpage()
pdf('raw/figure_1bcdef.pdf', width=3.75*6,  height=3.75*2)
grid.draw(g_final)
dev.off()

