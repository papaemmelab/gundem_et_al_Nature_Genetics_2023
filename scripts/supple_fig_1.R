library(cowplot)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)

WDIR = '~/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setw(WDIR)
source('scripts/dataset.R')
symnum_args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
my_comparisons1 = list(c("MYCN","TERT"), c("MYCN","ATRX"), c("MYCN","SEG_CNV"), c("TERT","ATRX"), c("TERT","SEG_CNV"), c("ATRX","SEG_CNV"))
my_comparisons2 = list(c("diagnosis","t-resection"), c("diagnosis","relapse"), c("diagnosis","frelapse"), c('t-resection','relapse'), c('t-resection','frelapse'), c("relapse","frelapse"))
my_comparisons3 = list( c("no_xrt","xrt_direct"), c("no_xrt","xrt_other"), c("xrt_direct","xrt_other"))
# Read Indel counts ratio
indel_data = read.table('data/indel_counts.txt', header=T, sep="\t", stringsAsFactors=F)
indel_data = indel_data[!is.na(indel_data$SAMPLE_TYPE),]
indel_data$XRT_TO_TUMOR = factor(indel_data$XRT_TO_TUMOR, levels=c("no_xrt","xrt_direct","xrt_other"))
indel_data$SAMPLE_TYPE[indel_data$SAMPLE_TYPE=='2ndlook'] = 't-resection'
indel_data$SAMPLE_TYPE = factor(indel_data$SAMPLE_TYPE, levels=c("diagnosis","t-resection","relapse","frelapse"))
indel_data$DISEASE_SUBTYPE = factor(indel_data$DISEASE_SUBTYPE, levels=c("MDM2-CDK4","SEG_CNV","ATRX","MYCN","TERT","NUMERIC_CNV"))
indel_data$DISEASE_SUBTYPE = factor(indel_data$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
# Check associations with indel burden
glm_nodel <- glm(NU_INDELS_DEL ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=indel_data)
plot_glm_nodel <- plot_model(glm_nodel, vline.color="grey", show.values=TRUE, title="Indels")
glm_noins <- glm(NU_INDELS_INS ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=indel_data)
plot_glm_noins <- plot_model(glm_noins, vline.color="grey", show.values=TRUE, title="Insertions")
glm_nodelmh <- glm(NU_INDELS_MH ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=indel_data)
plot_glm_nodelmh <- plot_model(glm_nodelmh, vline.color="grey", show.values=TRUE, title="MH-mediated deletions")
glm_nodelrep <- glm(NU_INDELS_REP ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=indel_data)
plot_glm_nodelrep <- plot_model(glm_nodelrep, vline.color="grey", show.values=TRUE, title="Repeat-mediated deletions")
# Read SV event count data. 3 WGS tumors don't have SVs.
sv_data = read.delim('data/sv_event_counts.txt', header=T, sep="\t", stringsAsFactors=F)
sv_data = merge(x=sv_data, y=mdata[mdata$SS_ANALYSIS=='yes',c('DNA_SAMPLE','XRT_TO_TUMOR','SAMPLE_TYPE','DISEASE_SUBTYPE')], by.x='sample', by.y='DNA_SAMPLE')
sv_data$SAMPLE_TYPE[sv_data$SAMPLE_TYPE=='2ndlook'] = 't-resection'
# Total SV number
data_plot_total = ddply(sv_data[sv_data$presence==1,], c("sample", "patient", "XRT_TO_TUMOR","SAMPLE_TYPE","DISEASE_SUBTYPE"), summarise, n=length(event_id))
data_plot_total = data_plot_total[!is.na(data_plot_total$XRT_TO_TUMOR),]
data_plot_total$XRT_TO_TUMOR = factor(data_plot_total$XRT_TO_TUMOR, levels=c("no_xrt","xrt_direct","xrt_other"))
data_plot_total$SAMPLE_TYPE = factor(data_plot_total$SAMPLE_TYPE, levels=sample_type_cols$sample_type)
data_plot_total$DISEASE_SUBTYPE = factor(data_plot_total$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
glm_sv_events <- glm(n ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_total)
plot_glm_sv_events <- plot_model(glm_sv_events, vline.color="grey", show.values=TRUE, title='SV events')
# Different types of SVs
data_plot_type = ddply(sv_data[sv_data$presence==1,], c("sample", "patient", "event_type", "XRT_TO_TUMOR","SAMPLE_TYPE","DISEASE_SUBTYPE"), summarise, n=length(event_id))
data_plot_type = droplevels(data_plot_type[!is.na(data_plot_type$XRT_TO_TUMOR),])
data_plot_type$XRT_TO_TUMOR = factor(data_plot_type$XRT_TO_TUMOR, levels=c("no_xrt","xrt_direct","xrt_other"))
data_plot_type$SAMPLE_TYPE = factor(data_plot_type$SAMPLE_TYPE, levels=c("diagnosis","t-resection","relapse","frelapse"))
data_plot_type$DISEASE_SUBTYPE = factor(data_plot_type$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)

data_plot_type_dcast = dcast(data_plot_type, sample + patient + DISEASE_SUBTYPE + SAMPLE_TYPE + XRT_TO_TUMOR ~ event_type, value.var='n')
colnames(data_plot_type_dcast)[8:9] = c('complex_lt_4','complex_gte_4')
data_plot_type_dcast$chromoplexy[is.na(data_plot_type_dcast$chromoplexy)] = 0
data_plot_type_dcast$chromothripsis[is.na(data_plot_type_dcast$chromothripsis)] = 0
data_plot_type_dcast$complex_lt_4[is.na(data_plot_type_dcast$complex_lt_4)] = 0
data_plot_type_dcast$complex_gte_4[is.na(data_plot_type_dcast$complex_gte_4)] = 0
data_plot_type_dcast$DM[is.na(data_plot_type_dcast$DM)] = 0
data_plot_type_dcast$DM_chromothripsis[is.na(data_plot_type_dcast$DM_chromothripsis)] = 0
data_plot_type_dcast$DM_simple[is.na(data_plot_type_dcast$DM_simple)] = 0
data_plot_type_dcast$DM_tandem_duplication[is.na(data_plot_type_dcast$DM_tandem_duplication)] = 0
data_plot_type_dcast$DM_translocation[is.na(data_plot_type_dcast$DM_translocation)] = 0
data_plot_type_dcast$inversion[is.na(data_plot_type_dcast$inversion)] = 0
data_plot_type_dcast$reciprocal_inversion[is.na(data_plot_type_dcast$reciprocal_inversion)] = 0
data_plot_type_dcast$reciprocal_translocation[is.na(data_plot_type_dcast$reciprocal_translocation)] = 0
data_plot_type_dcast$simple_deletion[is.na(data_plot_type_dcast$simple_deletion)] = 0
data_plot_type_dcast$tandem_duplication[is.na(data_plot_type_dcast$tandem_duplication)] = 0
data_plot_type_dcast$unbalanced_translocation[is.na(data_plot_type_dcast$unbalanced_translocation)] = 0

glm_sv_simple_del <- glm(simple_deletion ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_simple_del <- plot_model(glm_sv_simple_del, vline.color="grey", show.values=TRUE, title="simple deletion")

glm_sv_tandem_dup <- glm(tandem_duplication ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_tandem_dup <- plot_model(glm_sv_tandem_dup, vline.color="grey", show.values=TRUE, title="tandem duplication")

glm_sv_recip_trans <- glm(reciprocal_translocation ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_recip_trans <- plot_model(glm_sv_recip_trans, vline.color="grey", show.values=TRUE, title='reciprocal translocation')

glm_sv_unbal_trans <- glm(unbalanced_translocation ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_unbal_trans <- plot_model(glm_sv_unbal_trans, vline.color="grey", show.values=TRUE, title='unbalanced_translocation')

glm_sv_recip_inv <- glm(reciprocal_inversion ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_recip_inv <- plot_model(glm_sv_recip_inv, vline.color="grey", show.values=TRUE, title='Reciprocal translocation')

glm_sv_complex_lt_4 <- glm(complex_lt_4 ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_complex_lt_4 <- plot_model(glm_sv_complex_lt_4, vline.color="grey", show.values=TRUE, title='Complex n<4')

glm_sv_complex_gte_4 <- glm(complex_gte_4 ~ XRT_TO_TUMOR + SAMPLE_TYPE + DISEASE_SUBTYPE, data=data_plot_type_dcast)
plot_glm_sv_complex_gte_4 <- plot_model(glm_sv_complex_gte_4, vline.color="grey", show.values=TRUE, tile="Complex >=4")

pdf('raw/supple_fig_1.pdf', width=21, height=12.6)
cowplot::plot_grid(
	plot_glm_nodel, plot_glm_noins, plot_glm_nodelmh, plot_glm_nodelrep,
	plot_glm_sv_events, plot_glm_sv_simple_del, plot_glm_sv_tandem_dup, plot_glm_sv_recip_trans,
	plot_glm_sv_unbal_trans, plot_glm_sv_recip_inv, plot_glm_sv_complex_lt_4, plot_glm_sv_complex_gte_4,  nrow=3, align='h'
)
dev.off()

