# Read all subclone data
WDIR = '~/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final'
setwd(WDIR)
source('scripts/util_read_subclone_data.R')
no_wgs_tumors = as.data.frame(table(mdata$PATIENT[mdata$SEQ=='WGS']))
colnames(no_wgs_tumors) = c('patient','no_wgs_tumor_samples')
subs_clones = merge(x=subs_clones, y=no_wgs_tumors, by.x='patient', by.y='patient')
indel_clones = merge(x=indel_clones, y=no_wgs_tumors, by.x='patient', by.y='patient')

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
prop_subs_on_trunk = unique(subs_clones[subs_clones$descr=='trunk' & subs_clones$no_wgs_tumor_samples>1,c('patient','prop.muts.in.cluster','no.muts.in.cluster','no_wgs_tumor_samples')])
colnames(prop_subs_on_trunk)[2:4] = c('prop_subs_on_trunk','no_subs_on_trunk','no_sample')
prop_subs_on_trunk = merge(x=prop_subs_on_trunk, y=unique(mdata[,c('PATIENT','DISEASE_SUBTYPE','AADx','AGE_GROUP')]), by.x='patient', by.y='PATIENT')
prop_indels_on_trunk = unique(indel_clones[indel_clones$descr=='trunk' & indel_clones$no_wgs_tumor_samples>1,c('patient','prop.muts.in.cluster')])
colnames(prop_indels_on_trunk)[2] = 'prop_indels_on_trunk'
prop_muts_on_trunk = merge(x=prop_subs_on_trunk, y=prop_indels_on_trunk, by.x='patient', by.y='patient')
prop_muts_on_trunk = merge(x=prop_muts_on_trunk, y=as.data.frame(table(unique(subs_clones[,c('patient','clone')])$patient)), by.x='patient', by.y='Var1')
colnames(prop_muts_on_trunk)[ncol(prop_muts_on_trunk)] = 'no_subclones'
prop_muts_on_trunk$no_subclones_per_sample = prop_muts_on_trunk$no_subclones / prop_muts_on_trunk$no_sample
prop_muts_on_trunk$no_sample[prop_muts_on_trunk$no_sample>=4] = "4+"
prop_muts_on_trunk$no_sample = factor(prop_muts_on_trunk$no_sample, levels=c('2','3','4+'))

# Check the correlates of trunk length
prop_muts_on_trunk_evolution = merge(x=prop_muts_on_trunk, y=evolution[,c('Patient','X.Pri','stage')], by.x='patient', by.y='Patient')
prop_muts_on_trunk_subset = prop_muts_on_trunk_evolution[prop_muts_on_trunk_evolution$X.Pri>0,]
tertile_limits = quantile(prop_muts_on_trunk_subset$no_subs_on_trunk, c(0:3/3))
prop_muts_on_trunk_subset$tertile = with(prop_muts_on_trunk_subset, cut(no_subs_on_trunk, tertile_limits, include.lowest=T, labels=c("lower","middle","higher")))
glm1 <- glm(no_subs_on_trunk ~ DISEASE_SUBTYPE + AADx + no_sample + stage, data=prop_muts_on_trunk_subset)

library(sjPlot)
library(sjlabelled)
library(sjmisc)
plot_model(glm1, vline.color="grey", show.values=TRUE)
# The only thing that is significantly correlated with the number of SNVs on trunk is the age at diagnosis (p=0.00209)

# Check the prevalence of disease subtype in different trunk length tertiles
dtype_counts = table(prop_muts_on_trunk_subset$DISEASE_SUBTYPE, prop_muts_on_trunk_subset$tertile)
for(column in colnames(dtype_counts)) {dtype_counts[,column] = dtype_counts[,column] / colSums(dtype_counts)[column]}
dtype_counts = as.data.frame(dtype_counts)
colnames(dtype_counts) = c('disease_subtype','trunk_length_tertile','freq')
dtype_counts$disease_subtype = factor(dtype_counts$disease_subtype, levels=disease_subtype_cols$disease_subtype)
barplot_dtype_trunk_length <- ggbarplot(dtype_counts, x='trunk_length_tertile', y='freq', fill='disease_subtype', palette=disease_subtype_cols$color, legend='none')
# Check the prevalence of low vs high stage in different trunk length tertiles
stage_counts = table(prop_muts_on_trunk_subset$stage, prop_muts_on_trunk_subset$tertile)
for(column in colnames(stage_counts)) {stage_counts[,column] = stage_counts[,column] / colSums(stage_counts)[column]}
stage_counts = as.data.frame(stage_counts)
colnames(stage_counts) = c('stage','trunk_length_tertile','freq')
barplot_stage_trunk_length <- ggbarplot(stage_counts, x='trunk_length_tertile', y='freq', fill='stage', palette=c('#E1E5E5','#2F4F4F'), legend='none')
scatterplot_age_trunk_length <- ggscatter(prop_muts_on_trunk_subset, x='AADx', y='no_subs_on_trunk', color='DISEASE_SUBTYPE', palette=disease_subtype_cols$color, add='reg.line', conf.int=TRUE, add.params=list(color="blue",fill="lightgray"), legend='none') + stat_cor(method="spearman", label.x=100, label.y=5000)
pdf(paste('raw/figure_2a.pdf', sep=''), width=9, height=3)
cowplot::plot_grid(scatterplot_age_trunk_length, barplot_dtype_trunk_length, barplot_stage_trunk_length, nrow=1, align='h', axis='b')
dev.off()

################################################### MRCA Analysis
library(gthor)
library(cowplot)
library(lme4)
library(pander)

### Helper functions for MRCA analysis
# Function to estimate the time of MRCA
estimateMRCA <- function(input.df, lmer.model){
	# Bootstrapping function for MRCA
	mySumm <- function(.) {
		predict(., newdata=mrca.df, re.form=NULL)
	}
	# Collapse bootstrap into median, 95% prediction interval (PI)
	sumBoot <- function(merBoot) {
		return(
			data.frame(
				fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
				lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
				upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
			)
		)
	}
	# setup variable
	mrca.df <- input.df
	mrca.df <- mrca.df[!duplicated(mrca.df$patient),]
	# predict the age of MCRA
	mrca.df$MRCA.pred.age <- predict(lmer.model, newdata = mrca.df)
	# Now generate 95% CIs on the predictions for MRCA timing
	# Note that CIs for lme models are challenging - bootstrapping appears to be most robust
	## bootMer = Perform model-based (Semi-)parametric bootstrap for mixed models.
	boot1 <- bootMer(lmer.model, mySumm, nsim=1000, use.u=FALSE, type="parametric")
	PI.boot1 <- sumBoot(boot1)
	PI.boot1 <- PI.boot1[complete.cases(PI.boot1),]
	mrca.df2 <- mrca.df[complete.cases(mrca.df$MRCA.pred.age),]
	PI.boot1$AADx <- mrca.df2$AADx
	PI.boot1$pred.fit <- mrca.df2$MRCA.pred.age
	PI.boot1$Time.lag <- PI.boot1$AADx - PI.boot1$fit
	PI.boot1$Time.lag[PI.boot1$Time.lag < -9] <- -9 # Time lag < -9 are set to -9. Earlier than conception does not make sense.
	PI.boot1$lag.lwr <- PI.boot1$AADx - PI.boot1$upr
	PI.boot1$lag.lwr[PI.boot1$lag.lwr < -9] <- -9 # Time lag < -9 are set to -9. Earlier than conception does not make sense.
	PI.boot1$lag.upr <- PI.boot1$AADx - PI.boot1$lwr
	PI.boot1$lag.upr[PI.boot1$lag.upr < -9] <- -9 # Time lag < -9 are set to -9. Earlier than conception does not make sense.
	PI.boot1$ID <- mrca.df2$patient
	PI.boot1 <- PI.boot1[order(PI.boot1$AADx),]
	# order by time lag
	PI.boot1 <- PI.boot1[order(PI.boot1$Time.lag),]
	return(PI.boot1)
}
# plotting function for MRCA
plotMRCA <- function(input.df){
	PI.boot1 <- input.df
	PI.boot1$pch = 16
	PI.boot1$pch[(PI.boot1$AADx - PI.boot1$upr)>0] = 17
	custom_axis_range = c(-10, (0:31) * 10)
	par(mfrow=c(1, 1), xpd=T, mar=c(7, 12, 5, 5))
	plot(
		PI.boot1$Time.lag, 1:nrow(PI.boot1), pch = PI.boot1$pch, xlim = c(min(PI.boot1$lag.lwr), max(PI.boot1$lag.upr)), axes = FALSE,
		xlab = "Estimated time between MRCA and first sample collection (months)", ylab = "", col = PI.boot1$color, cex = 2, cex.lab = 2
	)
	segments(x0=PI.boot1$lag.lwr, x1=PI.boot1$lag.upr, y0=1:nrow(PI.boot1), col = PI.boot1$color)
	par(new = T)
	plot(
		PI.boot1$Time.lag, 1:nrow(PI.boot1), pch = PI.boot1$pch, xlim = c(min(PI.boot1$lag.lwr), max(PI.boot1$lag.upr)),
		axes = FALSE, xlab = "", ylab = "", col = PI.boot1$color, cex = 1.7
	)
	# bottom axis
	axis(side = 1, at = custom_axis_range, labels = custom_axis_range, cex.axis = 2)
	# top axis
	axis(side = 3, at = custom_axis_range, labels = custom_axis_range, cex.axis = 2)
	grid(ny = c(-10, (0:31) * 10), nx = NULL)
	text(x = -11, y = 1:nrow(PI.boot1), labels = PI.boot1$label, adj = 1, cex = 2)
	# Shaded aread for in-utero period
	polygon(x=c(0, -10, -10, 0), y=c(0, 0, nrow(PI.boot1), nrow(PI.boot1)), col=adjustcolor("grey",alpha.f=0.2), border=NA)
	# Diagnosis
	points(PI.boot1$AADx, 1:nrow(PI.boot1), pch = 4, cex = 1.7, lwd = 2)
}
############################################ Read the subclones and mutational signatures
subclones = unique(subs_clones[,c('patient','descr','cluster.no','no.muts.in.cluster','prop.muts.in.cluster')])
subclones_sigs = NULL
for(patient in subclones$patient) {
	sig_file = paste('data/patient_level_clone_files/', patient, '/signature_exposures.txt', sep='')
	if(length(sig_file)==0) {next}
	clones_file = paste('data/patient_level_clone_files/', patient, '/subs_clones.txt', sep='')
	clone_data = read.delim(clones_file, header=T, sep="\t", stringsAsFactors=F)
	sigs_data = read.table(sig_file, header=T, sep="\t", stringsAsFactors=F)
	sigs_data = merge(y=clone_data[,c('cluster.no','label')], x=sigs_data, by.y='label', by.x='cluster')
	sigs_data = merge(x=sigs_data, y=subclones[subclones$patient==patient,], by.x='cluster.no', by.y='cluster.no')
	subclones_sigs = rbind(subclones_sigs, sigs_data)
}
colnames(subclones_sigs)[1:2] = c('dp_cluster_no','cluster_order')
subclones_sigs = merge(x=subclones_sigs, y=unique(mdata[,c('PATIENT','AADx','DISEASE_SUBTYPE')]), by.x='patient', by.y='PATIENT')
selected_pts = subclones_sigs$patient %in% unique(mdata$PATIENT[mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','t-resection')]) & subclones_sigs$C_SBS40>0
subclones_sigs_sub = subclones_sigs[selected_pts,]
############################################ Check correlation between mutational signatures and age at diagnosis
plot_lm_regression_total_muts = ggplotRegression(lm(no.muts.in.cluster ~ AADx, data=subclones_sigs))
plot_lm_regression_sbs40 = ggplotRegression(lm(C_SBS40 ~ AADx, data=subclones_sigs))
plot_lm_regression_sbs18 = ggplotRegression(lm(A_SBS18 ~ AADx, data=subclones_sigs))
plot_lm_regression_sbs38 = ggplotRegression(lm(B_SBS38 ~ AADx, data=subclones_sigs))

plot_grid(plot_lm_regression_sbs18, plot_lm_regression_sbs38, plot_lm_regression_sbs40, plot_lm_regression_total_muts, nrow=2)
############################################ Estimate patient-specific mutation rates
### free intercepts
muts.per.year.lmer.with.pt.intercepts <- lmer(C_SBS40 ~ AADx + (1 + AADx | patient), data=subclones_sigs_sub, REML=FALSE)
### patient-specific intercept constrained to 0
muts.per.year.lmer.with.pop.intercept <- lmer(C_SBS40 ~ AADx + (0 + AADx | patient), data=subclones_sigs_sub, REML=FALSE)
### population and patient-specific intercepts constrained to 0
muts.per.year.lmer <- lmer(C_SBS40 ~ 0 + AADx + (0 + AADx | patient), data=subclones_sigs_sub, REML=FALSE)
### ANOVA for model differences
pandoc.table(anova(muts.per.year.lmer, muts.per.year.lmer.with.pop.intercept, muts.per.year.lmer.with.pt.intercepts), split.tables = "Inf")
#----------------------------------------------------------------------------------------------------------------
#                  &nbsp;                     npar    AIC     BIC    logLik   deviance   Chisq   Df   Pr(>Chisq) 
#------------------------------------------- ------ ------- ------- -------- ---------- ------- ---- ------------
#          **muts.per.year.lmer**              3     56694   56712   -28344    56688      NA     NA       NA     
#
# **muts.per.year.lmer.with.pop.intercept**    4     56674   56699   -28333    56666     21.41   1    3.705e-06  
#
# **muts.per.year.lmer.with.pt.intercepts**    6     56670   56707   -28329    56658     8.688   2     0.01298   
#----------------------------------------------------------------------------------------------------------------
### Choose the model with non-zero population intercept
fitted_number_of_muts = (fixef(muts.per.year.lmer)['AADx'] + ranef(muts.per.year.lmer)$patient[,'AADx']) * unique(subclones_sigs_sub[,c('patient','AADx')])$AADx
patients = cbind(unique(subclones_sigs_sub[,c('patient','AADx')]), fitted_number_of_muts)
subclones_sigs_sub = merge(x=subclones_sigs_sub, y=patients[,c('patient','fitted_number_of_muts')], by.x='patient', by.y='patient')
plot_mutation_rates <- ggplot(data=subclones_sigs_sub, aes(x=AADx, y=C_SBS40, color=DISEASE_SUBTYPE)) + geom_segment(aes(x=0, y=0, xend=AADx, yend=fitted_number_of_muts, color=DISEASE_SUBTYPE)) + geom_abline(0, slope=fixef(muts.per.year.lmer)['AADx'], size=1) + geom_point(shape=16, cex=1) + scale_color_manual(values=disease_subtype_cols$color)

relative_mutation_rate = (fixef(muts.per.year.lmer.with.pop.intercept)['AADx'] + ranef(muts.per.year.lmer.with.pop.intercept)$patient[,'AADx'])
pt_mutation_rate = relative_mutation_rate + abs(min(relative_mutation_rate))
patients = cbind(unique(subclones_sigs_sub[,c('patient','AADx','DISEASE_SUBTYPE')]), pt_mutation_rate)
patients$DISEASE_SUBTYPE = factor(patients$DISEASE_SUBTYPE, levels=disease_subtype_cols$disease_subtype)
############################################ Estimate time of emergence for MRCA
subclones_sigs_sub$C_SBS40_scaled = scale(subclones_sigs_sub$C_SBS40, center=FALSE)
scale_factor_sbs40 <- attr(subclones_sigs_sub$C_SBS40_scaled, which = "scaled:scale")
age_lmer_no_intercepts <- lmer(AADx ~ 0 + C_SBS40_scaled + (0 + C_SBS40_scaled | patient), data=subclones_sigs_sub, REML=FALSE)
age_lmer_pop_intercept <- lmer(AADx ~ C_SBS40_scaled + (0 + C_SBS40_scaled | patient), data=subclones_sigs_sub, REML=FALSE)
pandoc.table(anova(age_lmer_no_intercepts, age_lmer_pop_intercept), split.tables = "Inf")
#
#-------------------------------------------------------------------------------------------------
#           &nbsp;             npar    AIC     BIC    logLik   deviance   Chisq   Df   Pr(>Chisq) 
#---------------------------- ------ ------- ------- -------- ---------- ------- ---- ------------
# **age_lmer_no_intercepts**    3     40944   40962   -20469    40938      NA     NA       NA     
#
# **age_lmer_pop_intercept**    4     38119   38144   -19056    38111     2827    1        0      
#-------------------------------------------------------------------------------------------------

# Pick the model with non-zero population intersept
mrca_estimates <- estimateMRCA(subclones_sigs_sub, age_lmer_pop_intercept)
mrca_estimates = merge(x=mrca_estimates, y=unique(mdata[,c('PATIENT','DISEASE_SUBTYPE')]), by.x='ID', by.y='PATIENT')
mrca_estimates = merge(x=mrca_estimates, y=disease_subtype_cols, by.x='DISEASE_SUBTYPE', by.y='disease_subtype')
mrca_estimates_sub = mrca_estimates[mrca_estimates$ID %in% unique(mdata$PATIENT[mdata$SAMPLE_TYPE %in% c('diagnosis','reresection','t-resection')]),]
mrca_estimates_sub = mrca_estimates_sub[order(mrca_estimates_sub$DISEASE_SUBTYPE, mrca_estimates_sub$Time.lag),]
id_order = mrca_estimates_sub$ID
mrca_estimates_sub = merge(x=mrca_estimates_sub, y=unique(subclones_sigs_sub[subclones_sigs_sub$descr=='trunk',c('patient','C_SBS40')]), by.x='ID', by.y='patient')
mrca_estimates_sub$C_SBS40 = round(mrca_estimates_sub$C_SBS40)
mrca_estimates_sub$label = apply(cbind(apply(cbind('(',mrca_estimates_sub[,'C_SBS40'], ')'), 1, paste, collapse="") , mrca_estimates_sub[,'ID']), 1, paste, collapse=" ")
mrca_estimates_sub$C_SBS40 = NULL
mrca_estimates_sub = mrca_estimates_sub[match(id_order, mrca_estimates_sub$ID),]
pdf('raw/figure_2b.pdf', width=12, height=15)
plotMRCA(mrca_estimates_sub)
dev.off()
