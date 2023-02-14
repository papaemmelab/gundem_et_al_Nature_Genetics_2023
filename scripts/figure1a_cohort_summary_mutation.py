import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns; sns.set()
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import itertools
import matplotlib.lines as mlines
## path to data
path = sys.argv[1]
## import data
## 3 columns |sample_id | patient_id | info
data = pd.read_csv(path, header=None, sep=' ')
data.columns = ['sid', 'pid', 'info']
# define info groups in order 
dic = {'MYCN':1,'ATRX':2,'TERT':3,'MDM2-CDK4':4,'NUMERIC_CNV':5,'SEG_CNV':6}
sample_order = list(data['pid'].unique())
def create_cohort_summary(df, info_column):
    cohort_summary = df[['sid', 'pid','info']]
    cohort_summary = cohort_summary.groupby(['pid'])[info_column].agg(lambda x: ','.join(x)).reset_index()
    return cohort_summary
## Do data transformation
cohort_summary_hm = create_cohort_summary(data, 'info')
cohort_summary_expand = pd.concat([cohort_summary_hm['pid'], cohort_summary_hm['info'].str.split(',', expand=True)], axis=1).set_index('pid').T
## replace string values wiht int values using dic
cohort_summary_hm = cohort_summary_expand.replace(dic)
### reverse index
cohort_summary_hm = cohort_summary_hm.reindex(index=cohort_summary_hm.index[::-1]).reset_index(drop=True)
### sort each column individually and reverse
cohort_summary_hm = cohort_summary_hm.apply(lambda x: x.sort_values().values).reindex(index=cohort_summary_hm.index[::-1])
cohort_summary_bk = cohort_summary_hm.fillna(8).applymap(lambda x: (x/x)*8)
cohort_summary_hm = cohort_summary_hm.drop('-', axis=1, errors='ignore')
cohort_summary_bk = cohort_summary_bk.drop('-', axis=1, errors='ignore').fillna(8)
cohort_summary_hm = cohort_summary_hm[sample_order]
cohort_summary_bk = cohort_summary_bk[sample_order]
fig, ax = plt.subplots(figsize=(25,8))  
hm = sns.heatmap(cohort_summary_bk, cmap=['white'], cbar=False, linewidths=.5, square=True, xticklabels=True)
for (x,y), val in np.ndenumerate(cohort_summary_hm):        
    if val == 1:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#009E73', edgecolor='black',lw=.1))
    elif val == 2:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#737373', edgecolor='black',lw=.1))
    elif val == 3:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#9ad0f3', edgecolor='black',lw=.1))
    elif val == 4:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#CC79A7', edgecolor='black',lw=.1))
    elif val == 5:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#0072B2', edgecolor='black',lw=.1))
    elif val == 6:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#e79f00', edgecolor='black',lw=.1))
#ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
dic = {'MYCN':1,'ATRX':2,'TERT':3,'MDM2-CDK4':4,'NUMERIC_CNV':5,'SEG_CNV':6}
lg1 = mlines.Line2D([], [], color='#009E73', marker='s', linestyle='None', markersize=10, label='MYCN')
lg2 = mlines.Line2D([], [], color='#737373', marker='s', linestyle='None', markersize=10, label='ATRX')
lg3 = mlines.Line2D([], [], color='#9ad0f3', marker='s', linestyle='None', markersize=10, label='TERT')
lg4 = mlines.Line2D([], [], color='#CC79A7', marker='s', linestyle='None', markersize=10, label='MDM2-CDK4')
lg5 = mlines.Line2D([], [], color='#0072B2', marker='s', linestyle='None', markersize=10, label='SEG_CNV')
lg6 = mlines.Line2D([], [], color='#e79f00', marker='s', linestyle='None', markersize=10, label='NUMERIC_CNV')
plt.legend(handles=[lg1, lg2, lg3, lg4, lg5, lg6],loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=6)
plt.tight_layout()
plt.savefig('/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final/raw/figure_1a__cohort_characteristics_mutation.pdf')
