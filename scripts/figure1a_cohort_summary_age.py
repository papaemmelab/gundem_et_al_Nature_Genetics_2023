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
dic = {'A1':1,'A2':2,'A3':3}
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
cohort_summary_hm = cohort_summary_hm[sample_order]
cohort_summary_bk = cohort_summary_bk[sample_order]
cohort_summary_hm = cohort_summary_hm.drop('-', axis=1, errors='ignore')
cohort_summary_bk = cohort_summary_bk.drop('-', axis=1, errors='ignore').fillna(8)
fig, ax = plt.subplots(figsize=(25,8))  
hm = sns.heatmap(cohort_summary_bk, cmap=['white'], cbar=False, linewidths=.5, square=True, xticklabels=True)
for (x,y), val in np.ndenumerate(cohort_summary_hm):        
    if val == 1:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#E5E5E5', edgecolor='black',lw=.1))
    elif val == 2:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#A6A6A6', edgecolor='black',lw=.1))
    elif val == 3:
        ax.add_artist(plt.Rectangle((y+.05, x+.05), 0.85,0.85, facecolor='#4D4D4D', edgecolor='black',lw=.1))
dic = {'A1':'A1','2':'A2','3':'A3'}
lg1 = mlines.Line2D([], [], color='#E5E5E5', marker='s', linestyle='None', markersize=10, label='<18m')
lg2 = mlines.Line2D([], [], color='#A6A6A6', marker='s', linestyle='None', markersize=10, label='>=18m & <12y')
lg3 = mlines.Line2D([], [], color='#4D4D4D', marker='s', linestyle='None', markersize=10, label='>=12y')
plt.legend(handles=[lg1, lg2, lg3],loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3)
plt.tight_layout()
plt.savefig('/Users/gundemg/Documents/projects/neuroblastoma/triple_callers/analysis/scripts_final/raw/figure_1a__cohort_characteristics_age.pdf')
