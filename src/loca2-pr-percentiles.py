# %%
import itertools, glob, os
from datetime import timedelta

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mat
plt.style.use("seaborn")
import seaborn as sns
import pandas as pd
from sklearn.linear_model import LinearRegression

sns.set_theme(style="ticks")
mat.rcParams['figure.dpi']= 300

# %%
# files = glob.glob('../data/loca2-projections-daily-flow/*.csv')
files = glob.glob('../data/loca2-projections-daily-basin/*_01-19.csv')
region_title = 'shasta'

model_ssp_variant = pd.DataFrame({
    'm':[os.path.basename(i).split('_')[0] for i in files],
    'v':[os.path.basename(i).split('_')[2] for i in files],
    's':[os.path.basename(i).split('_')[1] for i in files]
})
model_ssp_variant = model_ssp_variant.loc[model_ssp_variant.s!='historical']
model_ssp_variant_counts = model_ssp_variant.groupby(['m','s'],as_index=False).count()
# df_all = model_ssp_variant_full.merge(model_ssp_variant.drop_duplicates(), on=['m','v','s'], 
#                    how='left', indicator=True)

# %%
quantiles=[0.5,0.8,0.9,0.95,0.99]

#%%
projections_data = pd.DataFrame()

for f in files:
    f_name = os.path.basename(f)
    model = f_name.split('_')[0]
    variant = f_name.split('_')[2] #1
    ssp = f_name.split('_')[1] #2

    df = pd.read_csv(f,usecols=['Year','Month','Day','Pr (mm)','Tasmax (degC)','Tasmin (degC)'])
    df['tavg'] = (df['Tasmax (degC)']+df['Tasmin (degC)'])/2
    df.columns = ['y','m','d','pr','tmax','tmin','tavg']
    df.insert(0,'model',value=model)
    df.insert(1,'variant',value=variant)

    if ssp =='historical':
        ssps = model_ssp_variant.loc[
                    (model_ssp_variant.m==model) & (model_ssp_variant.v==variant)
                ]
        for each in ssps['s']:
            df_copy = df.copy()
            df_copy.insert(2,'ssp',value=each)
            df_copy.insert(0,'mvs',value=model+'_'+variant+'_'+each)
            # df_copy.insert(0,'ms',value=model+'_'+each)
            projections_data = pd.concat([projections_data,df_copy],axis=0)
    else:
        df.insert(2,'ssp',value=ssp)
        df.insert(0,'mvs',value=model+'_'+variant+'_'+ssp)
        # df.insert(0,'ms',value=model+'_'+ssp)
        projections_data = pd.concat([projections_data,df],axis=0)

projections_data.reset_index(drop=True,inplace=True)
mvs_list = projections_data.mvs.unique()

# %% historical reference
ppt_hist_percentiles = pd.DataFrame()
for mvs in mvs_list:
    data = projections_data.loc[(projections_data.mvs==mvs) & (projections_data.y<=1999)]
    data_ppt_nonzero = data.loc[data.pr>=0.01]
    data_ppt_percentiles = pd.DataFrame({
        'mvs':mvs,
        'quantile':quantiles,
        'ppt':np.quantile(data_ppt_nonzero.pr,q=quantiles)
    })
    ppt_hist_percentiles = pd.concat([ppt_hist_percentiles,data_ppt_percentiles])
ppt_hist_percentiles.reset_index(inplace=True,drop=True)
ppt_hist_percentiles.to_csv('../processed/loca2_ppt_hist_percentiles_{}.csv'.format(region_title),index=False)

# %%
# ppt_hist_percentiles = pd.read_csv('../processed/loca2_ppt_hist_percentiles.csv')

# %% plot
f, ax = plt.subplots(figsize=(4,4),facecolor=None)
sns.boxplot(data=ppt_hist_percentiles,x='quantile',y='ppt',ax=ax,linewidth=1)
sns.despine(trim=True)
# ax.grid()
# plt.yscale('log')

# %%
ppt_hist_qtotals = pd.DataFrame()
for mvs in mvs_list:
    for q in quantiles:
            q_ppt_lvl = float(ppt_hist_percentiles.loc[(ppt_hist_percentiles.mvs==mvs) & 
                                                       (ppt_hist_percentiles['quantile']==q)]['ppt'])
            ppt_qwet = projections_data.loc[(projections_data.mvs==mvs) & 
                                            (projections_data.y<=1999) & 
                                            (projections_data.pr>=q_ppt_lvl)]['pr']
            ppt_qother = projections_data.loc[(projections_data.mvs==mvs) & 
                                              (projections_data.y<=1999) & 
                                              (projections_data.pr<q_ppt_lvl) & 
                                              (projections_data.pr>=0.01)]['pr']
            data_ppt = pd.DataFrame({
                'mvs':mvs,
                'quantile':q,
                'count_ppt':np.count_nonzero(ppt_qwet) + np.count_nonzero(ppt_qother),
                'ppt_qwet':np.sum(ppt_qwet),
                'ppt_qother':np.sum(ppt_qother),
                'ppt_qtotal':np.sum(ppt_qwet) + np.sum(ppt_qother)
            },index=[0])
            ppt_hist_qtotals = pd.concat([ppt_hist_qtotals,data_ppt])
ppt_hist_qtotals.reset_index(inplace=True,drop=True)
ppt_hist_qtotals.to_csv('../processed/loca2_ppt_hist_qtotals_{}.csv'.format(region_title),index=False)

# %%
# ppt_hist_qtotals = pd.read_csv('../processed/loca2_ppt_hist_qtotals.csv')

# %%
# yr_step = 5
# window = 50
# ppt_fut_qtotals = pd.DataFrame()
# for mvs in mvs_list:
#     for q in quantiles[-2:-1]:
#         q_ppt_lvl = float(ppt_hist_percentiles.loc[(ppt_hist_percentiles.mvs==mvs) & 
#                                                     (ppt_hist_percentiles['quantile']==q)]['ppt'])
#         for y in np.arange(2050,2100+yr_step,yr_step):
#             ppt_qwet = projections_data.loc[(projections_data.mvs==mvs) & 
#                                             (projections_data.y>=(y-window)) & 
#                                             (projections_data.y<y) & 
#                                             (projections_data.pr>=q_ppt_lvl)]['pr']
#             ppt_qother = projections_data.loc[(projections_data.mvs==mvs) & 
#                                             (projections_data.y>=(y-window)) & 
#                                             (projections_data.y<y) & 
#                                             (projections_data.pr<q_ppt_lvl) & 
#                                             (projections_data.pr>=0.01)]['pr']
#             data_ppt = pd.DataFrame({
#                 'mvs':mvs,
#                 'quantile':q,
#                 'window':y,
#                 'count_ppt':np.count_nonzero(ppt_qwet) + np.count_nonzero(ppt_qother),
#                 'count_qwet':np.count_nonzero(ppt_qwet),
#                 'count_qother':np.count_nonzero(ppt_qother),
#                 'ppt_qwet':np.sum(ppt_qwet),
#                 'ppt_qother':np.sum(ppt_qother),
#                 'ppt_qtotal':np.sum(ppt_qwet) + np.sum(ppt_qother)
#             },index=[0])
#             ppt_fut_qtotals = pd.concat([ppt_fut_qtotals,data_ppt])

# ppt_fut_qtotals.reset_index(inplace=True,drop=True)
# ppt_fut_qtotals.to_csv('../processed/loca2_ppt_fut_qtotals.csv',index=False)

# %%
# ppt_fut_qtotals = pd.read_csv('../processed/loca2_ppt_fut_qtotals.csv')

# %% plot
# sns.lineplot(data=ppt_fut_qtotals,x='window',y='count_qwet')
# sns.lineplot(data=ppt_fut_qtotals,x='window',y='count_qother')
# sns.lineplot(data=ppt_fut_qtotals,x='window',y='ppt_qtotal')
# sns.lineplot(data=ppt_fut_qtotals,x='window',y='ppt_qwet')

# %%
warming_levels = pd.read_csv('../processed/loca2_warming_levels.csv')

#%%
window=50
ppt_warming_levels_qtotals = {}
for t in np.arange(1,6,1):
    mvs_warming_levels = warming_levels.loc[warming_levels[str(t)]<2100]
    df = pd.DataFrame()

    for index, row in mvs_warming_levels.iterrows():
        for q in quantiles[-2:-1]:

            mvs = row['mvs']
            y = min(2100,row[str(t)]+10)

            try: 
                q_ppt_lvl = float(ppt_hist_percentiles.loc[(ppt_hist_percentiles.mvs==mvs) & 
                                                        (ppt_hist_percentiles['quantile']==q)]['ppt'])
            except:
                print(mvs)
                continue
            
            ppt_qwet = projections_data.loc[(projections_data.mvs==mvs) & 
                                            (projections_data.y>=(y-window)) & 
                                            (projections_data.y<y) & 
                                            (projections_data.pr>=q_ppt_lvl)]['pr']
            
            ppt_qother = projections_data.loc[(projections_data.mvs==mvs) & 
                                            (projections_data.y>=(y-window)) & 
                                            (projections_data.y<y) & 
                                            (projections_data.pr<q_ppt_lvl) & 
                                            (projections_data.pr>=0.01)]['pr']
            
            data_ppt = pd.DataFrame({
                'mvs':mvs,
                'quantile':q,
                't':t,
                'window':y,
                'count_ppt':np.count_nonzero(ppt_qwet) + np.count_nonzero(ppt_qother),
                'count_qwet':np.count_nonzero(ppt_qwet),
                'count_qother':np.count_nonzero(ppt_qother),
                'ppt_qwet':np.sum(ppt_qwet),
                'ppt_qother':np.sum(ppt_qother),
                'ppt_qtotal':np.sum(ppt_qwet) + np.sum(ppt_qother)
            },index=[0])

            df = pd.concat([df,data_ppt])

    ppt_warming_levels_qtotals[t] = df

# %%
q = quantiles[-2]
ppt_warming_levels_hist = {}
for t in np.arange(1,6,1):
    df = pd.merge(ppt_warming_levels_qtotals[t],
                                          ppt_hist_qtotals.loc[ppt_hist_qtotals['quantile']==q],
                                          how='inner',on='mvs',
                                        suffixes=['_fut','_hist'])
    df['ppt_change'] = (df['ppt_qtotal_fut'] - df['ppt_qtotal_hist'])/df['ppt_qtotal_hist']
    df['qwet_change'] = (df['ppt_qwet_fut'] - df['ppt_qwet_hist'])/df['ppt_qtotal_hist']
    df['qother_change'] = (df['ppt_qother_fut'] - df['ppt_qother_hist'])/df['ppt_qtotal_hist']
    ppt_warming_levels_hist[t] = df.copy()

# %%
# ppt_warming_levels_hist[4]

# %%
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,8), facecolor=None,sharex=True,sharey=True)

plt.subplots_adjust(hspace=0.15)
for j,ax in enumerate(axes.flatten()):
    t=j+1
    sns.regplot(data=ppt_warming_levels_hist[t],y='ppt_change',x='qwet_change',ax=ax,label='Wettest 5% Days')
    sns.regplot(data=ppt_warming_levels_hist[t],y='ppt_change',x='qother_change',marker='x',ax=ax,label='Other Days')
    ax.set_title('T Level = {}C'.format(t))
    ax.set_ylabel('Total Precipitation Change (%)') if j in [0,2] else ax.set_ylabel('')
    ax.set_xlabel('Quantile Precipitation Change\n(as percent of historical PPT)') if j in [2,3] else ax.set_xlabel('')
    ax.axvline(0,color='k',linestyle='dashed',lw=0.5)
    ax.axhline(0,color='k',linestyle='dashed',lw=0.5)
    ax.set_xticks(np.arange(-0.2,0.2+0.1,0.1))
    ax.set_yticks(np.arange(-0.2,0.2+0.1,0.1))
    ax.set_ylim(-0.2,0.2)
    ax.set_xlim(-0.2,0.2)
    ax.grid()
    if j==0: ax.legend(loc='lower right')
    # if j==3: ax.legend(loc='lower right', ncol=2, bbox_to_anchor=(0.4,-0.5,0,-1))

# %%
for t in np.arange(1,6,1):
    ppt_warming_levels_hist[t].to_csv('../processed/loca2_ppt_warming_levels_{}C_{}.csv'.format(t,region_title),index=False)

# %%
for t in np.arange(1,6,1):
    ppt_warming_levels_hist[t] = pd.read_csv('../processed/loca2_ppt_warming_levels_{}C_{}.csv'.format(t,'tuolumne'))
# %%
