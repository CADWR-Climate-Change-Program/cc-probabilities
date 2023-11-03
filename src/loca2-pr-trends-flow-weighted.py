# %%
import glob, os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

# %%
def linear_regression(x, y, plot=False, ylim=None,title=None, print=False):

    X = x.to_numpy().reshape(-1, 1) 
    reg = LinearRegression(fit_intercept=True)
    _ = reg.fit(X, y)
    a = reg.coef_[0]
    b = reg.intercept_
    
    if plot:
        y_hat = a * x + b
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(1, 1, 1)
        _ = plt.plot(x, y, 'o')
        _ = plt.plot(x, y_hat, '-', label='Linear regression')
        _ = ax.set_xlabel('Year (30y Start)')
        _ = ax.set_ylabel('Precip')
        _ = ax.legend()
        if ylim: _ = ax.set_ylim(ylim)
        if title: _ = ax.set_title(title)
        if print: plt.savefig('../figures/lm_fits/'+mvs+'_lm'+".png",dpi=300,bbox_inches='tight')
    return a, b

# %% diff
def difference(df, column, position=0):
    df.reset_index(inplace=True,drop=True)
    return df[column] - df[column].iloc[position]


# %%
files = glob.glob('../data/loca2-projections-annual-flow/*.csv')

model_ssp_variant = pd.DataFrame({
    'm':[os.path.basename(i).split('_')[0] for i in files],
    'v':[os.path.basename(i).split('_')[2] for i in files],
    's':[os.path.basename(i).split('_')[1] for i in files]
})
model_ssp_variant = model_ssp_variant.loc[model_ssp_variant.s!='historical']
model_ssp_variant_counts = model_ssp_variant.groupby(['m','s'],as_index=False).count()

#%%
projections_data = pd.DataFrame()
for f in files:
    f_name = os.path.basename(f)
    model = f_name.split('_')[0]
    variant = f_name.split('_')[2] #1
    ssp = f_name.split('_')[1] #2

    df = pd.read_csv(f)
    df = df.groupby('Year',as_index=False).aggregate({'Pr (mm)':np.sum,'Tave (degC)':np.mean})
    df.columns = ['y','pr','tavg']
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
            projections_data = pd.concat([projections_data,df_copy],axis=0)
    else:
        df.insert(2,'ssp',value=ssp)
        df.insert(0,'mvs',value=model+'_'+variant+'_'+ssp)
        projections_data = pd.concat([projections_data,df],axis=0)

# %% window avg
window=30
projections = projections_data.copy()
projections.reset_index(inplace=True,drop=True)
projections['pr_roll'] = projections.groupby('mvs',as_index=False)['pr'].rolling(window).mean()['pr']
projections['tavg_roll'] = projections.groupby('mvs',as_index=False)['tavg'].rolling(window).mean()['tavg']
projections.dropna(axis = 0, how = 'any', inplace = True)
projections = projections.loc[projections.y>=1981+window]
mvs_list = projections.mvs.unique()

# %%
fits = pd.DataFrame()
for mvs in mvs_list:
    data = projections.loc[projections['mvs']==mvs]
    a, b = linear_regression(data.y, data['pr_roll'], ylim=(850,1300),title=mvs,print=False)
    df = pd.DataFrame({'mvs':mvs,
                       'm':mvs.split('_')[0],
                       'v':mvs.split('_')[1],
                       's':mvs.split('_')[2],
                       'slope':a,'intercept':b,
                       'base':a*(1992+window) + b, 
                       '2043':a*(2028+window) + b, 
                       '2050':a*(2035+window) + b, 
                       '2070':a*(2055+window) + b},
                       index=[mvs])
    df.insert(6,'pr_change_2043',value=round((df['2043']-df['base'])/df['base']*100,1))
    df.insert(7,'pr_change_2070',value=round((df['2070']-df['base'])/df['base']*100,1))
    df.insert(8,'pr_change_2050',value=round((df['2050']-df['base'])/df['base']*100,1))
    fits = pd.concat([fits,df])

fits.sort_values(by=['m','s'],inplace=True)
fits = fits.loc[fits.s!='historical']

# %%
fits.to_csv('../processed/loca2_lmfits_1981_cv-flow-weighted.csv',index=False)

# %%
t_diffs = pd.DataFrame(projections.groupby(['mvs']).apply(lambda x: difference(x,'tavg_roll',0)))
t_diffs.columns = np.arange(1981+window,2101)
t_diffs = t_diffs.melt(ignore_index=False).reset_index()
t_diffs.columns = ['mvs','y','tavg_roll_diff']
projections = pd.merge(projections,t_diffs,on=['mvs','y'])

# %% 
for t in [1,2,3,4,5]:
    diff_levels = projections.groupby(['mvs']).apply(lambda x: x['y'].values[np.searchsorted(x['tavg_roll_diff'],t)-1])
    if t==1:
        warming_levels = pd.DataFrame(diff_levels)
        warming_levels.columns = [t]
    else:
        warming_levels[t] = diff_levels.values
warming_levels.reset_index(inplace=True)

warming_levels.to_csv('../processed/loca2_warming_levels_cv-flow-weighted.csv',index=False)
