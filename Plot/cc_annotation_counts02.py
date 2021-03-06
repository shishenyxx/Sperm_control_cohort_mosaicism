# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 14:03:21 2020

@author: Martin
"""

#Import modules

import sys, os, gzip
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import statsmodels.stats.api as sms
import scipy.stats as st

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()

#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions to annotate additional columns for categories

'''table is csv, open with regular pandas commands. a working table of mosaic
variants .IS_MOSAIC == True and of those not in Blastocysts .COHORT !=
'Blastocyst is generated as ususally done.'''

def annotate_sets(df):
    '''use information provided by XY to annotate categories that are needed
    for plotting and subdivision of the table. cutoffs will be explained at the
    annotations inplace.'''
    
    df['SET_GENE'] = df['NCBI_gene_merged']
    df['SET_EXON'] = df['NCBI_exons_merged']
    df['SET_IMPACT_NON_GENIC'] = (df.define_Pathogenic_prob > 0.95).apply(int)
    df['SET_IMPACT_GENE'] = ((df.CADD13_PHRED > 25) |
                             (df.Loss_of_function == True)).apply(int)
    df['SET_IMPACT'] = ((df.SET_IMPACT_NON_GENIC == True) |
                        (df.SET_IMPACT_GENE == True)).apply(int)

#Function to make tables for plotting (descriptive)

def make_summary_table_with_CI(df, category='DUMMY'):
    '''make a summary table for counts that are available for plotting and also
    allow to annotate the 95% CI. relies on restricting input by category, so
    that the numbers are plotted correctly. Dummy will include all.'''
    
    df_ = df.copy()
    
    df_['DUMMY'] = 1
    
    df_ = df_[df_[category] == True]
    
    df_ = df_.groupby(['INDIVIDUAL', 'MUTSIG_CLASS']).count().reset_index()
    
    
    df0 = df.groupby(['INDIVIDUAL', 'MUTSIG_CLASS']).count().reset_index()
    df0 = df0[['INDIVIDUAL', 'MUTSIG_CLASS']]
    df0['POS1'] = 0
    
    df_ = pd.concat([df_, df0], sort=True)
    df_ = df_.groupby(['INDIVIDUAL', 'MUTSIG_CLASS']).sum().reset_index()
    
    df_s = df_.copy()
    
    df_ = df_.groupby(['MUTSIG_CLASS']).POS1.describe()
    
    cis = st.t.interval(0.95, np.add(df_['count'], -1), loc=df_['mean'],
                        scale=df_['std']/df_['count']**0.5)
    
    df_['95CI_LOWER'] = cis[0]
    df_['95CI_UPPER'] = cis[1]
    
    return df_, df_s


def make_estimation_table(df):
    
    df_ = df[(df.MUTSIG_CLASS.isin(['1_sperm', '2_shared'])) &
             (df.SET_IMPACT_GENE == True)]
    
    df_ = df_.groupby('INDIVIDUAL').count().reset_index()
    
    df0 = df.groupby(['INDIVIDUAL']).count().reset_index()
    df0 = df0[['INDIVIDUAL']]
    df0['POS1'] = 0
    
    df_ = pd.concat([df_, df0], sort=True)
    df_ = df_.groupby(['INDIVIDUAL']).sum().reset_index()
    
    df_ = df_.POS1.describe()
    
    cis = st.t.interval(0.95, np.add(df_['count'], -1), loc=df_['mean'],
                        scale=df_['std']/df_['count']**0.5)
    
    df_['95CI_LOWER'] = cis[0]
    df_['95CI_UPPER'] = cis[1]
    
    df_ = df_[['mean', '95CI_LOWER', '95CI_UPPER']]
    
    #Data provided by XY: {ALL:68556364bp; HI:17651035; SFARI:3185225
    
    df_100 = df_ * 100
    df_100_hi = df_100 * 17651035 / 68556364
    df_100_sfari = df_100 * 3185225 / 68556364
    
    df_100['SUBGROUP'] = 'ALL_GENES'
    df_100_hi['SUBGROUP'] = 'HI_GENES'
    df_100_sfari['SUBGROUP'] = 'SFARI_HI'
    
    df_complete = pd.concat([df_100, df_100_hi, df_100_sfari], axis=1).T\
                                                                .reset_index()
    
    return df_complete

#Functions to prepare for permutation enrichment analysis

def assemble_simulations(path):
    '''new version where shuffle_summary for everything is already generated by
    XY and allows to load with a category.'''
    
    df = pd.read_csv(path, sep='\t')
    
    df.rename({'hg19_early_timing_PMID19966280': 'early_timing_PMID19966280',
               'hg19_late_timing_PMID19966280': 'late_timing_PMID19966280'},
              axis=1, inplace=True)
    
    df.sort_values(by='GROUP', inplace=True)
    
    groups = df.GROUP.unique().tolist()
    mscs = ['1_sperm', '2_shared', '3_soma_y', '4_soma_o']
    cats = ['cc', 'gnomAD', 'ssc']
    
    mscs_cats = [(cat, msc) for cat in cats for msc in mscs]
    
    dictionary = {grp: msc_cat for grp, msc_cat in zip(groups, mscs_cats)}
    
    df['CATEGORY_ShufgnoSSC'] = df.apply(lambda row: dictionary[row.GROUP][0],
                                         axis=1)
    
    df['MUTSIG_CLASS'] = df.apply(lambda row: dictionary[row.GROUP][1], axis=1)
    
    df = df.groupby(['CATEGORY_ShufgnoSSC', 'MUTSIG_CLASS'])\
           .describe(percentiles=[0.025, 0.975])
           
    df.reset_index(inplace=True)
    
    return df


def load_actual_fractions(path):
    '''percentage_of_features.csv; no changing required. use the PERCENTAGE
    column for the fraction readout.'''
    
    return pd.read_csv(path)


def process_lymph_semi(path):
    '''import .csv with all data from lymphona or seminoma.'''
    
    df = pd.read_csv(path)
    
    df.rename({'hg19_early_timing_PMID19966280': 'early_timing_PMID19966280',
               'hg19_late_timing_PMID19966280': 'late_timing_PMID19966280'},
              axis=1, inplace=True)
    
    df['One'] = 1
    df['Zonk'] = 'zonk'
    df_ = df.groupby('Zonk').sum()
    df_ = df_.apply(lambda row: row / row['One'], axis=1)
    df_ = df_.stack().reset_index().sort_values('level_1')
    
    return df_


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting
def plot_category_counts(df):
    '''simplified and optimized version to plot for this specific problem. Can
    be more customizable if not using a split first frame and giving the cats
    of interest as an argument.'''
    
    categories = ('IS_MOSAIC', 'SET_EXON', 'SET_IMPACT_GENE')
    ranges = (-0.3, -0.1, 0.1, 0.3)
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']    
    fig = plt.figure()
    
    axs = []
    
    for i in [1,4]: #suplots start at 1, not 0
        axs.append(fig.add_subplot(2,3,i))
    for i in range (2,4):
        axs.append(fig.add_subplot(1,3,i))
    
    for i in range(len(axs)):
        
        if i == 0:
            cat = categories[0]
        elif i == 1:
            cat = categories[0]
        else:
            cat = categories[i-1]
            
        df_, df_s = make_summary_table_with_CI(df, cat)
        df_['PLOTTER'] = cat
        df_s['PLOTTER'] = cat
        df_s.sort_values(by='MUTSIG_CLASS', inplace=True)
        
        sns.stripplot(x='PLOTTER', y='POS1', hue='MUTSIG_CLASS', data=df_s,
                      dodge=True, ax=axs[i], alpha=0.5, s=4, palette=colors)
        axs[i].errorbar(x=ranges, y=df_['mean'],
                        yerr=np.subtract(df_['95CI_UPPER'], df_['mean']),
                        marker='d', ms=10, mfc='k',mec ='k', ecolor='k', ls='',
                        capsize=5)
    
    plt.subplots_adjust(hspace=0.1)
    
    axs[0].set_ylim(50,500)
    axs[0].set_yticks([50, 250, 500])
    axs[0].set_xticks([])
    axs[0].set_xlabel('')
    axs[0].set_xticklabels('')
    axs[0].set_ylabel('Number of Mosaic Variants')
    
    axs[1].set_ylim(0,50)
    axs[1].set_yticks([0, 25, 50])
    axs[1].set_xlabel('All Variants')
    axs[2].set_ylim(-0.5, 20)
    axs[2].set_yticks([0, 5, 10, 15, 20])
    axs[2].set_xlabel('Exonic Variants')
    axs[3].set_ylim(-0.2, 6)
    axs[3].set_yticks([0, 2, 4, 6])
    axs[3].set_xlabel('CADD > 25/LoF')
    
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    
    for i in range(1,4):
        axs[i].get_legend().remove()
        axs[i].set_ylabel('')
        axs[i].set_xticks([])
        axs[i].set_xticklabels('')
        sns.despine(bottom=True, offset=5, trim=True, ax=axs[i])
    
    plt.show()


def plot_estimates(df):
    
    df = make_estimation_table(df)
    
    sns.pointplot()
    #fix, as sns.despine does not work with categorical other than sns
    #same is true for the post-labeling
    plt.errorbar(x=[0,1,2], y=df['mean'],
                 yerr=np.subtract(df['mean'], df['95CI_LOWER']), marker='d',
                 ms=10, mfc='k',mec ='k',ls='', capsize=5, ecolor='k')
    
    plt.ylim(0,50)
    plt.ylabel('Number of Variants per 100 Men')
    sns.despine(offset=5, trim=True)
    plt.xticks(ticks=[0,1,2], labels=df['SUBGROUP'].tolist(), rotation=45)
    
    plt.xticks()
    
    plt.show()

#------------------------------------------------------------------------------
def plot_anno_sims(df, sims, category='cc'):
    '''input is the mos_nob table and the assembled sims results.'''
    
    sims = sims[sims.CATEGORY_ShufgnoSSC == category]
    sims.drop('CATEGORY_ShufgnoSSC', axis=1, inplace=True)
    sims.set_index('MUTSIG_CLASS', inplace=True)
    
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']    
    
    columns = list(set([col[0] for col in sims.columns]))
    columns.sort()
    
    ranges = make_ranges(columns)
    
    df = df.copy()[columns + ['MUTSIG_CLASS', 'IS_MOSAIC']]
    
    df['wgEncodeRegTfbsClusteredV3'] = df.wgEncodeRegTfbsClusteredV3.fillna(0)\
                                         .apply(lambda x: x != 0)
    df['wgEncodeRegDnaseClusteredV3'] = df.apply(lambda row:
                                row.wgEncodeRegDnaseClusteredV3 == True,
                                                 axis=1)
    
    df_ = df.groupby('MUTSIG_CLASS').sum()
    
    #is_mosaic was all 1, so effectively replaces a SUM column
    df_ = df_.apply(lambda row: row / row['IS_MOSAIC'], axis=1)
    df_.reset_index(inplace=True)
    
    df_ = df_.groupby(['IS_MOSAIC', 'MUTSIG_CLASS']).sum().stack()\
                                                    .reset_index()
    df_.sort_values(by=['level_2', 'MUTSIG_CLASS'], inplace=True)
    
    sims_ = sims.unstack().reset_index()\
                          .sort_values(by=['level_0', 'MUTSIG_CLASS'])
    means = sims_[sims_.level_1 == 'mean'][0]
    errors = sims.stack().stack().unstack(1).reset_index()\
                         .sort_values(by=['level_1', 'MUTSIG_CLASS'])\
                         [['2.5%', '97.5%']]
    minus = abs(errors['2.5%'].values - means.values)
    plus = abs(means.values - errors['97.5%'].values)
    errors_plotting = [minus, plus]
    
    ecolors = ['0.5' if((value >= error_lst[0]) & (value <= error_lst[1]))
                     else 'r'
               for value, error_lst in zip(df_[0].values, errors.values)]
    
    #actual plotting
    sns.stripplot(x='level_2', y=0, hue='MUTSIG_CLASS', data=df_,
                  palette=colors, dodge=True, jitter=False)
    
    plt.errorbar(x=ranges, y=means, yerr=errors_plotting, linestyle='',
                 marker='_', mfc='0.5', mec='0.5', alpha=0.5, elinewidth=5,
                 ecolor=ecolors)
    
    plt.ylim(0,0.71)
    plt.xlabel('')
    plt.ylabel('Fraction of Variants')
    sns.despine(offset=5, trim=True)
    plt.xticks(rotation=45, ha='right')
    
    ax = plt.gca()
    
    ax.get_legend().remove()
    
    plt.show()
    
def make_ranges(columns):
    
    ranges = []
    
    for i in range(len(columns)):
        ranges.append(i - 0.3)
        ranges.append(i - 0.1)
        ranges.append(i + 0.1)
        ranges.append(i + 0.3)
        
    return ranges
    
#------------------------------------------------------------------------------

def plot_anno_fractions(df, fractions):
    '''input is the mos_nob anno df and the fractions from XY.'''
    
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']    
    
    columns = fractions.ANNO.sort_values().unique().tolist()
    
    ranges = [i for i in range(len(columns))]
    
    df = df.copy()[columns + ['MUTSIG_CLASS', 'IS_MOSAIC']]
    
    df['wgEncodeRegTfbsClusteredV3'] = df.wgEncodeRegTfbsClusteredV3.fillna(0)\
                                         .apply(lambda x: x != 0)
    df['wgEncodeRegDnaseClusteredV3'] = df.apply(lambda row:
                                row.wgEncodeRegDnaseClusteredV3 == True,
                                                 axis=1)
    
    df_ = df.groupby('MUTSIG_CLASS').sum()
    
    #is_mosaic was all 1, so effectively replaces a SUM column
    df_ = df_.apply(lambda row: row / row['IS_MOSAIC'], axis=1)
    df_.reset_index(inplace=True)
    
    fractions = fractions[~fractions.PERCENTAGE.duplicated()]
    fractions.sort_values(by='ANNO', inplace=True)
    
    df_ = df_.groupby(['IS_MOSAIC', 'MUTSIG_CLASS']).sum().stack()\
                                                    .reset_index()
    df_.sort_values(by=['level_2', 'MUTSIG_CLASS'], inplace=True)
    
    sns.stripplot(x='level_2', y=0, hue='MUTSIG_CLASS', data=df_,
                  palette=colors, dodge=True)
    
    plt.hlines(y=fractions.PERCENTAGE, xmin=np.subtract(ranges, 0.45),
               xmax=np.add(ranges, 0.45), colors='w')
    plt.hlines(y=fractions.PERCENTAGE, xmin=np.subtract(ranges, 0.45),
               xmax=np.add(ranges, 0.45), colors='k', linestyles='dotted')
    
    plt.xticks(rotation=45)
    
    plt.show()
    

def plot_anno_sims_lymph_semi(df, sims, lymph, semi, category='cc'):
    '''input is the mos_nob table and the assembled sims results. note that for
    exons only the input of the df has to be adjusted accordingly, i.e. the df
    has to be adjusted to include gene region only.'''
    
    sims = sims[sims.CATEGORY_ShufgnoSSC == category]
    sims.drop('CATEGORY_ShufgnoSSC', axis=1, inplace=True)
    sims.set_index('MUTSIG_CLASS', inplace=True)
    
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']    
    
    columns = list(set([col[0] for col in sims.columns]))
    columns.sort()
    
    ranges = make_ranges(columns)
    
    ranges_lines = [i for i in range(len(columns))]
    
    df = df.copy()[columns + ['MUTSIG_CLASS', 'IS_MOSAIC']]
    
    df['wgEncodeRegTfbsClusteredV3'] = df.wgEncodeRegTfbsClusteredV3.fillna(0)\
                                         .apply(lambda x: x != 0)
    df['wgEncodeRegDnaseClusteredV3'] = df.apply(lambda row:
                                row.wgEncodeRegDnaseClusteredV3 == True,
                                                 axis=1)
    
    df_ = df.groupby('MUTSIG_CLASS').sum()
    
    lymph = lymph[lymph.level_1.isin(columns)]
    semi = semi[semi.level_1.isin(columns)]
    
    #is_mosaic was all 1, so effectively replaces a SUM column
    df_ = df_.apply(lambda row: row / row['IS_MOSAIC'], axis=1)
    df_.reset_index(inplace=True)
    
    df_ = df_.groupby(['IS_MOSAIC', 'MUTSIG_CLASS']).sum().stack()\
                                                    .reset_index()
    df_.sort_values(by=['level_2', 'MUTSIG_CLASS'], inplace=True)
    
    sims_ = sims.unstack().reset_index()\
                          .sort_values(by=['level_0', 'MUTSIG_CLASS'])
    means = sims_[sims_.level_1 == 'mean'][0]
    errors = sims.stack().stack().unstack(1).reset_index()\
                         .sort_values(by=['level_1', 'MUTSIG_CLASS'])\
                         [['2.5%', '97.5%']]
    minus = abs(errors['2.5%'].values - means.values)
    plus = abs(means.values - errors['97.5%'].values)
    errors_plotting = [minus, plus]
    
    ecolors = ['0.5' if((value >= error_lst[0]) & (value <= error_lst[1]))
                     else 'r'
               for value, error_lst in zip(df_[0].values, errors.values)]
    
    #actual plotting
    sns.stripplot(x='level_2', y=0, hue='MUTSIG_CLASS', data=df_,
                  palette=colors, dodge=True, jitter=False)
    
    plt.errorbar(x=ranges, y=means, yerr=errors_plotting, linestyle='',
                 marker='_', mfc='0.5', mec='0.5', alpha=0.5, elinewidth=5,
                 ecolor=ecolors)
    
    plt.hlines(y=lymph[0], xmin=np.subtract(ranges_lines, 0.45),
               xmax=np.add(ranges_lines, 0.45), colors='w')
    plt.hlines(y=lymph[0], xmin=np.subtract(ranges_lines, 0.45),
               xmax=np.add(ranges_lines, 0.45), colors='r',
               linestyles='dotted')
    
    plt.hlines(y=semi[0], xmin=np.subtract(ranges_lines, 0.45),
               xmax=np.add(ranges_lines, 0.45), colors='w')
    plt.hlines(y=semi[0], xmin=np.subtract(ranges_lines, 0.45),
               xmax=np.add(ranges_lines, 0.45), colors='b',
               linestyles='dotted')
    
    plt.ylim(0,1.)
    plt.xlabel('')
    plt.ylabel('Fraction of Variants')
    sns.despine(offset=5, trim=True)
    plt.xticks(rotation=45, ha='right')
    
    ax = plt.gca()
    
    ax.get_legend().remove()
    
    plt.show()
    
    
    
    
    
    
    
    
    








