# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:00:58 2020

@author: Martin
"""

#Import modules

import sys, os, gzip
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()

#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions to annotate validation

'''use 20191010_validation1_cc.csv'''

#------------------------------------------------------------------------------
def validation_bool_by_tissue(df):
    '''check for each tissue whether the lower threshold is above both 0.005
    and the upper for JGG; also annotates a FLAG if control lower is above
    0.005.'''
    
    df['SET_BLOOD'] = df.apply(lambda row: valbool(row, 'BLOOD'), axis=1)
    df['SET_SALIVA'] = df.apply(lambda row: valbool(row, 'SALIVA'), axis=1)
    df['SET_SPERM_A'] = df.apply(lambda row: valbool(row, 'SPERM_A'), axis=1)
    df['SET_SPERM_B'] = df.apply(lambda row: valbool(row, 'SPERM_B'), axis=1)
    df['SET_SPERM_C'] = df.apply(lambda row: valbool(row, 'SPERM_C'), axis=1)
    
    df['FLAG_CTRL'] = df.apply(lambda row: row['UPPER_BLOOD_CONTROL'] > 0.005,
                               axis=1)

def valbool(row, tissue):
    '''use tissue information like BLOOD etc. to determine whether the AF
    fulfills the threshold of being called mosaic.'''
    
    if row['LOWER_' + tissue] > 0.005 and\
       row['LOWER_' + tissue] > row['UPPER_BLOOD_CONTROL']:
           return 1
      
    else:
        return 0

#------------------------------------------------------------------------------

def set_mosaic_et_al(df):
    '''use the bools to define mosaic variants.'''
    
    #mosaic in any tissue
    df['SET_MOSAIC'] = ((df.SET_BLOOD + df.SET_SALIVA + df.SET_SPERM_A +
                         df.SET_SPERM_B + df.SET_SPERM_C) > 0).apply(int)
    
    df['SET_SPERM_MOSAIC'] = ((df.SET_SPERM_A == True) |
                              (df.SET_SPERM_B == True) |
                              (df.SET_SPERM_C == True)).apply(int)
    
    df['SET_SOMA_MOSAIC'] = ((df.SET_BLOOD == True) |
                             (df.SET_SALIVA == True)).apply(int)
    
    df['SET_BOTH_MOSAIC'] = ((df.SET_SPERM_MOSAIC == True) &
                              (df.SET_SOMA_MOSAIC == True)).apply(int)
    
    df['SET_SPERM_ONLY'] = ((df.SET_SPERM_MOSAIC == True) &
                            (df.SET_SOMA_MOSAIC == False)).apply(int)
    
    df['SET_SOMA_ONLY'] = ((df.SET_SPERM_MOSAIC == False) &
                          (df.SET_SOMA_MOSAIC == True)).apply(int)
    
    df['FLAG_SALIVA'] = ((df.SET_SALIVA == True) &
                         (df.SET_BLOOD == False)).apply(int)

#------------------------------------------------------------------------------

def replace_INDIVIDUAL_PapID(df):
    '''replace the INDIVIDUAL (here actually SAMPLE) column with paper_IDs and
    save the original IDs in a separate column called SAMPLE_ORIGINAL. Is
    customized to be used for validation data.'''
    
    df['SAMPLE_ORIGINAL'] = df.SAMPLE
    
    dictionary = make_yo_dic(df)
    
    df['SAMPLE'] = df.apply(lambda row: dictionary.get(row.SAMPLE,
                                                       row.SAMPLE),
                                axis=1)
        
    df['SAMPLE'] = df.apply(add_0, axis=1)
    

def make_yo_dic(df):
    
    #not needed for validation table
    #yo = df[df.COHORT.isin(['Young', 'Old'])]
    
    ids = df.SAMPLE.sort_values().unique()
    
    dic = {ID:('S' + str(PapID)) for
               ID, PapID in
                   zip(ids, range(1, (len(ids) + 1)))}
    
    return dic

def add_0(row):
    
    ID = row.SAMPLE
    
    if len(str(ID)) == 2:
        ID = ID[0] + '0' + ID[1]
    
    return ID

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def annotate_samples_time_MAF(df):
    '''annotate samples with what is available.'''

    dic_samples = make_sample_dic()
    
    df['NUMBER_SPERM'] = df.apply(lambda row: dic_samples[row.SAMPLE][0],
                                  axis=1)
    
    df['SALIVA_BOOL']= df.apply(lambda row: dic_samples[row.SAMPLE][1],
                                axis=1)
    
    df['SPERM_ORDER'] = df.apply(lambda row: dic_samples[row.SAMPLE][2],
                                 axis=1)
    
    df['T1'] = df.apply(lambda row: row['MAF_SPERM_' + row['SPERM_ORDER'][0]],
                        axis=1)
    df['T2'] = df.apply(lambda row: row['MAF_SPERM_' + row['SPERM_ORDER'][1]],
                        axis=1)
    df['T3'] = df.apply(lambda row: row['MAF_SPERM_' + row['SPERM_ORDER'][2]],
                        axis=1)
    
    
def make_sample_dic():
    '''use the sample info to make a dictionary in the format:
        
        SAMPLE_ID : (Number of sperm samples, SALIVA_BOOL)'''
        
    dic = {'S01' : (1, 1, 'ABC'), 'S02' : (3, 1, 'BCA'), 'S03' : (3, 1, 'BAC'),
           'S04' : (3, 1, 'CBA'), 'S05' : (2, 1, 'BAC'), 'S06' : (3, 1, 'BAC'),
           'S07' : (3, 1, 'BCA'), 'S08' : (3, 1, 'BCA'), 'S09' : (1, 0, 'ABC'),
           'S10' : (1, 0, 'ABC'), 'S11' : (3, 1, 'BAC'), 'S12' : (3, 1, 'ABC')}
    
    return dic
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def add_max_delta_abs_rel(df):
    '''calculate the max delta and put it into a column, make the relative one
    following that.'''
    
    df['MAX_DELTA'] = df.apply(max_delta, axis=1)
    df['REL_MAX_DELTA'] = df.apply(lambda row: row['MAX_DELTA'] / row['T1'],
                                   axis=1)

def max_delta(row):
    
    if row.NUMBER_SPERM == 1:
        return 0
    
    elif row.NUMBER_SPERM == 2:
        return abs(row.T1 - row.T2)
    
    elif row.NUMBER_SPERM == 3:
        d1 = abs(row.T1 - row.T2)
        d2 = abs(row.T1 - row.T3)
        d3 = abs(row.T2 - row.T3)
        
        return max(d1, d2, d3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plotting

def home_made_paired_plot_t1_3_b(df):
    '''make four aligned pairplots for T1/T2/T3/B.
    
    input version with mosaic variants only (SET_MOSAIC).
    
    this version is using validation data.
    '''
    
    df = df[(df.SET_SPERM_MOSAIC == True)]
    
    df = df.sort_values(by='SAMPLE')
    
    df['T1'] = df.T1**0.5
    df['T2'] = df.T2**0.5
    df['T3'] = df.T3**0.5
    df['B'] = df.MAF_BLOOD**.5
    
    t1_2 = df[df.NUMBER_SPERM > 1]
    #exclude not only 1 and 2, but also this one position due to incomplete
    #data
    t1_3 = df[(df.NUMBER_SPERM == 3) & (df.POS != 63226011)]
    
    f, axs = plt.subplots(nrows=2, ncols=2)
    
    col_1_2 = ['xkcd:brown' if sbm == 1 else 'g'
               for sbm in t1_2.SET_BOTH_MOSAIC]
    col_1_3 = ['xkcd:brown' if sbm == 1 else 'g'
               for sbm in t1_3.SET_BOTH_MOSAIC]
    col = ['xkcd:brown' if sbm == 1 else 'g' for sbm in df.SET_BOTH_MOSAIC]
    
    axs[0,0].scatter(t1_2.T1, t1_2.T2, marker='o', color=col_1_2,
                     edgecolors='w', s=50)
    axs[1,0].scatter(t1_3.T1, t1_3.T3, marker='o', color=col_1_3,
                     edgecolors='w', s=50)
    axs[0,1].scatter(t1_3.T2, t1_3.T3, marker='o', color=col_1_3,
                     edgecolors='w', s=50)
    axs[1,1].scatter(df.T1, df.B, marker='o', color=col, edgecolors='w', s=50)
    
    for i in range(2):
        for j in range(2):
            axs[i,j].set(adjustable='box', aspect='equal')
            axs[i,j].set_xticks([0., 0.2**0.5, 1.])
            axs[i,j].set_xticklabels(['0.0', '0.2', '1.0'])
            axs[i,j].set_yticks([0., 0.2**0.5, 1.])
            axs[i,j].set_yticklabels(['0.0', '0.2', '1.0'])
            sns.despine(offset=5, trim=True, ax=axs[i,j])

    axs[0,0].set_xlabel('t1 Sperm AF (sqr-t)')
    axs[0,0].set_ylabel('t2 Sperm AF (sqr-t)')
    axs[1,0].set_xlabel('t1 Sperm AF (sqr-t)')
    axs[1,0].set_ylabel('t3 Sperm AF (sqr-t)')
    axs[0,1].set_xlabel('t2 Sperm AF (sqr-t)')
    axs[0,1].set_ylabel('t3 Sperm AF (sqr-t)')
    axs[1,1].set_xlabel('t1 Sperm AF (sqr-t)')
    axs[1,1].set_ylabel('Blood AF (sqr-t)')
    
    plt.show()


def home_made_paired_plot_b_s(df):
    '''make four aligned pairplots for b_s.
    
    this version is using validation data.
    '''
    
    df = df[(df.SET_SOMA_MOSAIC == True) & (df.SALIVA_BOOL == True)]
    
    df = df.sort_values(by='SAMPLE')
    
    df['B'] = df.MAF_BLOOD**0.5
    df['S'] = df.MAF_SALIVA**0.5
    
    col = ['xkcd:brown' if sbm == 1 else 'xkcd:bright orange'
           for sbm in df.SET_BOTH_MOSAIC]
    
    f, g = plt.subplots()
    
    plt.scatter(df.B, df.S, color=col, edgecolors='w', s=50)

    g.set_aspect(adjustable='box', aspect='equal')
    g.set_xticks([0., 0.2**0.5, 1.])
    g.set_xticklabels(['0.0', '0.2', '1.0'])
    g.set_yticks([0., 0.2**0.5, 1.])
    g.set_yticklabels(['0.0', '0.2', '1.0'])
    
    plt.xlabel('Blood AF (sqr-t)')
    plt.ylabel('Saliva AF (sqr-t)')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()



def lines_t1_t3(df):
    '''plot all t1_t2_t3 AFs.'''
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM > 1)]
    
    df = df.sort_values(by='SAMPLE')
    
    df['T1'] = df.T1**0.5
    df['T2'] = df.T2**0.5
    df['T3'] = df.T3**0.5

    
    for index, row in df.iterrows():
        
        if row.SET_BOTH_MOSAIC == 1:
            color = 'xkcd:brown'
        else:
            color = 'g'
        
        plt.plot([1,2,3], [row.T1, row.T2, row.T3], color=color, lw=0.2)
    
    plt.show()


def lines_t1_t3_rel(df):
    '''plot lines, but relative to t1'''
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM > 1)]
    
    df = df.sort_values(by='SAMPLE')

    df['T3'] = df.T3/df.T1
    df['T2'] = df.T2/df.T1
    df['T1'] = df.T1/df.T1

    
    for index, row in df.iterrows():
        
        if row.SET_BOTH_MOSAIC == 1:
            color = 'xkcd:brown'
        else:
            color = 'g'
        
        plt.plot([1,2,3], [row.T1, row.T2, row.T3], color=color, lw=0.2)
    
    plt.ylim(0, 2)
    plt.xticks(ticks=[1, 2, 3], labels=['t1', 't2', 't3'])
    plt.yticks(ticks=[0., 0.5, 1., 1.5, 2.])
    plt.xlabel('')
    plt.ylabel('Relative Sperm AF')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()
    

def lines_t1_t3_abs(df):
    '''plot lines, but absolute anchored at t1=0'''
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM > 1)]
    
    df = df.sort_values(by='SAMPLE')

    df['T3'] = df.T3-df.T1
    df['T2'] = df.T2-df.T1
    df['T1'] = df.T1-df.T1

    
    for index, row in df.iterrows():
        
        if row.SET_BOTH_MOSAIC == 1:
            color = 'xkcd:brown'
        else:
            color = 'g'
        
        plt.plot([1,2,3], [row.T1, row.T2, row.T3], color=color, lw=0.2)
    
    plt.xticks(ticks=[1, 2, 3], labels=['t1', 't2', 't3'])
    plt.yticks(ticks=[-0.04, -0.02, 0., .02, .04])
    plt.xlabel('')
    plt.ylabel('Variation of Sperm AF')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()


def max_violins(df):
    '''plot both violins for abs and relative delta max'''
    
    f, axs = plt.subplots(nrows=2, ncols=1)
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM > 1)]
    
    df = df.sort_values(by='SAMPLE')
    
    df['REL_MAX_DELTA'] = 1 + df.REL_MAX_DELTA
    
    col = ['g', 'xkcd:brown']

    
    sns.violinplot(x='SET_MOSAIC', y='MAX_DELTA', data=df, inner=None, cut=0.,
                   ax=axs[0], palette=['0.5'])
    sns.swarmplot(x='SET_MOSAIC', y='MAX_DELTA', data=df, ax=axs[0],
                  hue='SET_BOTH_MOSAIC', palette=col, alpha=0.5, size=4)
    
    sns.violinplot(x='SET_MOSAIC', y='REL_MAX_DELTA', data=df, inner=None,
                   cut=0., ax=axs[1], palette=['0.5'])
    sns.swarmplot(x='SET_MOSAIC', y='REL_MAX_DELTA', data=df, ax=axs[1],
                  hue='SET_BOTH_MOSAIC', palette=col, alpha=0.5, size=4)

    axs[0].set_yticks([-0.04, -0.02, 0., .02, .04])
    axs[0].set_xlabel('')
    axs[0].set_ylabel('')
    axs[0].get_legend().remove()

    axs[1].set_yticks([0., 0.5, 1., 1.5, 2.])
    axs[1].set_xlabel('')
    axs[1].set_ylabel('')
    axs[1].get_legend().remove()
    
    sns.despine(offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])

    plt.show()
    

def clustermap(df, absolute=True):
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM == 3) &
            (df.POS != 63226011)]
    df = df[['T1', 'T2', 'T3']]
    
    if absolute == True:
        
        df['T3'] = df.T3 - df.T1
        df['T2'] = df.T2 - df.T1
        df['T1'] = df.T1 - df.T1
        df = df.append({'T1' : 0., 'T2': 0.04, 'T3' : -0.04},
                       ignore_index=True)
    
    else:
        
        df['T3'] = df.T3 / df.T1
        df['T2'] = df.T2 / df.T1
        df['T1'] = df.T1 / df.T1
        df = df.append({'T1' : 1, 'T2': 0., 'T3' : 2.},
                       ignore_index=True)
    
    df.drop('T1', axis=1, inplace=True)
    
    sns.clustermap(df.T, cmap='vlag')
    
    plt.show()
    
    
def plot_AF_against_delta(df):
    
    f, axs = plt.subplots(nrows=1, ncols=2)
    
    df = df[(df.SET_SPERM_MOSAIC == True) & (df.NUMBER_SPERM > 1)]
    
    df['MAF_SPERM'] = df.apply(lambda row: max(row.T1, row.T2, row.T3),
                               axis=1)
    
    col = ['xkcd:brown' if sbm == 1 else 'xkcd:bright orange'
           for sbm in df.SET_BOTH_MOSAIC]
    
    axs[0].scatter(x=df.MAF_SPERM, y=df.MAX_DELTA, color=col, edgecolors='w',
                   s=50)
    axs[1].scatter(x=df.MAF_SPERM, y=df.REL_MAX_DELTA, color=col,
                   edgecolors='w', s=50)
    
    axs[0].set_xlim(0., 0.25)
    axs[0].set_ylim(0., 0.04)
    axs[0].set_yticks([0., 0.01, 0.02, 0.03, 0.04])
    axs[0].set_xlabel('Max Sperm AF')
    axs[0].set_ylabel('Sperm AF Dmax')

    axs[1].set_xlim(0., 0.25)
    axs[1].set_ylim(0., 1.)
    axs[1].set_xlabel('Max. Sperm AF')
    axs[1].set_ylabel('Relative Sperm AF Dmax')
    
    sns.despine(offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    