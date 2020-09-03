# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 17:24:46 2020

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
#Functions to determine mosaics (outupt table from callers)

'''provided as tab-separated. regular pd.read_table() is used for import.'''

def define_mosaic_firstlevel(df):
    '''uses gnomAD AF, repeat info, segdup info, vicinity to homopolymer/di-
    nucleotide repeat/indel, and a variant's presence in both, Strelka2 and
    MuTect2.'''
    
    df['SET_MOSAIC'] = (
                        (df.UCSC_RPMSK == 'pass') &
                        (df.REPEAT_MASKER == False) &
                        (df.SEGDUP == False) &
                        (df.HOMOPOLYMER == False) &
                        (df.DINUCLEOTIDE == False) &
                        (df.NEAR_INDEL == False) &
                        ~((df.MAF > 0.35) & (df.NORMAL_MAF > 0.35)) &
                        (df.ALT_COUNT > 2) &
                        #called by callers
                        (
                        #either M2 and S2 and gnomAD < 0.01
                        ((df.GNOMAD_FREQ < 0.01) &
                         (df.IN_STRELKA2 == True) & 
                         (df.IN_MUTECT2 == True)
                         ) |
                        #or in MF and not present in gnomAD
                        ((df.GNOMAD_FREQ == 0) &
                         (df.IN_MOSAICFORECAST == True))                       
                        )
                        )
                        

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions to transform output table to anlaysis table
                        
'''provided as comma-separated. regular import. this is a table that assesses
all mosaic variants in all available tissues. pivot can then be done on this
new table.'''

#------------------------------------------------------------------------------
#Make an analysis table with all data organized in rows

def make_analysis_table(df):
    '''master function for making the pivot table.'''
    
    df = df.copy()
    
    fill_validation(df)
    
    harmonize_gene(df)
    
    define_mosaic(df)
    rename_organ(df)
    add_depth(df)
    
    df = drop_unnecessary_colums(df)
    
    df = apply_organ_pivot(df)
    
    return df
    

#base functions
def fill_validation(df):
    '''N/A behaves weirdly for pivot.'''
    
    df.fillna('not_validated', inplace=True)


def harmonize_gene(df):
    '''necessary, as GENE might differ between the same variants by order. will
    use a uniqued table as reference.'''
    
    df_ = df[['CHR_POS_REF_ALT', 'GENE']]
    df_ = df_[~df_.CHR_POS_REF_ALT.duplicated()]
    
    gene_dic = {}
    
    for index, row in df_.iterrows():
        gene_dic[row.CHR_POS_REF_ALT] = row.GENE
    
    df['GENE'] = df.apply(lambda row: gene_dic[row.CHR_POS_REF_ALT], axis=1)


def define_mosaic(df):
    '''use a minimum of three reads to determine mosaicism.'''
    
    df['SET_MOSAIC'] = (df.ALT_COUNT > 2)


def rename_organ(df):
    '''replace - in column organs with _.'''
    
    df['ORGAN'] = df.apply(lambda row: '_'.join(row.ORGAN.upper().split('-')),
                           axis=1)


def drop_unnecessary_colums(df):
    '''drop columns that can be revisited later in master table.'''
    
    df_ = df.drop(['ID', 'REPEAT_MASKER', 'SEGDUP', 'HOMOPOLYMER',
                   'DINUCLEOTIDE', 'NEAR_INDEL', 'UCSC_RPMSK',
                   'NORMAL_REF_COUNT', 'NORMAL_ALT_COUNT', 'NORMAL_MAF',
                   'NORMAL_LOWER_CI', 'NORMAL_UPPER_CI',
                   'NORMAL_CI_IS_GREATER', 'TUMOR_IS_BLOOD_SALIVA',
                   'TUMOR_IS_SPERM', 'INDIVIDUAL-CHR_POS_REF_ALT'],
                  axis=1)
    
    return df_


def add_depth(df):
    '''not included in new version. obtain depth by adding up ref and alt.'''
    
    df['DEPTH'] = df.apply(lambda row: row['REF_COUNT'] + row['ALT_COUNT'],
                           axis=1)


def apply_organ_pivot(df):
    '''use pivot table function to combine all rows into one column for each
    variant.'''
    
    piv = pd.pivot_table(df, values=['DEPTH', 'MAF', 'LOWER_CI', 'UPPER_CI',
                                     'REF_COUNT', 'ALT_COUNT', 'SET_MOSAIC',
                                     'CI_IS_GREATER'],
                             index=['CHR_POS_REF_ALT', 'CHROM', 'POS', 'REF',
                                    'ALT', 'INDIVIDUAL', 'AGE', 'COHORT',
                                    'ANNO', 'GENE', 'GNOMAD_FREQ', 'REF_SEQ',
                                    'VARIANT_GROUP', 'VALIDATION',
                                    'BLOOD_SALIVA_M2_S2_MF',
                                    'SPERM_A_M2_S2_MF', 'SPERM_B_M2_S2_MF',
                                    'SPERM_C_M2_S2_MF'],
                             columns=['ORGAN'])
    
    piv.columns = ['_'.join(col) for col in piv.columns.values]
    
    piv.reset_index(inplace=True)
    
    return piv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up pivot table by flagging

def flag_hets(df):
    '''flag variants not on X or Y that are >0.35 in both blood/saliva and
    sperm A. not necessary for this iteration of the data 20200303.'''
    
    df['FLAG_HET'] = (
                      ((df.CHROM.isin(['X', 'Y'])) &
                       ((df.MAF_BLOOD > 0.7) | (df.MAF_SALIVA > 0.7)) &
                       (df.MAF_SPERM_A > 0.7))
                      |
                      ((~(df.CHROM.isin(['X', 'Y']))) &
                       ((df.MAF_BLOOD > 0.35) | (df.MAF_SALIVA > 0.35)) &
                       (df.MAF_SPERM_A > 0.35))
                      )
    
    
def flag_dups(df):
    '''flag all dups, assuming that they are the result of alignment errors.
    this will get rid of hotspots as well.'''
    
    df['FLAG_DUP'] = (df.CHR_POS_REF_ALT.duplicated(keep=False))
    

#------------------------------------------------------------------------------

def flag_coverage(df):
    '''flag if one of the available categories (i.e. blood, sperm, saliva) is
    outside the 5%-95% range, which was annotated by XY from heterozygous
    variants.'''
    
    df['FLAG_COV'] = df.apply(lambda row: cov_flag_row(row), axis=1)


def cov_flag_row(row):
    '''use information of what samples are available to categorize.'''
    
    if row['INDIVIDUAL'] in [6517, 6762]: #those with B and C
        
        if row.CHROM in ['X', 'Y']:
            
            a = (row.SPA_XY_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_XY_DEPTH_CI95)
            b = (row.SPB_XY_DEPTH_CI5 < row.DEPTH_SPERM_B <
                 row.SPB_XY_DEPTH_CI95)
            c = (row.SPC_XY_DEPTH_CI5 < row.DEPTH_SPERM_C <
                 row.SPC_XY_DEPTH_CI95)
            bl = (row.BL_XY_DEPTH_CI5 < row.DEPTH_BLOOD <
                  row.BL_XY_DEPTH_CI95)
        
        else:
            
            a = (row.SPA_AUTO_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_AUTO_DEPTH_CI95)
            b = (row.SPB_AUTO_DEPTH_CI5 < row.DEPTH_SPERM_B <
                 row.SPB_AUTO_DEPTH_CI95)
            c = (row.SPC_AUTO_DEPTH_CI5 < row.DEPTH_SPERM_C <
                 row.SPC_AUTO_DEPTH_CI95)
            bl = (row.BL_AUTO_DEPTH_CI5 < row.DEPTH_BLOOD <
                  row.BL_AUTO_DEPTH_CI95)
        
        return not (a and b and c and bl)
    
    elif row['INDIVIDUAL'] in [7670, 8090]: #those with saliva

        if row.CHROM in ['X', 'Y']:
            
            a = (row.SPA_XY_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_XY_DEPTH_CI95)
            sl = (row.SL_XY_DEPTH_CI5 < row.DEPTH_SALIVA <
                  row.SL_XY_DEPTH_CI95)
        
        else:
            
            a = (row.SPA_AUTO_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_AUTO_DEPTH_CI95)
            sl = (row.SL_AUTO_DEPTH_CI5 < row.DEPTH_SALIVA <
                  row.SL_AUTO_DEPTH_CI95)
        
        return not (a and sl)
    
    else:
        
        if row.CHROM in ['X', 'Y']:
            
            a = (row.SPA_XY_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_XY_DEPTH_CI95)
            bl = (row.BL_XY_DEPTH_CI5 < row.DEPTH_BLOOD <
                  row.BL_XY_DEPTH_CI95)
        
        else:
            
            a = (row.SPA_AUTO_DEPTH_CI5 < row.DEPTH_SPERM_A <
                 row.SPA_AUTO_DEPTH_CI95)
            bl = (row.BL_AUTO_DEPTH_CI5 < row.DEPTH_BLOOD <
                  row.BL_AUTO_DEPTH_CI95)
        
        return not (a and bl)
        
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotate table

def add_indel_info(df):
    '''add info whether a variant is an indel or a long indel (more than 1bp
    difference).'''
    
    df['INDEL'] = (df.REF.str.len() != df.ALT.str.len())
    df['LONG_INDEL'] = ((df.REF.str.len() > 2) | (df.ALT.str.len() > 2))


def calculate_CI_deltas(df):
    '''delta of lower, upper and the maf.'''
    
    for organ in ['_BLOOD', '_SALIVA', '_SPERM_A', '_SPERM_B', '_SPERM_C']:
    
        df['LOWER_ERROR' + organ] = df['MAF' + organ] - df['LOWER_CI' + organ]
        df['UPPER_ERROR' + organ] = df['UPPER_CI' + organ] - df['MAF' + organ]

#------------------------------------------------------------------------------

def add_mut_cats(df):
    '''extract the mutational category from the sequence snippet.'''
    
    df['TRI_REF'] = df.REF_SEQ.str[7:10]
    
    df['TRI_CAT'] = df.apply(lambda row: categorize_mut(row)[0], axis=1)
    df['CAT'] = df.apply(lambda row: categorize_mut(row)[1], axis=1)


def categorize_mut(row):
    '''use standard approach to make mutationl signature categories.'''
    
    if row['INDEL'] == True:
        
        return ('INDEL', 'INDEL')
    
    elif row['REF'] in set('CT'):
        
        tri_cat = row['TRI_REF']
        cat = row['REF'] + '>' + row['ALT']
        
        return (tri_cat, cat)
    
    else:
        
        tri_cat = rev_comp(row['TRI_REF'])
        cat = rev_comp(row['REF']) + '>' + rev_comp(row['ALT'])
        
        return (tri_cat, cat)

    
def rev_comp(seq):
    '''make reverse complement of sequence. simple, assumes that all are upper
    case and regular nucleotides. left without failsaves to make sure that the
    data passed on is as expected.'''
    reverse_comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_seq = ''
    
    for n in seq[::-1]:
        rev_seq += reverse_comp[n]
    
    return rev_seq

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def replace_INDIVIDUAL_PapID(df):
    '''replace the INDIVIDUAL column with paper_IDs and save the original IDs
    in a separate column called INDIVIDUAL_ORIGINAL. Is customized to be used
    for young/old cohort.'''
    
    df['INDIVIDUAL_ORIGINAL'] = df.INDIVIDUAL
    
    dictionary = make_yo_dic(df)
    
    df['INDIVIDUAL'] = df.apply(lambda row: dictionary.get(row.INDIVIDUAL,
                                                           row.INDIVIDUAL),
                                axis=1)
        
    df['INDIVIDUAL'] = df.apply(add_0, axis=1)
    

def make_yo_dic(df):
    
    yo = df[df.COHORT.isin(['Young', 'Old'])]
    
    ids = yo.INDIVIDUAL.sort_values().unique()
    
    dic = {ID:('ID' + str(PapID)) for
               ID, PapID in
                   zip(ids, range(1, (len(ids) + 1)))}
    
    return dic

def add_0(row):
    
    ID = row.INDIVIDUAL
    
    if len(str(ID)) == 2:
        ID = ID[0] + '0' + ID[1]
    
    return ID

#added for ASD
    
def asd_ids(df):
    '''replace the ID of ASD cases with the NatMed ID.'''
    
    dictionary = {'6058': 'F06', '6091': 'F05', '6099': 'F01', '6107': 'F02',
                  '6314': 'F03', '6451': 'F07', '6463': 'F08', '6490': 'F04'}
    
    df['INDIVIDUAL'] = df.apply(lambda row: dictionary.get(row['INDIVIDUAL'],
                                                           row['INDIVIDUAL']),
                                axis=1)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def correct_t1_3(df):
    '''annotate sperm a, b, c for all with t1, 2, 3, as appropriate.'''
    
    df['T1'] = df.apply(t1, axis=1)
    df['T2'] = df.apply(t2, axis=1)
    df['T3'] = df.apply(t3, axis=1)

def t1(row):
    
    if row['INDIVIDUAL'] in ['ID01', 'ID09', 'ID10', 'ID12']:
        return row['MAF_SPERM_A']
    
    elif row['INDIVIDUAL'] in ['ID02', 'ID03', 'ID05',
                               'ID06', 'ID07', 'ID08', 'ID11']:
        return row['MAF_SPERM_B']
    
    elif row['INDIVIDUAL'] in ['ID04']:
        return row['MAF_SPERM_C']
    
    else:
        return 'zonk'

def t2(row):
    
    if row['INDIVIDUAL'] in ['ID04', 'ID12']:
        return row['MAF_SPERM_B']
    
    elif row['INDIVIDUAL'] in ['ID03', 'ID05', 'ID06', 'ID11']:
        return row['MAF_SPERM_A']
    
    elif row['INDIVIDUAL'] in ['ID02', 'ID07', 'ID08']:
        return row['MAF_SPERM_C']
    
    else:
        return 'zonk'
    
def t3(row):
    
    if row['INDIVIDUAL'] in ['ID03', 'ID06', 'ID11', 'ID12']:
        return row['MAF_SPERM_C']
    
    elif row['INDIVIDUAL'] in ['ID02', 'ID04', 'ID07', 'ID08']:
        return row['MAF_SPERM_A']
    
    else:
        return 'zonk'
    
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define sets

#------------------------------------------------------------------------------
def add_category_alpha(df):
    '''use categorize_row to get category and alpha.'''
    
    df['MOSAIC_CAT'] = df.apply(lambda row: categorize_row(row)[0], axis=1)
    df['ALPHA'] = df.apply(lambda row: categorize_row(row)[1], axis=1)


def categorize_row(row):
    '''categorize each row based on the values present for MAF_SPERM/BLOOD or
    saliva.'''
    
    category = 0
    alpha = 'no_mosaicism_blood_A'
    
    sp_bool = row['ALT_COUNT_SPERM_A'] > 2
    blsal_bool = (row['ALT_COUNT_BLOOD'] > 2 or row['ALT_COUNT_SALIVA'] > 2)
    sal_bool = row['ALT_COUNT_SALIVA'] > 2
        
        
    if sp_bool and not blsal_bool:
        alpha = -1
        category = 1
    elif not sp_bool and blsal_bool:
        alpha = -1
        category = 5
    
    elif sp_bool and blsal_bool:
    
        if sal_bool:
            alpha = row['MAF_SPERM_A']/row['MAF_SALIVA']
        else:
            alpha = row['MAF_SPERM_A']/row['MAF_BLOOD']
        category = 3
            
        if alpha >= 3:
            category = 2
        
        elif alpha <= 1/3:
            category = 4
    
    return category, alpha

#------------------------------------------------------------------------------

def set_define_allmosaics(df):
    '''define sets across all the mosaic variants. use blood and saliva as
    well and call it soma!'''
    
    df['SET_SPERM_MOSAIC'] = ((df.MOSAIC_CAT < 5) &
                              (df.MOSAIC_CAT > 0)).apply(int)
                          
    df['SET_SOMA_MOSAIC'] = (df.MOSAIC_CAT > 1).apply(int)
    
    df['SET_BOTH_MOSAIC'] = ((df.SET_SPERM_MOSAIC == True) &
                             (df.SET_SOMA_MOSAIC == True)).apply(int)
    
    df['SET_SPERM_ENRICHED'] = (df.MOSAIC_CAT == 2).apply(int)
    
    df['SET_SOMA_ENRICHED'] = (df.MOSAIC_CAT == 4).apply(int)
    
    df['SET_EQUALLY_ENRICHED'] = (df.MOSAIC_CAT == 3).apply(int)
    
    df['SET_SPERM_ONLY'] = (df.MOSAIC_CAT == 1).apply(int)
    
    df['SET_SOMA_ONLY'] = (df.MOSAIC_CAT == 5).apply(int)
    
#------------------------------------------------------------------------------
    
def set_mosaic_class(df):
    '''define mosaic class with 1_sperm, 2_shared, 3_soma.'''
    
    df['MOSAIC_CLASS_SSS'] = df.apply(lambda row: extract_mosaic_class(row),
                                      axis=1)

def extract_mosaic_class(row):
    '''actuall function extracting the class.'''
    
    sss = 'zonk'
    
    if row.SET_SPERM_ONLY == True:
        sss = '1_sperm'
    elif row.SET_BOTH_MOSAIC == True:
        sss = '2_shared'
    elif row.SET_SOMA_ONLY == True:
        sss = '3_soma'
    
    return sss

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------    
def make_anno_mutsigs(df):
    '''annotation does two things: 1) function to properly characterize sperm
    calls that are from sperm b and c only, but not in sperm a. essentially
    exchanges 'zonk' and one wrongly assigned blood call. 2) used to assig
    the categories used for mut sig analysis 1_sperm, 2_shared, 3_soma_y,
    4_soma_o.
    
    Note that Blastocyst data has to be removed before analysis!
    Note that Indels have to be removed before analysis!
    
    (probably should be done in function though)
    
    '''
    
    df['MUTSIG_CLASS'] = df.apply(lambda row: row_anno_mutsig(row), axis=1)
    
def row_anno_mutsig(row):
    
    if row['MOSAIC_CLASS_SSS'] == '1_sperm':
        return '1_sperm'
    elif row['MOSAIC_CLASS_SSS'] == 'zonk':
        return '1_sperm'
    
    elif row['MOSAIC_CLASS_SSS'] == '2_shared':
        return '2_shared'
    elif row['POS'] == 35906599:
        return '2_shared'
    
    elif row['MOSAIC_CLASS_SSS'] == '3_soma':
        if row['COHORT'] == 'Young':
            return '3_soma_y'
        else:
            return '4_soma_o'

    else:
        return 'zonk'

#------------------------------------------------------------------------------

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def make_working_table_mos(df):
    '''drop unnecessary columns, remove non-mosaic variants.'''
    
    mos = df[(df.FLAG_DUP == False) & (df.FLAG_HET == False) &
             (df.FLAG_COV == False)]
    
    mos.drop(['FLAG_DUP', 'FLAG_HET', 'FLAG_COV', 'SET_MOSAIC_BLOOD',
              'SET_MOSAIC_SALIVA', 'SET_MOSAIC_SPERM_A', 'SET_MOSAIC_SPERM_B',
              'SET_MOSAIC_SPERM_C', 'BL_AUTO_DEPTH_CI5', 'BL_AUTO_DEPTH_CI95',
              'SL_AUTO_DEPTH_CI5', 'SL_AUTO_DEPTH_CI95', 'SPA_AUTO_DEPTH_CI5',
              'SPA_AUTO_DEPTH_CI95', 'SPB_AUTO_DEPTH_CI5',
              'SPB_AUTO_DEPTH_CI95', 'SPC_AUTO_DEPTH_CI5',
              'SPC_AUTO_DEPTH_CI95', 'BL_XY_DEPTH_CI5', 'BL_XY_DEPTH_CI95',
              'SL_XY_DEPTH_CI5', 'SL_XY_DEPTH_CI95', 'SPA_XY_DEPTH_CI5',
              'SPA_XY_DEPTH_CI95', 'SPB_XY_DEPTH_CI5', 'SPB_XY_DEPTH_CI95',
              'SPC_XY_DEPTH_CI5', 'SPC_XY_DEPTH_CI95'], inplace=True, axis=1)
    
    return mos

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting Functions
    
def plot_with_custum_ci_SP_highnumber(df, category='SPERM_A'):
    
    '''use category as defined by the _XXX. e.g. SPERM_A'''
    
    df = df.sort_values(by=['MAF_' + category], ascending=False)
    df = df.reset_index(drop=True).reset_index()
    plt.plot((df.iloc[:,0] + 1), df['MAF_' + category], color='g')
    plt.plot((df.iloc[:,0] + 1), df['UPPER_CI_' + category], color='0.5')
    plt.plot((df.iloc[:,0] + 1), df['LOWER_CI_' + category], color='0.5')
    plt.xlabel('Ranked Mosaic Variants')
    plt.ylabel('Allelic Fraction (Sperm)')
    #plt.xlim(0,70)
    #plt.ylim(0,0.25)
    sns.despine(offset=5, trim=True)
    plt.show()


def plot_with_custom_ci_sp_sh_bl(df, x_s=300, x_sh=150, x_b=500,
                                 y_s=0.21, y_sh=0.51, y_b=0.21):
    
    df = df.copy()
    
    dict_div = {str(n) : 1 for n in range(1,24)}
    dict_div['X'] = 2
    dict_div['Y'] = 2
    
    for cat in ['MAF_SPERM_A', 'MAF_BLOOD', 'LOWER_ERROR_BLOOD',
                'UPPER_ERROR_BLOOD', 'UPPER_CI_SPERM_A', 'LOWER_CI_SPERM_A',
                'UPPER_CI_BLOOD', 'LOWER_CI_BLOOD']:
        
        df[cat] = df.apply(lambda row: row[cat]/dict_div[str(row['CHROM'])],
                           axis=1)
    
    df_s = df[df.MOSAIC_CLASS_SSS == '1_sperm']\
             .sort_values(by=['MAF_SPERM_A'], ascending=False)
    df_s = df_s.reset_index(drop=True).reset_index()
    
    df_sh = df[df.MOSAIC_CLASS_SSS == '2_shared']\
             .sort_values(by=['MAF_SPERM_A'], ascending=False)
    df_sh = df_sh.reset_index(drop=True).reset_index()
    
    df_b = df[df.MOSAIC_CLASS_SSS == '3_soma']\
             .sort_values(by=['MAF_BLOOD'], ascending=False)
    df_b = df_b.reset_index(drop=True).reset_index()
    
    f, axs = plt.subplots(nrows=3)
    
    axs[0].plot((df_s.iloc[:,0] + 1), df_s['MAF_SPERM_A'], color='g')
    axs[0].plot((df_s.iloc[:,0] + 1), df_s['UPPER_CI_SPERM_A'], color='0.5')
    axs[0].plot((df_s.iloc[:,0] + 1), df_s['LOWER_CI_SPERM_A'], color='0.5')
    
    axs[1].plot((df_sh.iloc[:,0] + 1), df_sh['MAF_SPERM_A'],
                color='xkcd:khaki')
    axs[1].plot((df_sh.iloc[:,0] + 1), df_sh['UPPER_CI_SPERM_A'], color='0.5')
    axs[1].plot((df_sh.iloc[:,0] + 1), df_sh['LOWER_CI_SPERM_A'], color='0.5')
    axs[1].errorbar(x=(df_sh.iloc[:,0] + 1), y=df_sh['MAF_BLOOD'],
               yerr=df_sh[['LOWER_ERROR_BLOOD', 'UPPER_ERROR_BLOOD']].T.values,
               ecolor='xkcd:burnt orange',
               capsize=0, elinewidth=2, marker='o', alpha=0.25,
               linestyle='None', mec='None', mfc='xkcd:black', markersize=4)
    
    axs[2].plot((df_b.iloc[:,0] + 1), df_b['MAF_BLOOD'],
                color='xkcd:bright orange')
    axs[2].plot((df_b.iloc[:,0] + 1), df_b['UPPER_CI_BLOOD'], color='0.5')
    axs[2].plot((df_b.iloc[:,0] + 1), df_b['LOWER_CI_BLOOD'], color='0.5')
    
    axs[0].set_xlim(0,x_s)
    axs[0].set_ylim(0,y_s)
    sns.despine(offset=5, trim=True, ax=axs[0])
    axs[0].set_xlabel('Ranked Mosaic Variants')
    axs[0].set_ylabel('Norm. Sperm AF')
    
    axs[1].set_xlim(0,x_sh)
    axs[1].set_ylim(0,y_sh)
    axs[1].set_yticks([0, 0.25, 0.5])
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_xlabel('Ranked Mosaic Variants')
    axs[1].set_ylabel('Norm. Sperm/Blood AF')
    
    axs[2].set_xlim(0,x_b)
    axs[2].set_ylim(0,y_b)
    sns.despine(offset=5, trim=True, ax=axs[2])
    axs[2].set_xlabel('Ranked Mosaic Variants')
    axs[2].set_ylabel('Norm. Blood AF')
    
    plt.show()

    

def ratio_overlay_combine(df, x=1000, n=20):
    
    f, axs = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3,1]},
                          sharex=True)
    
    df = df.copy()
    
    dict_div = {str(n) : 1 for n in range(1,24)}
    dict_div['X'] = 2
    dict_div['Y'] = 2
    
    for cat in ['MAF_SPERM_A', 'MAF_BLOOD', 'LOWER_ERROR_BLOOD',
                'UPPER_ERROR_BLOOD', 'UPPER_CI_SPERM_A', 'LOWER_CI_SPERM_A']:
        
        df[cat] = df.apply(lambda row: row[cat]/dict_div[str(row['CHROM'])],
                           axis=1)
    
    df['RATIO'] = df.apply(lambda row:
                                math.log((row['MAF_SPERM_A'] + 10**(-8))/
                                         (row['MAF_BLOOD'] + 10**(-8))),
                           axis=1)
    df = df.sort_values(by=['MAF_SPERM_A'], ascending=False)
    df = df.reset_index(drop=True).reset_index()
    
    axs[0].errorbar(x=(df.iloc[:,0] + 1), y=df['MAF_BLOOD'],
                 yerr=df[['LOWER_ERROR_BLOOD', 'UPPER_ERROR_BLOOD']].T.values,
                 ecolor='xkcd:bright orange',
                 capsize=0, elinewidth=2, marker='o', alpha=0.25,
                 linestyle='None', mec='None', mfc='xkcd:black', markersize=4)
    axs[0].plot((df.iloc[:,0] + 1), df['MAF_SPERM_A'], color='g')
    axs[0].plot((df.iloc[:,0] + 1), df['UPPER_CI_SPERM_A'], color='0.5')
    axs[0].plot((df.iloc[:,0] + 1), df['LOWER_CI_SPERM_A'],color='0.5')
    #axs[0].set_ylim(-0.01,0.6)
    axs[0].set_xlim(0,x)
    sns.despine(offset=5, trim=True, ax=axs[0])
    axs[0].set_ylabel('Norm. Sperm/Blood AF')
    
    axs[1].plot((df.iloc[:,0] + 1),df['RATIO'], color='k', linestyle='None',
             marker='o', markersize=3, mfc='None')
    roav = np.convolve(df['RATIO'], np.ones((n,))/n, mode='same')
    axs[1].plot((df.iloc[:,0] + 1), roav, color='0.5')
    axs[1].set_ylim(-5, 21)
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('log(AF Ratio)')
    axs[1].set_xlabel('Ranked Mosaic Variants')
    
    plt.show()


def plot_all_lines(df):
    '''plot to visualize the AF for all 4 categories in two plots (sperm af and
    blood af).'''
    
    f, axs = plt.subplots(ncols=2, sharey=True)
    
    df = df.copy()
    
    dict_div = {str(n) : 1 for n in range(1,24)}
    dict_div['X'] = 2
    dict_div['Y'] = 2
    
    for cat in ['MAF_SPERM_A', 'MAF_BLOOD', 'LOWER_ERROR_BLOOD',
                'UPPER_ERROR_BLOOD', 'UPPER_CI_SPERM_A', 'LOWER_CI_SPERM_A']:
        
        df[cat] = df.apply(lambda row: row[cat]/dict_div[str(row['CHROM'])],
                           axis=1)
    
    df_1 = df[df.MUTSIG_CLASS == '1_sperm']
    df_1 = df_1.sort_values(by='MAF_SPERM_A', ascending=False)
    df_1 = df_1.reset_index(drop=True).reset_index()
    
    df_2 = df[df.MUTSIG_CLASS == '2_shared']
    df_2_s = df_2.sort_values(by='MAF_SPERM_A', ascending=False)
    df_2_s = df_2_s.reset_index(drop=True).reset_index()
    df_2_b = df_2.sort_values(by='MAF_BLOOD', ascending=False)
    df_2_b = df_2_b.reset_index(drop=True).reset_index()
    
    df_3 = df[df.MUTSIG_CLASS == '3_soma_y']
    df_3 = df_3.sort_values(by='MAF_BLOOD', ascending=False)
    df_3 = df_3.reset_index(drop=True).reset_index()
    
    df_4 = df[df.MUTSIG_CLASS == '4_soma_o']
    df_4 = df_4.sort_values(by='MAF_BLOOD', ascending=False)
    df_4 = df_4.reset_index(drop=True).reset_index()

    axs[0].plot((df_1.iloc[:,0] + 1), df_1.MAF_SPERM_A, color='g')
    axs[0].plot((df_2_s.iloc[:,0] + 1), df_2_s.MAF_SPERM_A, color='xkcd:brown')
    
    axs[1].plot((df_2_b.iloc[:,0] + 1), df_2_b.MAF_BLOOD, color='xkcd:brown')
    axs[1].plot((df_3.iloc[:,0] + 1), df_3.MAF_BLOOD, color='xkcd:golden')
    axs[1].plot((df_4.iloc[:,0] + 1), df_4.MAF_BLOOD,
                color='xkcd:orangish red')
    
    axs[0].set_xlabel('Ranked Mosaic Variants')
    axs[0].set_ylabel('Norm. Sperm AF')
    #axs[0].set_xlim(0,550)
    axs[0].set_ylim(0,0.61)
    axs[0].set_xticks([0, 100, 200, 300, 400, 500])
   
    axs[1].set_xlabel('Ranked Mosaic Variants')
    axs[1].set_ylabel('Norm. Blood AF')
    #axs[0].set_xlim(0,550)
    axs[0].set_ylim(0,0.61)
    axs[1].set_xticks([0, 400, 800, 1200, 1600])

    sns.despine(offset=5, trim=True,ax=axs[0])    
    sns.despine(offset=5, trim=True,ax=axs[1])
    
    plt.show()
    

def mosaic_number_per_individual(df, y=500):
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['INDIVIDUAL', 'MOSAIC_CLASS_SSS'])
    
    sns.countplot(x='INDIVIDUAL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors)
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')

    plt.show()


def mosaic_number_per_individual_broken_axis(df):
    
    colors = ['g', 'xkcd:brown', 'xkcd:bright orange']
    
    df = df.sort_values(by=['INDIVIDUAL', 'MOSAIC_CLASS_SSS'])
    
    f, axs = plt.subplots(nrows=2, ncols=1,
                          gridspec_kw = {'height_ratios':[1,3]})
    
    sns.countplot(x='INDIVIDUAL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[0])
    sns.countplot(x='INDIVIDUAL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[1])
    
    axs[0].set_ylim(60,500)
    axs[0].set_yscale('log')
    axs[1].set_ylim(0,60)
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('Number of Mosaic Variants')
    axs[1].set_xlabel('')
    axs[0].set_xlabel('')
    axs[0].set_ylabel('')
    axs[0].set_xticks([])

    plt.show()
    
    
def mosaic_number_per_individual_broken_axis_asd(df):
    
    colors = ['g', 'xkcd:brown', 'xkcd:bright orange']
    
    df = df.sort_values(by=['INDIVIDUAL', 'MOSAIC_CLASS_SSS'])
    
    f, axs = plt.subplots(nrows=2, ncols=1,
                          gridspec_kw = {'height_ratios':[1,3]})
    
    sns.countplot(x='INDIVIDUAL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[0])
    sns.countplot(x='INDIVIDUAL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[1])
    
    axs[0].set_ylim(60,500)
    axs[1].set_ylim(0,60)
    axs[0].set_yticks([60, 250, 500])
    axs[1].set_yticks([0, 20, 40, 60])
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('Number of Mosaic Variants')
    axs[1].set_xlabel('')
    axs[0].set_xlabel('')
    axs[0].set_ylabel('')
    axs[0].set_xticks([])
    
    axs[0].get_legend().remove()
    axs[1].get_legend().remove()

    plt.show()

def mosaic_number_per_cohort(df, y=500):
    '''used for cohort-specific plot by only feeding one cohort at a time.'''
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['COHORT', 'MOSAIC_CLASS_SSS'])
    
    sns.countplot(x='COHORT', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors)
    
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')

    plt.show()


def mosaic_number_indel_onecohort(df, y=500, cohort='Young'):
    '''used for cohort-specific plot by keyword above.
    this version has a handler for the cohort and only prints INDELs
    separately.'''
    
    df = df[df.COHORT == cohort]
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'INDEL'])
    
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors)
    
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')

    plt.show()


def mosaic_number_indel_onecohort_oldbrokenaxis(df):
    '''has to be used b/c of representation issues.'''
    
    df = df[df.COHORT == 'Old']
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'INDEL'])
    
    f, axs = plt.subplots(nrows=2, ncols=1,
                          gridspec_kw = {'height_ratios':[1,3]})
    
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[0])
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[1])
    
    axs[0].set_ylim(150,1200)
    axs[1].set_ylim(0,150)
    axs[0].set_yticks([150, 600, 1200])
    axs[1].set_yticks([0, 50, 100, 150])
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('Number of Mosaic Variants')
    axs[1].set_xlabel('')
    axs[0].set_xlabel('')
    axs[0].set_ylabel('')
    axs[0].set_xticks([])
    axs[0].get_legend().remove()
    axs[1].get_legend().remove()
    
    plt.show()
    
    
def mosaic_number_indel_onecohort_ASDbrokenaxis(df):
    '''has to be used b/c of representation issues. actually returns and svg
    that cannot be placed in Illustrator, b/c there seems to be an issue with
    0 value bars. Was fixed by adding a fake 2_shared that is deleted in post-
    processing. not ideal, but could not find another work around.'''
    
    df = df[df.COHORT == 'ASD']
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'INDEL'])
    
    
    
    f, axs = plt.subplots(nrows=2, ncols=1,
                          gridspec_kw = {'height_ratios':[1,3]})
    
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[0])
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors, ax=axs[1])
    
    axs[0].set_ylim(75,600)
    axs[1].set_ylim(0,75)
    axs[0].set_yticks([75, 300, 600])
    axs[1].set_yticks([0, 25, 50, 75])
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('Number of Mosaic Variants')
    axs[1].set_xlabel('')
    axs[0].set_xlabel('')
    axs[0].set_ylabel('')
    axs[0].set_xticks([])
    axs[0].get_legend().remove()
    axs[1].get_legend().remove()
    
    plt.show()
    
def blood_and_sperm_violins(df, spermA_flag=True):
    
    ''' plot two violin plots for sperm af and blood af for both sperm/shared
    and blood shared. use spermA_flag to exclude variants that are only called
    in b/c to not bias towards the two multi-sample samples.
    
    Note that this is currently optimized for sperm and blood, will need slight
    changes for saliva.'''
    
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'COHORT'], ascending=False)
    
    if spermA_flag == True:
        df = df[df.VARIANT_GROUP == 'A']
    
    sperm = df[df.MOSAIC_CLASS_SSS != '3_soma']
    blood = df[df.MOSAIC_CLASS_SSS != '1_sperm']
    
    f, axs = plt.subplots(nrows=1, ncols=2)
    
    sns.violinplot(x='MOSAIC_CLASS_SSS', y='MAF_SPERM_A', data=sperm,
                   hue='COHORT', ax=axs[0], cut=0)
    
    sns.violinplot(x='MOSAIC_CLASS_SSS', y='MAF_BLOOD', data=blood,
                   hue='COHORT', ax=axs[1], cut=0)
    
    axs[0].set_ylabel('Allelic Fraction (Sperm)')
    axs[1].set_ylabel('Allelic Fraction (Blood)')
    
    axs[0].set_ylim(0, 0.51)
    axs[1].set_ylim(0, 0.51)
    
    sns.despine(offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    
    axs[0].get_legend().remove()
    axs[1].legend(edgecolor='None')
    
    plt.show()
    
    
def blood_sperm_comparison_box(df, spermA_flag=True):
    
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'COHORT'], ascending=False)
    
    if spermA_flag == True:
        df = df[df.VARIANT_GROUP == 'A']
    
    sns.boxplot(x='MOSAIC_CLASS_SSS', y='POS',
                data=df.groupby(['INDIVIDUAL', 'COHORT', 'MOSAIC_CLASS_SSS'])\
                                                        .count().reset_index(),
                hue='COHORT')

    sns.swarmplot(x='MOSAIC_CLASS_SSS', y='POS',
                  data=df.groupby(['INDIVIDUAL', 'COHORT', 'MOSAIC_CLASS_SSS'])\
                                                        .count().reset_index(),
                  hue='COHORT', color='w', edgecolor='k', linewidth=1,
                  alpha=0.5)
    
    sns.despine(offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')
    plt.xlabel('')

    plt.show()


def mini_violins(df):
    '''make af violins for all individuals for sperm, shared-sperm,
    shared-blood, blood.'''
    
    df = df.copy()
    
    #make squre root-transformed column for use in violin
    df['MAF_SPERM_A_'] = df.MAF_SPERM_A**0.5
    df['MAF_BLOOD_'] = df.MAF_BLOOD**0.5
    
    ids = df.INDIVIDUAL.sort_values().unique()
    n = len(ids)
    
    f, axs = plt.subplots(nrows=4, ncols=len(ids))
    plt.subplots_adjust(wspace=0.)
    
    for i, sid in enumerate(ids):
        
        #sperm-specific
        sns.violinplot(x='INDIVIDUAL', y='MAF_SPERM_A_',
                       data=df[(df.INDIVIDUAL == sid) &
                               (df.MOSAIC_CLASS_SSS == '1_sperm')],
                       palette=['g'], ax=axs[0,i], cut=0., inner=None)
        
        #shared sp-af
        sns.violinplot(x='INDIVIDUAL', y='MAF_SPERM_A_',
                       data=df[(df.INDIVIDUAL == sid) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:khaki'], ax=axs[1,i], cut=0., inner=None)
    
        #shared bl-af
        sns.violinplot(x='INDIVIDUAL', y='MAF_BLOOD_',
                       data=df[(df.INDIVIDUAL == sid) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:burnt orange'], ax=axs[2,i], cut=0.,
                       inner=None)
        
        #blood-specific
        sns.violinplot(x='INDIVIDUAL', y='MAF_BLOOD_',
                       data=df[(df.INDIVIDUAL == sid) &
                               (df.MOSAIC_CLASS_SSS == '3_soma')],
                       palette=['xkcd:bright orange'], ax=axs[3,i], cut=0.,
                       inner=None)
    
    for i in range(4):
        for j in range(n):
        
            axs[i,j].set_ylim(0.,1.)
            axs[i,j].set_xlabel('')
            axs[i,j].set_ylabel('')
            sns.despine(bottom=True, trim=True, offset=5, ax=axs[i,j])
    
    for i in range(3):
        for j in range(n):
        
            axs[i,j].set_xticks([])
    
    for i in range(4):
        for j in range(1, n):
        
            sns.despine(left=True, trim=True, ax=axs[i,j])
            axs[i,j].set_yticks([])
    
    for i in range(4):
        
        axs[i,0].set_yticks([0., 0.2**0.5, 1.])
        axs[i,0].set_yticklabels(['0.0', '0.2', '1.0'])
    
    axs[2,0].set_ylabel('Sperm/Blood AF (square root-transformed)')
    
    plt.show()

    
def mini_violins_cohort(df):
    '''make af violins for entire cohort'''
    
    df = df.copy()
    
    #make squre root-transformed column for use in violin
    df['MAF_SPERM_A_'] = df.MAF_SPERM_A**0.5
    df['MAF_BLOOD_'] = df.MAF_BLOOD**0.5
    
    cohorts = df.COHORT.sort_values().unique()
    n = len(cohorts)
    
    f, axs = plt.subplots(nrows=4, ncols=len(cohorts))
    plt.subplots_adjust(wspace=0.)
    
    for i, cohort in enumerate(cohorts):
        
        #sperm-specific
        sns.violinplot(x='COHORT', y='MAF_SPERM_A_',
                       data=df[(df.COHORT == cohort) &
                               (df.MOSAIC_CLASS_SSS == '1_sperm')],
                       palette=['g'], ax=axs[0,i], cut=0., inner=None)
        
        #shared sp-af
        sns.violinplot(x='COHORT', y='MAF_SPERM_A_',
                       data=df[(df.COHORT == cohort) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:khaki'], ax=axs[1,i], cut=0., inner=None)
    
        #shared bl-af
        sns.violinplot(x='COHORT', y='MAF_BLOOD_',
                       data=df[(df.COHORT == cohort) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:burnt orange'], ax=axs[2,i], cut=0.,
                       inner=None)
        
        #blood-specific
        sns.violinplot(x='COHORT', y='MAF_BLOOD_',
                       data=df[(df.COHORT == cohort) &
                               (df.MOSAIC_CLASS_SSS == '3_soma')],
                       palette=['xkcd:bright orange'], ax=axs[3,i], cut=0.,
                       inner=None)
    
    for i in range(4):
        for j in range(n):
        
            axs[i,j].set_ylim(0.,1.)
            axs[i,j].set_xlabel('')
            axs[i,j].set_ylabel('')
            sns.despine(bottom=True, trim=True, offset=5, ax=axs[i,j])
    
    for i in range(3):
        for j in range(n):
        
            axs[i,j].set_xticks([])
    
    for i in range(4):
        for j in range(1, n):
        
            sns.despine(left=True, trim=True, ax=axs[i,j])
            axs[i,j].set_yticks([])
    
    for i in range(4):
        
        axs[i,0].set_yticks([0., 0.2**0.5, 1.])
        axs[i,0].set_yticklabels(['0.0', '0.2', '1.0'])
    
    axs[2,0].set_ylabel('Sperm/Blood AF (square root-transformed)')
    
    plt.show()


def home_made_paired_plot_t1_3_b(df):
    '''make four  aligned pairplots for T1/T2/T3/B. has to have a list of
    sperm-present variants as input. use table that contains all variants that
    are mosaic and filters are set for specific data set within function.'''
    
    df = df[df.INDIVIDUAL.isin(['ID04', 'ID12'])] #only three WGS individuals
    df = df[(df.SET_SPERM_MOSAIC == True) | (df.VARIANT_GROUP == 'BC')]
    
    df = df.sort_values(by='INDIVIDUAL')
    
    df['T1'] = df.T1**0.5
    df['T2'] = df.T2**0.5
    df['T3'] = df.T3**0.5
    df['B'] = df.MAF_BLOOD**.5
    
    f, axs = plt.subplots(nrows=2, ncols=2)
    
    colors = ['xkcd:brown' if sss == '2_shared' else 'g'
             for sss in df.MOSAIC_CLASS_SSS]
    
    axs[0,0].scatter(df.T1, df.T2, marker='o', color=colors, edgecolors='w',
                     s=50)
    axs[1,0].scatter(df.T1, df.T3, marker='o', color=colors, edgecolors='w',
                     s=50)
    axs[0,1].scatter(df.T2, df.T3, marker='o', color=colors, edgecolors='w',
                     s=50)
    axs[1,1].scatter(df.T1, df.B, marker='o', color=colors, edgecolors='w',
                     s=50)
    
    for i in range(2):
        for j in range(2):
            axs[i,j].set(adjustable='box-forced', aspect='equal')
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


def age_counts(df):
    '''use mosaic variants of young and old only as input.'''
    
    counts = df.groupby(['INDIVIDUAL', 'MOSAIC_CLASS_SSS']).count()\
                                                           .reset_index()
    
    counts = counts[['INDIVIDUAL', 'MOSAIC_CLASS_SSS', 'POS']]
    
    counts['AGE'] = df.groupby(['INDIVIDUAL', 'MOSAIC_CLASS_SSS']).mean()\
                                                        .reset_index()['AGE']
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    
    sns.lmplot(x='AGE', y='POS', data=counts, hue='MOSAIC_CLASS_SSS',
               palette=colors, markers=['o', 's', '^'])
    
    plt.xlim(0,70)
    plt.ylim(-10, 500)
    plt.xlabel('Age [years]')
    plt.ylabel('Number of Mosaic Variants')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()

def age_counts_boxes(df):
    '''make thin boxplot to add to figure.'''
    
    counts = df.groupby(['COHORT', 'INDIVIDUAL', 'MOSAIC_CLASS_SSS']).count()\
                                                           .reset_index()
    
    counts = counts[['COHORT', 'INDIVIDUAL', 'MOSAIC_CLASS_SSS', 'POS']]
    
    counts.sort_values(by=['MOSAIC_CLASS_SSS', 'INDIVIDUAL'], inplace=True)
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    
    sns.boxplot(x='COHORT', y='POS', hue='MOSAIC_CLASS_SSS', data=counts,
                palette=colors)
    
    plt.ylim(-10, 500)
    plt.ylabel('')
    
    sns.despine(offset=5, trim=True, bottom=True)
    
    plt.show()
    

def age_afs(df):
    '''use mosaic variants of young and old only as input.'''
    
    afs = df.groupby(['INDIVIDUAL', 'MOSAIC_CLASS_SSS']).mean()\
                                                           .reset_index()
    
    afs = afs[['INDIVIDUAL', 'MOSAIC_CLASS_SSS', 'MAF_SPERM_A', 'MAF_BLOOD']]
    
    afs['AGE'] = df.groupby(['INDIVIDUAL', 'MOSAIC_CLASS_SSS']).mean()\
                                                        .reset_index()['AGE']
                                                        
    afs['MAF_SPERM_A'] = afs.MAF_SPERM_A**0.5
    afs['MAF_BLOOD'] = afs.MAF_BLOOD**0.5
    
    #colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    
    sns.regplot(x='AGE', y='MAF_SPERM_A',
               data=afs[afs.MOSAIC_CLASS_SSS == '1_sperm'],
               color='g', marker='o')
    
    sns.regplot(x='AGE', y='MAF_SPERM_A',
               data=afs[afs.MOSAIC_CLASS_SSS == '2_shared'],
               color='xkcd:khaki', marker='+')
    
    sns.regplot(x='AGE', y='MAF_BLOOD',
               data=afs[afs.MOSAIC_CLASS_SSS == '2_shared'],
               color='xkcd:burnt orange', marker='x')
    
    sns.regplot(x='AGE', y='MAF_BLOOD',
               data=afs[afs.MOSAIC_CLASS_SSS == '3_soma'],
               color='xkcd:bright orange', marker='^')
    
    plt.xlim(0,70)
    plt.ylim(-0.005, 0.2**0.5)
    
    plt.yticks(ticks=[0.,0.1**0.5, 0.2**0.5],
               labels=['0.0', '0.1', '0.2'])
    
    plt.xlabel('Age [years]')
    plt.ylabel('Sperm/Blood AF')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()
    

def age_afs_individual_values(df):
    '''make individual value plots for each of the four categories.'''
    
    afs = df.copy()
    afs['MAF_SPERM_A'] = afs.MAF_SPERM_A**0.5
    afs['MAF_BLOOD'] = afs.MAF_BLOOD**0.5
    
    fig, axs = plt.subplots(nrows=4, ncols=1)
    
    sns.regplot(x='AGE', y='MAF_SPERM_A',
               data=afs[afs.MOSAIC_CLASS_SSS == '1_sperm'],
               color='g', marker='o', ax=axs[0])
    
    sns.regplot(x='AGE', y='MAF_SPERM_A',
               data=afs[afs.MOSAIC_CLASS_SSS == '2_shared'],
               color='xkcd:khaki', marker='+', ax=axs[1])
    
    sns.regplot(x='AGE', y='MAF_BLOOD',
               data=afs[afs.MOSAIC_CLASS_SSS == '2_shared'],
               color='xkcd:burnt orange', marker='x', ax=axs[2])
    
    sns.regplot(x='AGE', y='MAF_BLOOD',
               data=afs[afs.MOSAIC_CLASS_SSS == '3_soma'],
               color='xkcd:bright orange', marker='^', ax=axs[3])
    
    for ax in axs:
        
        ax.set_xlim(10, 70)
        ax.set_ylim(0., 1.)
        ax.set_yticks([0.,0.2**0.5, 1.])
        ax.set_yticklabels(['0.0', '0.2', '1.0'])
        ax.set_xlabel('')
        ax.set_ylabel('')
        sns.despine(offset=5, trim=True, ax=ax)
        
    axs[3].set_xlabel('Age [yrs]')
    axs[0].set_ylabel('Sperm AF (sqrt-t)')
    axs[1].set_ylabel('Sperm AF (sqrt-t)')
    axs[2].set_ylabel('Blood AF (sqrt-t)')
    axs[3].set_ylabel('Blood AF (sqrt-t)')
    
    plt.show()


def all_regs(df):
    
    df = df.sort_values(by='INDIVIDUAL')
    
    df['MAF_SPERM_A'] = df.MAF_SPERM_A**0.5
    df['MAF_BLOOD'] = df.MAF_BLOOD**0.5

    inds = df.INDIVIDUAL.unique()
    
    f, axs = plt.subplots(nrows=3, ncols=6, sharex=True, sharey=True)
    
    for i, ind in enumerate(inds):
        
        sns.regplot(x='MAF_SPERM_A', y='MAF_BLOOD',
                    data=df[(df.MOSAIC_CLASS_SSS == '2_shared') &
                            (df.INDIVIDUAL == ind)],
                    ax=axs[i//6, i-(i//6)*6])
        axs[i//6, i-(i//6)*6].set_title(ind)
        axs[i//6, i-(i//6)*6].set_xlim(0,1.)
        axs[i//6, i-(i//6)*6].set_ylim(0,1.)
        
        axs[i//6, i-(i//6)*6].set_xticks([0., 0.2**0.5, 1.])
        axs[i//6, i-(i//6)*6].set_xticklabels(['0.0', '0.2', '1.0'])
        axs[i//6, i-(i//6)*6].set_yticks([0., 0.2**0.5, 1.])
        axs[i//6, i-(i//6)*6].set_yticklabels(['0.0', '0.2', '1.0'])
    
    plt.show()
        
        
def plot_kde_hist(df):
    '''distplots for yo split by young, old-no ch, old-ch; optimized for data
    etc. needs the young/old mosaics as input.'''
    
    f, axs = plt.subplots(ncols=2, sharey=True)
    
    
    sns.distplot(df[(df.COHORT == 'Young') &
                    (df.MOSAIC_CLASS_SSS == '3_soma')].MAF_BLOOD**0.5,
                 bins=[(n*0.01)**0.5 for n in range(21)], color='xkcd:salmon',
                 kde=False, norm_hist=True, ax=axs[0])
    
    sns.distplot(df[(df.COHORT == 'Young') &
                    (df.MOSAIC_CLASS_SSS == '3_soma')].MAF_BLOOD**0.5,
                 bins=[(n*0.01)**0.5 for n in range(21)], color='xkcd:salmon',
                 kde=False, norm_hist=True, ax=axs[1])
    
    sns.distplot(df[(df.INDIVIDUAL.isin(['ID13', 'ID15', 'ID16'])) &
                    (df.MOSAIC_CLASS_SSS == '3_soma')].MAF_BLOOD**0.5,
                 bins=[(n*0.01)**0.5 for n in range(21)],
                 color='xkcd:cerulean', kde=False, norm_hist=True, ax=axs[0])
    
    sns.distplot(df[(df.INDIVIDUAL.isin(['ID14', 'ID17'])) &
                    (df.MOSAIC_CLASS_SSS == '3_soma')].MAF_BLOOD**0.5,
                 bins=[(n*0.01)**0.5 for n in range(21)],
                 color='xkcd:jade', kde=False, norm_hist=True, ax=axs[1])
    
    axs[0].set_xlim(0, 0.2**0.5)
    axs[1].set_xlim(0, 0.2**0.5)

    axs[0].set_xticks([0., 0.05**0.5, 0.10**0.5, 0.15**0.5, 0.20**0.5])
    axs[1].set_xticks([0., 0.05**0.5, 0.10**0.5, 0.15**0.5, 0.20**0.5])

    axs[0].set_xticklabels(['0.0', '0.05', '0.10', '0.15', '0.20'])
    axs[1].set_xticklabels(['0.0', '0.05', '0.10', '0.15', '0.20'])
    
    axs[0].set_ylabel('Normalized Density')
    axs[0].set_xlabel('')
    axs[1].set_xlabel('Blood AF (sqrt-t)')
    
    sns.despine(offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    
    plt.show()
    

def plot_bars_kde_with95pop(df):
    
    df = df[df.MUTSIG_CLASS.isin(['1_sperm', '2_shared'])]
    
    df['MAF_SPERM'] = df.apply(lambda row: max(row.MAF_SPERM_A,
                                               row.MAF_SPERM_B,
                                               row.MAF_SPERM_C),
                               axis=1)
    
    #needed for kdeplot
    df_ = df.copy()
    
    #barplot-preparation
    df['BINS'] = pd.cut(df.MAF_SPERM, bins=[0., 0.05, 0.10, 0.15,
                                            0.20, 0.25, 1.])
    
    df = df.groupby(['BINS', 'MUTSIG_CLASS']).count()
    df['SUM'] = df['POS'].sum()
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    df.sort_values(by=['BINS', 'MUTSIG_CLASS'], inplace=True)
    df.REL_COUNT.fillna(0, inplace=True)
    
    pos = [0, 1, 2, 3, 4, 5]
    width = 0.9
    names = ['0.00-0.05', '0.05-0.10', '0.10-0.15', '0.15-0.20', '0.20-0.25',
             '0.25+']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0,0,0,0,0,0]
    
    for msc in df.MUTSIG_CLASS.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.MUTSIG_CLASS == msc].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.MUTSIG_CLASS == msc].REL_COUNT)
        
    colors = ['g', 'xkcd:brown']
    
    #kde with pop-interval preparation
    desc = df_.MAF_SPERM.describe(percentiles=[0.025, 0.975])
    perc = desc['2.5%'], desc['97.5%']
    
    f, axs = plt.subplots(ncols=2)
    
    for i in range(2):
        
        axs[0].bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                   edgecolor='None', width=width)
    
    sns.kdeplot(df_.MAF_SPERM, shade=True, color='0.8', ax=axs[1])
    axs[1].vlines(x=perc, ymin=0, ymax=28, color='r', linestyle='--')
    
    axs[0].set_xlabel('Sperm AF Bin')
    axs[0].set_ylabel('Relative Contribution')
    
    axs[1].set_xlabel('Sperm AF')
    axs[1].set_ylabel('Normalized Density')
    axs[1].get_legend().remove()
    
    axs[0].set_ylim(0,1)
    axs[1].set_ylim(-1, 30)
    axs[1].set_xlim(0, 0.61)
    
    sns.despine(bottom=True, offset=5, trim=True, ax=axs[0])
    sns.despine(offset=5, trim=True, ax=axs[1])
    
    axs[0].set_xticks(pos)
    axs[0].set_xticklabels(names, rotation=45)
    plt.show()


#------------------------------------------------------------------------------
def scatter_blood_sperm(df, four_classes=False):
    '''sperm vs blood scatter, scaled with sqrt-t, so no additional arguments
    are needed for it. color different points by MUTSIG_CLASS.'''
    
    df = df.copy()
    
    df['MAF_SPERM'] = df.apply(lambda row: max(row.MAF_SPERM_A,
                                               row.MAF_SPERM_B,
                                               row.MAF_SPERM_C),
                               axis=1)
            
    df['MAF_SPERM'] = df.MAF_SPERM**0.5
    df['MAF_BLOOD'] = df.MAF_BLOOD**0.5
    
    colors = assign_color(df, four_classes=four_classes)
    
    f, g = plt.subplots()
    
    plt.scatter(df.MAF_SPERM, df.MAF_BLOOD, color=colors, edgecolors='w', s=50)

    g.set_aspect(adjustable='box-forced', aspect='equal')
    g.set_xticks([0., 0.2**0.5, 1.])
    g.set_xticklabels(['0.0', '0.2', '1.0'])
    g.set_yticks([0., 0.2**0.5, 1.])
    g.set_yticklabels(['0.0', '0.2', '1.0'])
    
    plt.xlabel('Sperm AF (sqr-t)')
    plt.ylabel('Blood AF (sqr-t)')
    
    sns.despine(offset=5, trim=True)
    
    plt.show()
    
def assign_color(df, four_classes=False):
    
    colors = []
    
    for msc in df.MUTSIG_CLASS:
        
        if msc == '1_sperm':
            colors.append('g')
        elif msc == '2_shared':
            colors.append('xkcd:brown')
        elif msc == '3_soma_y':
            if four_classes == False:
                colors.append('xkcd:bright orange')
            else:
                colors.append('xkcd:golden')
        elif msc == '4_soma_o':
            if four_classes == False:
                colors.append('xkcd:bright orange')
            else:
                colors.append('xkcd:orangish red')
        else:
            colors.append('zonk')
    
    return colors
#------------------------------------------------------------------------------


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#------------------------------------------------------------------------------
#mut sig functions

#make an appropriate table
def make_0_cat_tricat():
    '''makes all 96 combinations with a 0 at the POS column to ensure that all
    categories are used.'''
    
    refs = ['C', 'T']
    nts = ['A', 'T', 'G', 'C']
    
    cats = [ref + '>' + alt for ref in refs for alt in nts]
    cats.remove('T>T')
    cats.remove('C>C')
    
    pairs = []
    
    for cat in cats:
        for N1 in nts:
            for N2 in nts:
                tri = N1 + cat[0] + N2
                pairs.append((cat, tri))
    
    df = pd.DataFrame(data=pairs, columns=['CAT', 'TRI_CAT'])
    df['POS'] = 0
    
    return df

def make_table_TRIs_CATs(df, divide_category=False, cat_name=None):
    '''wrangles the df to summarize and prepare for normalization in the
    plotting function. categories have been assigned before e.g. MOSAIC_CLASS_
    SSS.'''
    
    if divide_category == True:
        cat = [cat_name, 'CAT']
        tri = [cat_name, 'CAT', 'TRI_CAT']
        df_pre = df[df.INDEL == False][['INDIVIDUAL', 'CHROM', 'POS',
                                        'TRI_CAT','CAT', cat_name]]
    else:
        cat = 'CAT'
        tri = ['CAT', 'TRI_CAT']
        df_pre = df[df.INDEL == False][['INDIVIDUAL', 'CHROM', 'POS',
                                        'TRI_CAT', 'CAT']]
    
    df_tri = df_pre.groupby(tri).count().reset_index()
    df_cat = df_pre.groupby(cat).count().reset_index()
    
    if divide_category == True:
        df0_ = make_0_cat_tricat()
        df0_[cat_name] = 'zonk'
        df0 = pd.DataFrame(columns=df0_.columns)
        
        for category in set(df[cat_name]):
            df0_[cat_name] = category
            df0 = df0.append(df0_)
            
    else:
        df0 = make_0_cat_tricat()
    
    df_tri = pd.concat([df_tri, df0], sort=True)
    df_tri.POS = df_tri.POS.astype('int64')
    df_cat = pd.concat([df_cat, df0], sort=True)
    df_cat.POS = df_cat.POS.astype('int64')
    
    df_tri = df_tri.groupby(tri).sum().reset_index()
    df_cat = df_cat.groupby(cat).sum().reset_index()
    
    df_all = df_tri.append(df_cat, ignore_index=True, sort=True)
    df_all.TRI_CAT.fillna(value='6_CAT', inplace = True)
    
    df_all['CATEGORY'] = df_all.apply(lambda row: annotate_category(row),
                                      axis=1)
    
    return df_all

def annotate_category(row):
    '''annotate the category as for the 95_CI files.'''
    
    if row['TRI_CAT'] == '6_CAT':
        return row['CAT']
    
    else:
        tri = row['TRI_CAT']
        cat = row['CAT']
        
        return (tri[0] + '[' + cat + ']' + tri[2])

#plot the data
def cat6_relative_plot(df, divide_category=False, cat_name=None, y=1.):
    '''plot the 6 categories, without information from gnomAD or SSC; can be
    added later as ranges/lines etc. use catplot from sns and the four
    as defined in the figure. use these colors as indicated below.'''
    
    #generates table with all the categories, so nothing is left out
    df = make_table_TRIs_CATs(df, divide_category, cat_name)
    #removes all the TRI_CAT counts   
    df = df[df.TRI_CAT == '6_CAT']
    
    if divide_category == True:
        df['SUM'] = df.groupby(cat_name)['POS'].transform('sum')
    else:
        df['SUM'] = df['POS'].sum()

    df['REL_COUNT'] = df.POS/df.SUM
    
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']
    
    sns.catplot(x='CAT', y='REL_COUNT', hue=cat_name, data=df,
                    kind='bar', palette=colors)
    
    plt.xlabel('Categories')
    plt.ylabel('Relative Contribution')
    plt.ylim(0,y)
    sns.despine(offset=5, trim=True)
    plt.show()

    
def cat96_relative_plot(df, divide_category=False, cat_name=None, y=1.):
    '''plot the 96 categories, without information from gnomAD or SSC; most
    likely this will be integrated in different plot that this one similar to
    what was done for the NMed paper.'''
    
    #generates table with all the categories, so nothing is left out
    df = make_table_TRIs_CATs(df, divide_category, cat_name)
    #removes all the CAT counts   
    df = df[df.TRI_CAT != '6_CAT']
    
    #sum all categories or only one depending on input
    if divide_category == True:
        df['SUM'] = df.groupby(cat_name)['POS'].transform('sum')
    else:
        df['SUM'] = df['POS'].sum()

    df['REL_COUNT'] = df.POS/df.SUM
    
    #make several plots or only one depending on categories
    if divide_category == True:
        f, axs = plt.subplots(nrows=len(df[cat_name].unique()), ncols=6,
                              sharex=True)
    else:
        f, axs = plt.subplots(nrows=2, ncols=6, sharex=True)
    
    
    poss_vars = list(set(df['CAT'].values.tolist()))
    poss_vars.sort()
    
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink']
    
    for i, cat in enumerate(df[cat_name].sort_values().unique()):

        for j, var in enumerate(poss_vars):
            
            g = sns.catplot(x='TRI_CAT', y='REL_COUNT',
                            data=df[(df.CAT == var) & (df[cat_name] == cat)],
                            kind='bar', color=colors[j], ax=axs[i,j])
            
            plt.close(g.fig)
            
            if j == 0:
                axs[i,0].set_ylim(0, y)
                sns.despine(offset=5, trim=True, ax=axs[i,0])
                axs[i,0].set_ylabel('Relative Contribution')
                axs[i,j].set_xlabel(var)
            else:
                sns.despine(offset=5, trim=True, ax=axs[i,j], left=True)
                axs[i,j].set_ylim(0, y)
                axs[i,j].set_yticks([])
                axs[i,j].set_ylabel('')
                    
                axs[i,j].set_xlabel(var)
                axs[i,j].set_xticks([])
                axs[i,j].set_title('')
        
        plt.show()

def tg_plot_huebars(df):
    
    df = df[(df.MUTSIG_CLASS == '1_sperm') & (df.INDEL == False)]
    
    df['BINS'] = pd.cut(df.MAF_SPERM_A, bins=[0., 0.02, 0.03, 0.05, 0.2])
    
    df = df.groupby(['BINS', 'CAT']).count()
    
    df['SUM'] = df.groupby('BINS')['POS'].transform('sum')
    df['REL_COUNT'] = df.POS/df.SUM
    df = df.reset_index()
    #df = df[df.CAT == 'T>G']
    
    sns.catplot(x='CAT', y='REL_COUNT', hue='BINS', data=df, kind='bar')
    plt.show()
    
    
def relative_contribution_AFbins(df):
    
    df = df[(df.MUTSIG_CLASS == '1_sperm') & (df.INDEL == False)]
    
    df['BINS'] = pd.cut(df.MAF_SPERM_A, bins=[0., 0.02, 0.03, 0.05, 0.2])
    
    df = df.groupby(['BINS', 'CAT']).count()
    
    df['SUM'] = df.groupby('BINS')['POS'].transform('sum')
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    df.sort_values(by=['BINS', 'CAT'], ascending=False, inplace=True)
    
    pos = [0, 1, 2, 3]
    width = 1
    names = ['0.20-0.05', '0.05-0.03', '0.03-0.02', '0.02-0.00']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0,0,0,0]
    
    for cat in df.CAT.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.CAT == cat].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.CAT == cat].REL_COUNT)
        
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink'][::-1]
    
    for i in range(6):
        
        plt.bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                edgecolor='white', width=width)
    
    plt.xlabel('AF Bin')
    plt.ylabel('Relative Contribution')
    
    sns.despine(bottom=True, offset=5, trim=True)
    plt.xticks(pos, names, rotation=45)
    plt.show()


def relative_contribution_cohort_sperm(df):
    
    df = df[(df.MUTSIG_CLASS == '1_sperm') & (df.INDEL == False)]
    
    df = df.groupby(['COHORT', 'CAT']).count()
    
    df['SUM'] = df.groupby('COHORT')['POS'].transform('sum')
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    
    #Old does not have T>A
    df = df.append(pd.DataFrame(data={'COHORT': ['Old'], 'CAT': ['T>A'],
                                      'REL_COUNT': [0.]}), ignore_index=True)
    
    df.sort_values(by=['COHORT', 'CAT'], ascending=False, inplace=True)
    
    pos = [0, 1, 2]
    width = 1
    names = ['ID01-12', 'ID13-17', 'F01-08']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0,0,0]
    
    for cat in df.CAT.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.CAT == cat].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.CAT == cat].REL_COUNT)
        
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink'][::-1]
    
    for i in range(6):
        
        plt.bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                edgecolor='white', width=width)
    
    plt.xlabel('Cohort')
    plt.ylabel('Relative Contribution')
    
    sns.despine(bottom=True, offset=5, trim=True)
    plt.xticks(pos, names, rotation=45)
    plt.show()
        

#---
def relative_contribution_blood_y_o_clocol(df):
    
    df = df[(df.MUTSIG_CLASS.isin(['3_soma_y', '4_soma_o'])) &
            (df.INDEL == False)]
    
    df['PLOT_CAT'] = df.apply(plot_cat, axis=1)
    
    df = df.groupby(['PLOT_CAT', 'CAT']).count()
    
    df['SUM'] = df.groupby('PLOT_CAT')['POS'].transform('sum')
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    
    df.sort_values(by=['PLOT_CAT', 'CAT'], ascending=False, inplace=True)
    
    pos = [0, 1, 2]
    width = 1
    names = ['Blood-Y', 'Blood-A (ncc)', 'Blood-A (cc)']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0,0,0]
    
    for cat in df.CAT.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.CAT == cat].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.CAT == cat].REL_COUNT)
        
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink'][::-1]
    
    for i in range(6):
        
        plt.bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                edgecolor='white', width=width)
    
    plt.xlabel('Cohort')
    plt.ylabel('Relative Contribution')
    
    sns.despine(bottom=True, offset=5, trim=True)
    plt.xticks(pos, names, rotation=45)
    plt.show()
    
    
def relative_contribution_blood_y_o_clocol_individual(df):
    
    df = df[(df.MUTSIG_CLASS.isin(['3_soma_y', '4_soma_o'])) &
            (df.INDEL == False)]
    
    df['PLOT_CAT'] = df.apply(plot_cat, axis=1)
    
    df0 = make_0_cat_indplot(df)
    
    df = df.groupby(['PLOT_CAT', 'INDIVIDUAL', 'AGE', 'CAT']).count()
    df.reset_index(inplace=True)
    df = df.append(df0, ignore_index=True, sort=False)
    df = df.groupby(['PLOT_CAT', 'INDIVIDUAL', 'AGE', 'CAT']).sum()
    df.reset_index(inplace=True)
    
    df['SUM'] = df.groupby(['PLOT_CAT', 'INDIVIDUAL'])['POS'].transform('sum')
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    
    df['INVERSE_AGE'] = df.apply(lambda row: - row['AGE'], axis=1)
    
    df.sort_values(by=['PLOT_CAT', 'INVERSE_AGE', 'CAT'], ascending=False,
                   inplace=True)
    
    pos = [i for i in range(25)]
    width = 1
    #names = ['Blood-Y', 'Blood-A (ncc)', 'Blood-A (cc)']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0 for i in range(25)]
    
    for cat in df.CAT.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.CAT == cat].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.CAT == cat].REL_COUNT)
        
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink'][::-1]
    
    for i in range(6):
        
        plt.bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                edgecolor='white', width=width)
    
    plt.xlabel('Cohort')
    plt.ylabel('Relative Contribution')
    
    sns.despine(bottom=True, offset=5, trim=True)
    #plt.xticks(pos, names, rotation=45)
    
    plt.show()
    
def plot_cat(row):
    
    if row['INDIVIDUAL'] in ['ID01', 'ID02', 'ID03', 'ID04', 'ID05', 'ID06',
                             'ID07', 'ID08', 'ID09', 'ID10', 'ID11', 'ID12']:
        
        return '3_young'
    
    elif row['INDIVIDUAL'] in ['F01', 'F03', 'F04', 'F05', 'F06', 'F07', 'F08',
                               'ID13', 'ID15', 'ID16']:
        
        return '2_old_pre'
    
    elif row['INDIVIDUAL'] in ['F02', 'ID14', 'ID17']:
        
        return '1_old_post'
    
    else:
        
        return 'zonk'

def make_0_cat_indplot(df):
    '''makes all 6 cats for Individual/Age etc. after plot_cat is added'''
    
    df = df[['INDIVIDUAL', 'AGE', 'PLOT_CAT']].drop_duplicates()
    
    refs = ['C', 'T']
    nts = ['A', 'T', 'G', 'C']
    
    cats = [ref + '>' + alt for ref in refs for alt in nts]
    cats.remove('T>T')
    cats.remove('C>C')
    
    columns = []
    
    for cat in cats:
        for ind, age, pcat in zip(df.INDIVIDUAL, df.AGE, df.PLOT_CAT):
            columns.append((ind, age, pcat, cat))
    
    df = pd.DataFrame(data=columns, columns=['INDIVIDUAL', 'AGE', 'PLOT_CAT',
                                             'CAT'])
    df['POS'] = 0
    
    return df
#---


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#plot functions for mutsig, but including significances; customized for the 
#problem at hand

def merge_SSC_gnomAD(path_to_SSC, path_to_gnomAD):
    '''function to merge the relative counts for SSC and gnomAD, so the table
    is ready to be merged with the data table. Note that the use is intended
    to gnerate a ssc_gno table with all the annotations here and save it as an
    object that can be passed directly to the function.'''
    
    ssc = pd.read_csv(path_to_SSC)
    gno = pd.read_csv(path_to_gnomAD)
    
    ssc_gno = pd.merge(ssc, gno,
                       on=['MUTSIG_CLASS', 'CAT', 'TRI_CAT', 'CATEGORY'],
                       suffixes=('_ssc', '_gno'))
    
    return ssc_gno

def make_ranges_MUTSIG_CLASS():
    '''makes the postional range for the vlines building on the distribution
    gleaned from ax.axes[0][0].patches[N].get_x, get_width etc. note that the
    order of patches is for this CAT, MUTSIG_CLASS.
    
    Note that it is always 0.8 for each element (i.e. 1 total with 0.1 spacing
    between the ends; total of 0.2 between two elements. divide the 0.8 by n,
    which is the number of columns, to get the spacing. Then adjust to get the
    center.'''
    
    ranges = []
    
    for i in range (6):
        ranges.append(i - 0.3)
        ranges.append(i - 0.1)
        ranges.append(i + 0.1)
        ranges.append(i + 0.3)
        
    return ranges

def cat6_relative_plot_SSCgno(df, ssc_gno, y=1., suffix='_ssc'):
    '''plot the 6 categories, with information from gnomAD or SSC; use catplot
    from sns and the four as defined in the figure. use these colors as
    indicated below. gnomAD ssc is provided as a df generated as indicated
    above. suffix determines whether plot is for gno or ssc. default is ssc,
    but both plots will be needed.'''
    
    #generates table with all the categories, so nothing is left out
    df = make_table_TRIs_CATs(df, True, 'MUTSIG_CLASS')
    #removes all the TRI_CAT counts   
    df = df[df.TRI_CAT == '6_CAT']
    
    df['SUM'] = df.groupby('MUTSIG_CLASS')['POS'].transform('sum')

    df['REL_COUNT'] = df.POS/df.SUM
    
    df_merge = pd.merge(df, ssc_gno)
    
    colors=['g', 'xkcd:brown', 'xkcd:golden', 'xkcd:orangish red']
    
    ranges = make_ranges_MUTSIG_CLASS()
    
    sns.catplot(x='CAT', y='REL_COUNT', hue='MUTSIG_CLASS', data=df_merge,
                    kind='bar', palette=colors)
    
    df_vlines = df_merge.sort_values(by=['CAT', 'MUTSIG_CLASS'])
    
    plt.vlines(ranges, ymin=df_vlines['2.5%'+suffix],
                       ymax=df_vlines['97.5%'+suffix],
               alpha=0.5, linewidths=5, colors='0.8')
    
    plt.xlabel('Categories')
    plt.ylabel('Relative Contribution')
    plt.ylim(0,y)
    sns.despine(offset=5, trim=True)
    plt.show()

#------------------------------------------------------------------------------




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Extra plots, not directly related
    
#------------------------------------------------------------------------------
def het_load(path):
    '''provide path to het variants, most likely
    '2020_03_22_het_variants_for_coverage_final_summary_depth_with_group.vcf'
    and add depth before violin plots are made. use this output to plot violins
    as per function below.'''
    
    hets = pd.read_csv(path, sep='\t')
    hets = hets[hets.COHORT.isin(['Young', 'Old'])]
    
    hets['DEPTH'] = hets.REF_COUNT + hets.ALT_COUNT
    hets['INDIVIDUAL'] = hets.apply(lambda row: row.ID.split('-')[0], axis=1)
    hets['SAMPLE'] = hets.apply(lambda row: '_'.join(row.ID.split('-')[1:]),
                                axis=1)
    replace_INDIVIDUAL_PapID(hets)
    
    return hets
    
    
def het_counts_plot(hets):
    
    sns.violinplot(x='INDIVIDUAL', y='DEPTH', hue='SAMPLE', data=hets,
                   palette=['xkcd:bright orange', 'g', 'g', 'g'])
    
    plt.ylim()
    plt.ylabel('Depth of Coverage')
    sns.despine(offset=5, trim=True)
    
    plt.show()









