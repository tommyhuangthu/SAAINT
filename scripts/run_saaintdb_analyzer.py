#!/usr/bin/env python3
import os, sys, datetime, pandas, subprocess, argparse, math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator

sns.color_palette('tab10')

def get_classification_counts(df):
    df_new = df[['PDB_ID', 'Classification']].drop_duplicates(ignore_index=True)
    classification_counts = df_new['Classification'].value_counts()
    print(f'counts for {classification_counts}')
    return

def get_date_range(df):
    df_new = df[['PDB_ID', 'Deposit_date', 'Release_date']].drop_duplicates(ignore_index=True)
    minimum_deposit_date = df['Deposit_date'].min()
    maximum_deposit_date = df['Deposit_date'].max()
    minimum_release_date = df['Release_date'].min()
    maximum_release_date = df['Release_date'].max()
    print(f'minimum_deposit_date: {minimum_deposit_date}, maximum_deposit_date: {maximum_deposit_date}, minimum_release_date: {minimum_release_date}, maximum_release_date: {maximum_release_date}')
    return

def plot_method(df):
    df_new = df[['PDB_ID', 'Method']].drop_duplicates()
    method_counts = df_new['Method'].value_counts()
    print(f'counts for {method_counts}')
    methods, counts = ['X-ray', 'EM', 'Other'], [0, 0, 0]
    for i, method in enumerate(method_counts.index):
        if method == 'x-ray diffraction':
            counts[0] += int(method_counts.iloc[i])
        elif method == 'electron microscopy':
            counts[1] += int(method_counts.iloc[i])
        else:
            counts[2] += int(method_counts.iloc[i])
    print(methods, counts)

    total = sum(counts)
    fig, ax = plt.subplots()
    ax.pie(counts, labels=methods, textprops={'fontsize': 10}, autopct=lambda p: '{:.0f} ({:.1f}%)'.format(p*total/100, p))
    ax.set_title('Experimental method', fontsize=12)
    plt.savefig('method.png', dpi=600, transparent=True)
    plt.close()
    return

def plot_resolution(df):
    df_new = df[['PDB_ID', 'Method', 'Resolution', 'R_free', 'R_work']].drop_duplicates(ignore_index=True)
    for i in range(len(df_new['Method'])):
        if df_new.loc[i, 'Method'] == 'x-ray diffraction':
            df_new.loc[i, 'Method'] = 'X-ray'
        elif df_new.loc[i, 'Method'] == 'electron microscopy':
            df_new.loc[i, 'Method'] = 'EM'
    df_new['Resolution'] = df_new['Resolution'].astype(str)
    df_new['Resolution'] = pd.to_numeric(df_new['Resolution'], errors='coerce')
    df_new['R_free'] = pd.to_numeric(df_new['R_free'], errors='coerce')
    df_new['R_work'] = pd.to_numeric(df_new['R_work'], errors='coerce')
    filter1 = df_new['Resolution'] > 0
    filter11  = df_new['Resolution'] <= 2.5
    
    filter2 = df_new['Method'] == 'X-ray'
    df_new1 = df_new[filter1 & filter2]
    min_resol, max_resol, med_resol, min_rfree, max_rfree, med_rfree, min_rwork, max_rwork, med_rwork = \
            df_new1['Resolution'].min(), df_new1['Resolution'].max(), df_new1['Resolution'].median(), \
            df_new1['R_free'].min(), df_new1['R_free'].max(), df_new1['R_free'].median(), \
            df_new1['R_work'].min(), df_new1['R_work'].max(), df_new1['R_work'].median()
    print(f'for x-ray structures, min_resol: {min_resol}, max_resol: {max_resol}, median_resol: {med_resol}')
    print(f'min_rfree: {min_rfree}, max_rfree: {max_rfree}, median_rfree: {med_rfree}')
    print(f'min_rwork: {min_rwork}, max_rwork: {max_rwork}, median_rwork: {med_rwork}')

    df_new11 = df_new[filter1 & filter11 & filter2]
    n_xray_res_le25 = df_new11['PDB_ID'].nunique()
    print(f'num_xray_resolution_le2.5: {n_xray_res_le25}')

    filter3 = df_new['Method'] == 'EM'
    df_new2 = df_new[filter1 & filter3]
    min_resol, max_resol, med_resol = df_new2['Resolution'].min(), df_new2['Resolution'].max(), df_new2['Resolution'].median()
    print(f'for electron microscopy structures, min_resol: {min_resol}, max_resol: {max_resol}, median_resol: {med_resol}')

    filter12 = df_new['Resolution'] <= 4.5
    df_new21 = df_new[filter1 & filter12 & filter3]
    n_em_res_le45 = df_new21['PDB_ID'].nunique()
    print(f'num_em_resolution_le4.5: {n_em_res_le45}')


    df_new3 = df_new[(filter2 | filter3) & filter1]

    # create resolution violin plot
    plt.figure(figsize=(2,3))
    sns.violinplot(data=df_new3, x='Method', y='Resolution', hue='Method', hue_order=['X-ray', 'EM'], saturation=1, linewidth=1)
    sns.despine()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Method', fontsize=12)
    plt.ylabel(r'Resolution ($\AA$)', fontsize=12)
    plt.tight_layout()
    plt.savefig('resolution_method.png', dpi=600, transparent=True)
    plt.close()
    
    # create rfactor violin plot
    plt.figure(figsize=(2,3))
    df_rfac = df_new3[['R_free', 'R_work']]
    df_rfac = df_rfac.melt(var_name='R-factor', value_name='Value')
    sns.violinplot(data=df_rfac, x='R-factor', y='Value', hue='R-factor', hue_order=['R_free', 'R_work'], saturation=1, linewidth=1)
    sns.despine()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('R-factor', fontsize=12)
    plt.ylabel('Value', fontsize=12)
    plt.tight_layout()
    plt.savefig('resolution_rfactor.png', dpi=600, transparent=True)
    plt.close()

    return

def get_publication_statistics(df):
    df_new = df[['PDB_ID', 'PMID', 'DOI']].drop_duplicates(ignore_index=True)
    filter1 = df_new['PMID'] != 'N.A.'
    filter2 = df_new['DOI'] != 'N.A.'
    df_pmid = df_new[filter1]
    num_uniq_pdbid_has_pmid = df_pmid['PDB_ID'].nunique()
    df_doi = df_new[filter2]
    num_uniq_pdbid_has_doi = df_doi['PDB_ID'].nunique()
    df_both = df_new[filter1 & filter2]
    num_uniq_pdbid_has_both = df_both['PDB_ID'].nunique()
    print(f'num_uniq_pdbid_has_pmid: {num_uniq_pdbid_has_pmid}, num_uniq_pdbid_has_doi: {num_uniq_pdbid_has_doi}, num_uniq_pdbid_has_both: {num_uniq_pdbid_has_both}')
    return

def get_asym_id_statistics(df):
    df_new = df[['PDB_ID', 'Asym_ID_type']].drop_duplicates(ignore_index=True)
    asymid_counts = df['Asym_ID_type'].value_counts()
    return

##################################
# Plot Ab-related properties
##################################
def plot_ab_species(df):
    filter1 = df['H_species'] != 'N.A.'
    df_new1 = df[filter1]
    heavy_species = df_new1['H_species'].value_counts()
    print(f'counts for {heavy_species}')
    H_species, H_counts = [], []
    for i, spe in enumerate(heavy_species.index):
        if i<5:
            H_species.append(spe)
            H_counts.append(int(heavy_species.iloc[i]))
        elif i==5:
            H_species.append('Other')
            H_counts.append(int(heavy_species.iloc[i]))
        else:
            H_counts[5] += int(heavy_species.iloc[i])
    print(H_species, H_counts)
    for i, spe in enumerate(H_species):
        if ' ' in spe: 
            part1 = spe.split(' ')[0]
            part2 = spe.split(' ')[1]
            new_spe = part1[0].upper()+'. '+part2
            H_species[i] = new_spe

    df_new11 = pd.DataFrame({'spe': H_species, 'count': H_counts})
    fig, ax1 = plt.subplots(figsize=(4,3))
    sns.barplot(data=df_new11, ax=ax1, x='spe', y='count', saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel('Heavy-chain species', fontsize=12)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('H_species.png', dpi=600, transparent=True)
    plt.close()

    filter2 = df['L_species'] != 'N.A.'
    df_new2 = df[filter2]
    light_species = df_new2['L_species'].value_counts()
    print(f'counts for {light_species}')
    L_species, L_counts = [], []
    for i, spe in enumerate(light_species.index):
        if i<5:
            L_species.append(spe)
            L_counts.append(int(light_species.iloc[i]))
        elif i==5:
            L_species.append('Other')
            L_counts.append(int(light_species.iloc[i]))
        else:
            L_counts[5] += int(light_species.iloc[i])
    print(L_species, L_counts)
    for i, spe in enumerate(L_species):
        if ' ' in spe: 
            part1 = spe.split(' ')[0]
            part2 = spe.split(' ')[1]
            new_spe = part1[0].upper()+'. '+ part2
            L_species[i] = new_spe

    df_new21 = pd.DataFrame({'spe': L_species, 'count': L_counts})
    fig, ax1 = plt.subplots(figsize=(4,3))
    sns.barplot(data=df_new21, ax=ax1, x='spe', y='count', saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel('Light-chain species', fontsize=12)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('L_species.png', dpi=600, transparent=True)
    plt.close()

    return


def plot_ab_type(df):
    fig, ax1 = plt.subplots(figsize=(13,4))
    sns.countplot(data=df, ax=ax1, x='Ab_type', order=df['Ab_type'].value_counts().index, saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel('Antibody type', fontsize=12)
    ax1.yaxis.minorticks_on()
    ax1.tick_params('y')
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ab_type.png', dpi=600, transparent=True)
    plt.close()

    return



def plot_HL_chain_len(df):
    H_types, H_seq_lens, H_seq_types = [], [], []
    L_types, L_seq_lens, L_seq_types = [], [], []
    for i in range(len(df['H_pdb_seq_len'])):
        ab_type = df.loc[i, 'Ab_type']
        H_fas_seq_len = int(df.loc[i, 'H_fas_seq_len'])
        H_pdb_seq_len = int(df.loc[i, 'H_pdb_seq_len'])
        H_filled_seq_len = int(df.loc[i, 'H_filled_seq_len'])
        L_fas_seq_len = int(df.loc[i, 'L_fas_seq_len'])
        L_pdb_seq_len = int(df.loc[i, 'L_pdb_seq_len'])
        L_filled_seq_len = int(df.loc[i, 'L_filled_seq_len'])
        if H_pdb_seq_len>0 and L_pdb_seq_len>0:
            H_type = ab_type.split(':')[0]
            H_types.append(H_type)
            H_types.append(H_type)
            H_types.append(H_type)
            H_seq_lens.append(H_fas_seq_len)
            H_seq_lens.append(H_pdb_seq_len)
            H_seq_lens.append(H_filled_seq_len)
            H_seq_types.append('FASTA-seq')
            H_seq_types.append('PDB-seq')
            H_seq_types.append('Filled-seq')
            L_type = ab_type.split(':')[1]
            L_types.append(L_type)
            L_types.append(L_type)
            L_types.append(L_type)
            L_seq_lens.append(L_fas_seq_len)
            L_seq_lens.append(L_pdb_seq_len)
            L_seq_lens.append(L_filled_seq_len)
            L_seq_types.append('FASTA-seq')
            L_seq_types.append('PDB-seq')
            L_seq_types.append('Filled-seq')
        elif H_pdb_seq_len>0 and L_pdb_seq_len==0:
            H_type = ab_type
            H_types.append(H_type)
            H_types.append(H_type)
            H_types.append(H_type)
            H_seq_lens.append(H_fas_seq_len)
            H_seq_lens.append(H_pdb_seq_len)
            H_seq_lens.append(H_filled_seq_len)
            H_seq_types.append('FASTA-seq')
            H_seq_types.append('PDB-seq')
            H_seq_types.append('Filled-seq')
        elif L_pdb_seq_len>0 and H_pdb_seq_len==0:
            L_type = ab_type
            L_types.append(L_type)
            L_types.append(L_type)
            L_types.append(L_type)
            L_seq_lens.append(L_fas_seq_len)
            L_seq_lens.append(L_pdb_seq_len)
            L_seq_lens.append(L_filled_seq_len)
            L_seq_types.append('FASTA-seq')
            L_seq_types.append('PDB-seq')
            L_seq_types.append('Filled-seq')

    flierprops=dict(marker='x', markeredgecolor='gray', markersize=5, markeredgewidth=0.5)
    sns.set(rc={'legend.fontsize': 9})
    #sns.set_style('darkgrid')
    sns.set_style('ticks')
    
    df_new = pd.DataFrame({'H_types': H_types, 'H_seq_lens': H_seq_lens, 'Sequence type': H_seq_types})
    fig, ax1 = plt.subplots(figsize=(8,8))
    ax1.yaxis.minorticks_on()
    ax1.yaxis.grid(True, linestyle='--')
    #plt.figure(figsize=(8,8))
    sns.boxplot(data=df_new, x='H_types', y='H_seq_lens', palette='Accent', hue='Sequence type', 
            hue_order = ['PDB-seq', 'Filled-seq', 'FASTA-seq'], saturation=1, dodge=True, gap=0.15, linewidth=1, flierprops=flierprops)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.xlabel('Heavy-chain antibody type', fontsize=12)
    plt.ylabel('Sequence length', fontsize=12)
    plt.tight_layout()
    plt.savefig('H_chain_len.png', dpi=600)
    plt.close()

    df_new2 = pd.DataFrame({'L_types': L_types, 'L_seq_lens': L_seq_lens, 'Sequence type': L_seq_types})
    fig, ax2 = plt.subplots(figsize=(5,8))
    ax2.yaxis.minorticks_on()
    ax2.yaxis.grid(True, linestyle='--')
    #plt.figure(figsize=(5,8))
    sns.boxplot(data=df_new2, x='L_types', y='L_seq_lens', palette='Accent', hue='Sequence type', 
            hue_order = ['PDB-seq', 'Filled-seq', 'FASTA-seq'], saturation=1, dodge=True, gap=0.15, linewidth=1, flierprops=flierprops)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.xlabel('Light-chain antibody type', fontsize=12)
    plt.ylabel('Sequence length', fontsize=12)
    plt.tight_layout()
    plt.savefig('L_chain_len.png', dpi=600)
    plt.close()
    return


def plot_HL_inf_res_num(df):
    filter1 = df['HL_inf_res_num'] > 0
    #counts = df['Ab_type'].value_counts()
    #filter2 = df['Ab_type'].isin(counts[counts>10].index)
    
    flierprops=dict(marker='x', markeredgecolor='gray', markersize=5, markeredgewidth=0.5)
    #sns.set_style('darkgrid')
    sns.set_style('ticks')
    
    df_new = df[filter1]
    df_res = df_new[['Ab_type', 'HL_inf_res_num']]
    fig, ax1 = plt.subplots(figsize=(8,8))
    ax1.yaxis.minorticks_on()
    ax1.yaxis.grid(True, linestyle='--')
    #plt.figure(figsize=(6,7))
    sns.boxplot(data=df_res, x='Ab_type', y='HL_inf_res_num', color='tab:blue', saturation=1, flierprops=flierprops, linewidth=1)
    #sns.despine()
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.xlabel('Antibody type', fontsize=12)
    plt.ylabel('Number of heavy chain-light chain interface residues', fontsize=12)
    plt.tight_layout()
    plt.savefig('HL_inf_res_num.png', dpi=600)
    plt.close()
    return


def plot_radius(df):
    abtypes, radii = [], []
    for i in range(len(df['Ab_type'])):
        abtype = df.loc[i, 'Ab_type']
        H_pdb_seq_len = df.loc[i, 'H_pdb_seq_len']
        L_pdb_seq_len = df.loc[i, 'L_pdb_seq_len']
        H_radius = float(df.loc[i, 'H_mean_radius'])
        L_radius = float(df.loc[i, 'L_mean_radius'])
        if H_pdb_seq_len>0 and L_pdb_seq_len>0:
            H_type = abtype.split(':')[0]
            L_type = abtype.split(':')[1]
        elif H_pdb_seq_len>0 and L_pdb_seq_len==0:
            H_type = abtype
            L_type = 'N.A.'
        elif H_pdb_seq_len==0 and L_pdb_seq_len>0:
            L_type = abtype
            H_type = 'N.A.'
        if 'scFv' in H_type or 'VHVL' in H_type:
            abtypes.append(H_type)
            radii.append(H_radius)
        if 'scFv' in L_type or 'VHVL' in L_type:
            abtypes.append(L_type)
            radii.append(L_radius)
    
    flierprops=dict(marker='x', markeredgecolor='gray', markersize=5, markeredgewidth=0.5)
    #sns.set_style('darkgrid')
    sns.set_style('ticks')
    
    df_new = pd.DataFrame({'type': abtypes, 'radii': radii})
    fig, ax1 = plt.subplots(figsize=(2,4))
    ax1.yaxis.minorticks_on()
    ax1.yaxis.grid(True, linestyle='--')
    #plt.figure(figsize=(2, 2))
    sns.boxplot(data=df_new, x='type', y='radii', color='tab:blue', saturation=1, dodge=True, gap=0.15, linewidth=1, flierprops=flierprops, order=sorted(df_new['type'].unique()))
    #sns.despine()
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.xlabel('Type', fontsize=12)
    plt.ylabel(r'Radius ($\AA$)', fontsize=12)
    plt.tight_layout()
    plt.savefig('radius.png', dpi=600)
    plt.close()

    sns.reset_defaults()
    return

########################################
# Plot Ag-related properties
########################################
def plot_ag_species(df):
    col_name = 'Ag_species'
    filter1 = df[col_name] != 'N.A.'
    df_new = df[filter1].drop_duplicates(ignore_index=True)
    ag_species = []
    for i in range(len(df_new[col_name])):
        for spe in df_new.loc[i, col_name].split(';'):
            ag_species.append(spe)
    df_new1 = pd.DataFrame({'ag_species': ag_species})
    species_counts = df_new1['ag_species'].value_counts()
    species, counts = [], []
    for i, spe in enumerate(species_counts.index):
        if i<5:
            species.append(spe)
            counts.append(int(species_counts.iloc[i]))
        elif i==5:
            species.append('Other')
            counts.append(int(species_counts.iloc[i]))
        else:
            counts[5] += int(species_counts.iloc[i])
    #print(f'species: {species}, counts: {counts}')
    for i, spe in enumerate(species):
        if spe == 'homo sapiens':
            species[i] = 'H. sapiens'
        elif spe == 'severe acute respiratory syndrome coronavirus 2':
            species[i] = 'SARS-CoV-2'
        elif spe == 'human immunodeficiency virus 1':
            species[i] = 'HIV-1'
        elif spe == 'influenza a virus':
            species[i] = 'Influenza A'
        elif species[i] == 'plasmodium falciparum':
            species[i] = 'P. falciparum'
    
    df_new2 = pd.DataFrame({'species': species, 'count': counts})
    fig, ax1 = plt.subplots(figsize=(3.5,3))
    ax1 = sns.barplot(data=df_new2, x='species', y='count', saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel('Ag species', fontsize=12)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ag_species.png', dpi=600, transparent=True)
    plt.close()
    return

def plot_ag_type(df):
    col_name = 'Ag_type(s)'
    filter1 = df[col_name] != 'N.A.'
    df_new = df[filter1].drop_duplicates(ignore_index=True)
    ag_type = ['protein', 'peptide', 'DNA', 'RNA']
    ag_type_counts = [0 for i in range(4)]
    for i in range(len(df_new[col_name])):
        if ag_type[0] in df_new.loc[i, col_name]:
            ag_type_counts[0] += 1
        if ag_type[1] in df_new.loc[i, col_name]:
            ag_type_counts[1] += 1
        if ag_type[2] in df_new.loc[i, col_name]:
            ag_type_counts[2] += 1
        if ag_type[3] in df_new.loc[i, col_name]:
            ag_type_counts[3] += 1
    print(f'counts for {col_name}:\n{ag_type}: {ag_type_counts}')

    df_new1 = pd.DataFrame({'ag_type': ag_type, 'count': ag_type_counts})
    fig, ax1 = plt.subplots(figsize=(3.5,3))
    ax1 = sns.barplot(data=df_new1, x='ag_type', y='count', saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel('Ag type', fontsize=12)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ag_type.png', dpi=600, transparent=True)
    plt.close()
    return

def plot_ag_chain_num(df):
    col_name = 'Ag_chain_ID(s)'
    filter1 = df[col_name] != 'N.A.'
    df_new = df[filter1].drop_duplicates(ignore_index=True)
    ag_nums = ['1', '2', '3', r'$\geq$4']
    counts = [0, 0, 0, 0]
    for i in range(len(df_new[col_name])):
        num = len(df_new.loc[i, col_name].split(';'))
        if num == 1:
            counts[0] += 1
        elif num == 2:
            counts[1] += 1
        elif num == 3:
            counts[2] += 1
        else:
            counts[3] += 1

    df_new1 = pd.DataFrame({'ag_nums': ag_nums, 'count': counts})
    fig, ax1 = plt.subplots(figsize=(3.5,3))
    ax1 = sns.barplot(data=df_new1, x='ag_nums', y='count', saturation=1)
    sns.despine()
    ax1.bar_label(ax1.containers[0], fontsize=9)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xlabel(r'$N_\mathrm{Ag\ chains}$', fontsize=12)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ag_chain_num.png', dpi=600, transparent=True)
    plt.close()
    return

def plot_ab_ag_inf_res_num(df):
    for i in range(len(df['Ab_ag_inf_res_num'])):
        strs = df['Ab_ag_inf_res_num'][i].split(';')
        ints = [int(x) for x in strs]
        tot_inf_res_num = int(np.sum(np.array(ints)))
        df.loc[i, 'Ab_ag_inf_res_num'] = tot_inf_res_num
    
    filter1 = df['Ab_ag_inf_res_num'] > 0
    df_new = df[filter1]
    df_new1 = df_new['Ab_ag_inf_res_num']
    df_new1_int = df_new1.astype(int)
    min_val = np.min(df_new1_int)
    max_val = np.max(df_new1_int)
    bins = np.arange(min_val - 0.5, max_val + 1.5, 1)
    fig, ax1 = plt.subplots(figsize=(12,3))
    plt.hist(df_new1_int, bins=bins, color='tab:blue', edgecolor='black')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$N_\mathrm{Ab-Ag\ interface\ residues}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ab_ag_inf_res_num.png', dpi=600, transparent=True)
    plt.close()
    
    filter2 = df['Ab_ag_inf_res_num'] > 100
    df_new = df[filter2]
    df_new2 = df_new['Ab_ag_inf_res_num']
    df_new2_int = df_new2.astype(int)
    min_val = np.min(df_new2_int)
    max_val = np.max(df_new2_int)
    bins = np.arange(min_val - 0.5, max_val + 1.5, 1)
    fig, ax2 = plt.subplots(figsize=(4,2))
    plt.hist(df_new2_int, bins=bins, color='tab:orange', edgecolor='black')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$N_\mathrm{Ab-Ag\ interface\ residues}$', fontsize=12)
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('ab_ag_inf_res_num_zoomin.png', dpi=600, transparent=True)
    plt.close()
    return df

def plot_cdr_inf_res_num(df):
    for i in range(len(df['CDR_inf_res_num'])):
        strs = df['CDR_inf_res_num'][i].split(';')
        ints = [int(x) for x in strs]
        mean_cdr_inf_res_num = float(np.mean(np.array(ints)))
        df.loc[i, 'CDR_inf_res_num'] = mean_cdr_inf_res_num
    
    filter1 = df['CDR_inf_res_num'] > 0
    df_new = df[filter1]
    df_new1 = df_new['CDR_inf_res_num']
    df_new1_int = df_new1.astype(int)
    min_val = np.min(df_new1_int)
    max_val = np.max(df_new1_int)
    bins = np.arange(min_val - 0.5, max_val + 1.5, 1)
    fig, ax1 = plt.subplots(figsize=(12, 3))
    plt.hist(df_new1_int, bins=bins, color='tab:blue', edgecolor='black')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$N_\mathrm{CDR\ residues\ at\ Ab-Ag\ interface}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('cdr_inf_res_num.png', dpi=600, transparent=True)
    plt.close()

    filter2 = df['CDR_inf_res_num'] > 30
    df_new = df[filter2]
    df_new2 = df_new['CDR_inf_res_num']
    df_new2_int = df_new2.astype(int)
    min_val = np.min(df_new2_int)
    max_val = np.max(df_new2_int)
    bins = np.arange(min_val - 0.5, max_val + 1.5, 1)
    fig, ax2 = plt.subplots(figsize=(4, 2))
    plt.hist(df_new2_int, bins=bins, color='tab:orange', edgecolor='black')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$N_\mathrm{CDR\ residues\ at\ Ab-Ag\ interface}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('cdr_inf_res_num_zoomin.png', dpi=600, transparent=True)
    plt.close()

    return df

def plot_cdr_inf_res_ratio(df):
    for i in range(len(df['CDR_inf_res_ratio'])):
        strs = df['CDR_inf_res_ratio'][i].split(';')
        ints = [float(x) for x in strs]
        mean_cdr_inf_res_ratio = float(np.mean(np.array(ints)))
        df.loc[i, 'CDR_inf_res_ratio'] = mean_cdr_inf_res_ratio
    
    filter1 = df['CDR_inf_res_ratio'] > 0
    df_new = df[filter1]
    df_new1 = df_new['CDR_inf_res_ratio']
    df_new1_float = df_new1.astype(float)
    min_val = np.min(df_new1_float)
    max_val = np.max(df_new1_float)
    bins = np.arange(min_val, max_val+0.01, 0.025)
    fig, ax1 = plt.subplots(figsize=(12, 3))
    plt.hist(df_new1_float, bins=bins, color='tab:blue', edgecolor='black')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$R_\mathrm{CDR\ residues\ at\ Ab-Ag\ interface}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('cdr_inf_res_ratio.png', dpi=600, transparent=True)
    plt.close()

    filter2 = df['CDR_inf_res_ratio'] <= 0.4
    df_new = df[filter1 & filter2]
    df_new2 = df_new['CDR_inf_res_ratio']
    df_new2_float = df_new2.astype(float)
    min_val = np.min(df_new2_float)
    max_val = np.max(df_new2_float)
    bins = np.arange(min_val, max_val+0.01, 0.025)
    fig, ax2 = plt.subplots(figsize=(4, 2))
    plt.hist(df_new2_float, bins=bins, color='tab:orange', edgecolor='black')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'$R_\mathrm{CDR\ residues\ at\ Ab-Ag\ interface}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('cdr_inf_res_ratio_zoomin.png', dpi=600, transparent=True)
    plt.close()



    return df

def plot_pdb_num(df):
    df_new = df[['PDB_ID', 'Deposit_date']].drop_duplicates(ignore_index=True)
    year_count = dict()
    for i in range(len(df_new['PDB_ID'])):
        year = df_new.loc[i, 'Deposit_date'][0:4]
        if year not in year_count:
            year_count[year] = 1
        else:
            year_count[year] += 1
    year_count = dict(sorted(year_count.items()))
    
    years, counts, counts_acc = [], [], []
    for key in sorted(year_count.keys()):
        val = year_count[key]
        years.append(int(key))
        counts.append(val)
        counts_acc.append(int(np.sum(np.array(counts))))
    print(f'years: {years}, counts: {counts}, counts_acc: {counts_acc}')

    min_year = min(years)
    max_year = max(years)

    fig, ax1 = plt.subplots(figsize=(18,6))
    ax1.set_xlim(min_year-0.3, max_year+0.7)
    ax1.bar(years, counts_acc, color='white', edgecolor='black')
    ax1.set_ylabel('Cumulative PDB count in SAAINT-DB', fontsize=14, color='black')
    ax1.spines[['right']].set_visible(False)
    #ax1.spines['left'].set_color('tab:blue')
    ax1.tick_params('y', labelsize=12, colors='black')
    ax1.yaxis.minorticks_on()
    import matplotlib.ticker as ticker
    tick_spacing = 1
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax1.tick_params('x', labelsize=12, rotation=45)
    ax1.set_xlabel('Year of PDB deposition', fontsize=14)
    ax2 = ax1.twinx()
    ax2.plot(years, counts, color='tab:orange')
    ax2.set_ylabel('Yearly PDB count in SAAINT-DB', fontsize=12, color='tab:orange')
    ax2.spines[['left']].set_visible(False)
    ax2.spines['right'].set_color('tab:orange')
    ax2.tick_params('y', labelsize=12, colors='tab:orange')
    ax2.yaxis.minorticks_on()
    ax2.tick_params('y', which='minor', color='tab:orange')
    plt.tight_layout()
    plt.savefig('pdb_num.png', dpi=600, transparent=True)
    plt.close()
    return

def plot_affinity(df):
    df_new1 = df[df['Affinity_KD(nM)'] != 'N.A.'].drop_duplicates(ignore_index=True)
    kds_all, kds_le6, kds_ge11 = [], [], []
    for i in range(len(df_new1['Affinity_KD(nM)'])):
        kd = df_new1.loc[i, 'Affinity_KD(nM)']
        kd = kd.replace('>', '').replace('<', '').replace('*', '')
        val = 9-math.log10(float(kd))
        kds_all.append(val)
        if val <= 6: kds_le6.append(val)
        elif val >= 11: kds_ge11.append(val)
    print(f'num_of_affinity_data: {len(kds_all)}')

    fig, ax1 = plt.subplots(figsize=(12,3))
    bins=np.arange(4, 14.1, 0.25)
    plt.hist(kds_all, bins=bins, color='tab:blue', edgecolor='black')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'p$K_\mathrm{D}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('affinity_pkd.png', dpi=600, transparent=True)
    plt.close()

    fig, ax2 = plt.subplots(figsize=(2,2))
    bins=np.arange(4, 6.1, 0.25)
    plt.hist(kds_le6, bins=bins, color='tab:orange', edgecolor='black')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'p$K_\mathrm{D}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('affinity_pkd_le6.png', dpi=600, transparent=True)
    plt.close()

    fig, ax3 = plt.subplots(figsize=(3,2))
    bins=np.arange(11, 14.1, 0.25)
    plt.hist(kds_ge11, bins=bins, color='tab:orange', edgecolor='black')
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    plt.ylabel('Count', fontsize=12)
    plt.xlabel(r'p$K_\mathrm{D}$', fontsize=12)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('affinity_pkd_ge11.png', dpi=600, transparent=True)
    plt.close()






    df_new2 = df[df['Affinity_method'] != 'N.A.'].drop_duplicates(ignore_index=True)
    method_counts = df_new2['Affinity_method'].value_counts()
    #print(method_counts)
    methods, counts = [], []
    for i, method in enumerate(method_counts.index):
        if i<5:
            methods.append(method)
            counts.append(int(method_counts.iloc[i]))
        elif i==5:
            methods.append('Other')
            counts.append(int(method_counts.iloc[i]))
        else:
            counts[5] += int(method_counts.iloc[i])
    for i, method in enumerate(methods):
        if method == 'microscale thermophoresis':
            methods[i] = 'MST'
        elif method == 'grating coupled interferometry':
            methods[i] = 'GCI'
        elif method == 'kinetic exclusion assay':
            methods[i] = 'KEA'
        elif method == 'flow cytometry':
            methods[i] = 'FC'
        elif method == 'fluorescence anisotropy titration':
            methods[i] = 'FAT'
        elif method == 'fluorescence polarization':
            methods[i] = 'FP'
        elif method == 'mass photometry':
            methods[i] = 'MP'
        elif method == 'solution equilibrium titration':
            methods[i] = 'SET'
        elif method == 'mass photometry':
            methods[i] = 'MP'
        elif method == 'electrochemiluminescence multiplex assay':
            methods[i] = 'ECIMA'
        elif method == 'yeast surface display':
            methods[i] = 'YSD'
        elif method == 'scatchard analysis':
            methods[i] = 'SA'

    df_new21 = pd.DataFrame({'methods': methods, 'count': counts})
    fig, ax2 = plt.subplots(figsize=(3.5,3))
    ax2 = sns.barplot(data=df_new21, x='methods', y='count', color='tab:blue', saturation=1)
    sns.despine()
    ax2.bar_label(ax2.containers[0], fontsize=9)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_xlabel('Affinity Method', fontsize=12)
    plt.xticks(rotation=30, fontsize=10, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig('affinity_methods.png', dpi=600, transparent=True)
    plt.close()
    return


def count_num_entries(df):
    unique_pdbs, data_num = [], 0
    for i in range(len(df['PDB_ID'])):
        data_num += 1
        pdbid = df.loc[i, 'PDB_ID']
        if pdbid not in unique_pdbs:
            unique_pdbs.append(pdbid)
    print(f'num_of_unique_pdbs: {len(unique_pdbs)}, num_of_data_entries: {data_num}')
    return


def count_num_entries_with_paired_vhvl(df):
    df_new = df[df['Ab_type'] != 'N.A.'].drop_duplicates(ignore_index=True)
    pdbs, num_vhvl = [], 0
    for i in range(len(df_new['Ab_type'])):
        ab_type = df_new.loc[i, 'Ab_type']
        HL_inf_res_num = df_new.loc[i, 'HL_inf_res_num']
        if ':' in ab_type or 'scFv' in ab_type or 'VHVL+' in ab_type or HL_inf_res_num>0:
            num_vhvl += 1
            pdb = df_new.loc[i, 'PDB_ID']
            if pdb not in pdbs:
                pdbs.append(pdb)
    print(f'num_pdbs_with_paired_vhvl: {len(pdbs)}, num_entries_with_paired_vhvl: {num_vhvl}')
    return

def count_num_entries_with_ag(df):
    df_new = df[df['Ag_type(s)'] != 'N.A.'].drop_duplicates(ignore_index=True)
    pdbs, num_ag_entries = [], 0
    for i in range(len(df_new['Ag_type(s)'])):
        num_ag_entries += 1
        pdb = df_new.loc[i, 'PDB_ID']
        if pdb not in pdbs:
            pdbs.append(pdb)
    print(f'num_pdbs_with_ag: {len(pdbs)}, num_entries_with_ag: {num_ag_entries}')
    return


def count_num_entries_with_affinity(df):
    df_new = df[df['Affinity_KD(nM)'] != 'N.A.'].drop_duplicates(ignore_index=True)
    pdbs, num_aff = [], 0
    for i in range(len(df_new['Affinity_KD(nM)'])):
        num_aff += 1
        pdb = df_new.loc[i, 'PDB_ID']
        if pdb not in pdbs:
            pdbs.append(pdb)
    print(f'num_pdbs_with_affinity: {len(pdbs)}, num_entries_with_affinity: {num_aff}')
    return

def check_ab_spe(df):
    H_types, H_seq_lens, H_seq_types = [], [], []
    L_types, L_seq_lens, L_seq_types = [], [], []
    for i in range(len(df['H_species'])):
        pdbid = df.loc[i, 'PDB_ID']
        H_chain_id = df.loc[i, 'H_chain_ID']
        L_chain_id = df.loc[i, 'L_chain_ID']
        H_spe = df.loc[i, 'H_species']
        L_spe = df.loc[i, 'L_species']
        if H_spe != 'N.A.' and L_spe != 'N.A.' and H_spe != L_spe:
            string = '\t'.join([pdbid, H_chain_id, L_chain_id, H_spe, L_spe])
            print(string)

    return




if __name__ == '__main__':
    usage = f'Usage: python {sys.argv[0]} [options] saaint_database.xlsx'
    # add arguments and options
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('saaint_file', type=str, help='Input a saaint database summary file')
    parser.add_argument('-j', '--job', type=str, choices=[
        '', 'date', 'classification', 'method', 'resolution', 'publication', 'asym_id', 'plot_pdb_num', 'num_entries', 
        'ab_spe', 'ab_type', 'HL_inf_res_num', 'HL_chain_len', 'radius',
        'ag_spe', 'ag_type', 'ab_ag_inf_res_num', 'cdr_inf_res_num', 'cdr_inf_res_ratio', 'ag_chain_num', 
        'num_entries_with_paired_vhvl', 'num_entries_with_ag', 'check_ab_spe', 'plot_all', 
        'affinity', 'num_entries_with_affinity'], default='', help='Choose a job type for analysis')
    args = parser.parse_args()
    saaint_file = ''

    # parse arguments and options
    if args.saaint_file:
        saaint_file = args.saaint_file
    else:
        quit()

    job = args.job
    if job == '':
        print(f'please specify a job type for analysis')
        parser.print_help()
        quit()

    if job == 'plot_all':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        get_date_range(df)
        plot_pdb_num(df)
        plot_method(df)
        plot_resolution(df)
        plot_ab_species(df)
        plot_ab_type(df)
        plot_HL_chain_len(df)
        plot_HL_inf_res_num(df)
        plot_radius(df)
        plot_ag_chain_num(df)
        plot_ag_type(df)
        plot_ag_species(df)
        plot_ab_ag_inf_res_num(df)
        plot_cdr_inf_res_num(df)
        plot_cdr_inf_res_ratio(df)
    
    elif job == 'plot_pdb_num':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_pdb_num(df)
    elif job == 'num_entries':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        count_num_entries(df)
    elif job == 'method':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_method(df)
    elif job == 'resolution':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_resolution(df)
    elif job == 'classification':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        get_classification_counts(df)
    elif job == 'publication':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        get_publication_statistics(df)
    elif job == 'date':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        get_date_range(df)
    elif job == 'HL_inf_res_num':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_HL_inf_res_num(df)
    elif job == 'asym_id':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        get_asym_id_counts(df)
    elif job == 'ab_spe':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ab_species(df)
    elif job == 'ag_spe':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ag_species(df)
    elif job == 'ab_type':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ab_type(df)
    elif job == 'ag_type':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ag_type(df)
    elif job == 'ag_chain_num':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ag_chain_num(df)
    elif job == 'ab_ag_inf_res_num':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_ab_ag_inf_res_num(df)
    elif job == 'cdr_inf_res_num':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_cdr_inf_res_num(df)
    elif job == 'cdr_inf_res_ratio':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_cdr_inf_res_ratio(df)
    elif job == 'HL_chain_len':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_HL_chain_len(df)
    elif job == 'radius':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        plot_radius(df)
    
    elif job == 'num_entries_with_paired_vhvl':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        count_num_entries_with_paired_vhvl(df)
    elif job == 'num_entries_with_ag':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        count_num_entries_with_ag(df)
    elif job == 'check_ab_spe':
        df = pd.read_excel(saaint_file, index_col=None, keep_default_na=False)
        check_ab_spe(df)
    
    # deal with SAAINT affinity data
    elif job == 'affinity':
        df = pd.read_csv(saaint_file, sep='\t', index_col=None, keep_default_na=False)
        plot_affinity(df)
    elif job == 'num_entries_with_affinity':
        df = pd.read_csv(saaint_file, sep='\t', index_col=None, keep_default_na=False)
        count_num_entries_with_affinity(df)
    

