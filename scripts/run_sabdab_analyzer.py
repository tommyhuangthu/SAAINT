#!/usr/bin/env python3
import os, sys, datetime, pandas, subprocess, argparse, math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import read_list


def count_pdb_num(df, exclude_pdbs):
    unique_pdbs = []
    for i in range(len(df['pdb'])):
        pdbid = df.loc[i, 'pdb']
        if pdbid not in exclude_pdbs:
            if pdbid not in unique_pdbs:
                unique_pdbs.append(pdbid)
    print(f'num_of_unique_pdbs: {len(unique_pdbs)}')
    return


def count_data_num(df, exclude_pdbs):
    num_data = 0
    for i in range(len(df['pdb'])):
        pdbid = df.loc[i, 'pdb']
        model = df.loc[i, 'model']
        if pdbid not in exclude_pdbs:
            if model == 0 and True:
            #if True:
                num_data += 1
    print(f'num_of_data_entries: {num_data}')
    return


def count_num_of_pdbs_with_paired_vhvl(df, exclude_pdbs):
    unique_pdbs, num_vhvl = [], 0
    for i in range(len(df['pdb'])):
        pdbid = df.loc[i, 'pdb']
        model = df.loc[i, 'model']
        if pdbid not in exclude_pdbs:
            #if model == 0 and df.loc[i, 'Hchain'] != 'NA' and df.loc[i, 'Lchain'] != 'NA':
            if df.loc[i, 'Hchain'] != 'NA' and df.loc[i, 'Lchain'] != 'NA':
                num_vhvl += 1
                if pdbid not in unique_pdbs:
                    unique_pdbs.append(pdbid)
    print(f'num_pdbs_with_paired_vhvl: {len(unique_pdbs)}, num_VHVLs: {num_vhvl}')
    return

def count_num_of_pdbs_with_ag(df, exclude_pdbs):
    unique_pdbs, num_ag_entries = [], 0
    for i in range(len(df['pdb'])):
        model = df.loc[i, 'model']
        ag = df.loc[i, 'antigen_chain']
        pdbid = df.loc[i, 'pdb']
        if pdbid not in exclude_pdbs: 
            if model == 0 and ag != '' and ag != 'NA':
            #if ag != '' and ag != 'NA':
                num_ag_entries += 1
                if pdbid not in unique_pdbs:
                    unique_pdbs.append(pdbid)
    print(f'num_pdbs_with_ag: {len(unique_pdbs)}, num_ag_entries: {num_ag_entries}')
    return


def count_num_of_pdbs_with_affinity(df, exclude_pdbs):
    unique_pdbs, unique_affs, num_aff_all, num_aff_uniq = [], [], 0, 0
    for i in range(len(df['pdb'])):
        model = df.loc[i, 'model']
        aff = df.loc[i, 'affinity']
        pdbid = df.loc[i, 'pdb']
        if pdbid not in exclude_pdbs: 
            if model == 0 and aff and aff != 'None' and aff != 'none':
            #if aff and aff != 'None' and aff != 'none':
                num_aff_all += 1
                if pdbid not in unique_pdbs:
                    unique_pdbs.append(pdbid)
                    unique_affs.append(aff)
                    num_aff_uniq += 1
                elif aff not in unique_affs:
                    unique_affs.append(aff)
                    num_aff_uniq += 1
    print(f'num_pdbs_with_affinity: {len(unique_pdbs)}, num_aff_all: {num_aff_all}, num_aff_uniq: {num_aff_uniq}')
    return


if __name__ == '__main__':
    usage = f'Usage: python {sys.argv[0]} [options] sabdab_summary_all.tsv'
    # add arguments and options
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('sabdab_tsv', type=str, help='Input a sabdab database summary file')
    parser.add_argument('-j', '--job', type=str, choices=[
        '', 'date', 'classification', 'method', 'resolution', 'publication', 'asym_id', 'pdb_num', 'data_num', 
        'ab_spe', 'ab_type', 'HL_inf_res_num', 'HL_chain_len', 'radius',
        'ag_spe', 'ag_type', 'ab_ag_inf_res_num', 'cdr_inf_res_num', 'cdr_inf_res_ratio', 'ag_chain_num', 
        'num_pdbs_with_paired_vhvl', 'num_pdbs_with_ag', 'check_ab_spe', 
        'affinity', 'num_pdbs_with_affinity'], default='', help='Choose a job type for analysis')
    parser.add_argument('-x', '--exclude', type=str, default='', help='Input a list of pdbs to be excluded')
    args = parser.parse_args()
    sabdab_tsv = ''

    # parse arguments and options
    if args.sabdab_tsv:
        sabdab_tsv = args.sabdab_tsv
    else:
        quit()

    job = args.job
    exclude_file = args.exclude
    if exclude_file:
        exclude_pdbs = read_list(exclude_file)
    else:
        exclude_pdbs = []

    if job == '':
        print(f'please specify a job type for analysis')
        parser.print_help()
        quit()

    if job == 'pdb_num':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        count_pdb_num(df, exclude_pdbs)
    elif job == 'data_num':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        count_data_num(df, exclude_pdbs)
    elif job == 'ab_spe':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        plot_ab_species(df)
    elif job == 'ag_spe':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        plot_ag_species(df)
    elif job == 'ab_type':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        plot_ab_type(df)
    elif job == 'ag_type':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        count_ag_type(df)
    elif job == 'ag_chain_num':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        plot_ag_chain_num(df)
    elif job == 'num_pdbs_with_paired_vhvl':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        count_num_of_pdbs_with_paired_vhvl(df, exclude_pdbs)
    elif job == 'num_pdbs_with_ag':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        count_num_of_pdbs_with_ag(df, exclude_pdbs)
    elif job == 'check_ab_spe':
        df = pd.read_excel(sabdab_tsv, index_col=None, keep_default_na=False)
        check_ab_spe(df)
    
    # deal with SAAINT affinity data
    elif job == 'affinity':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        plot_affinity(df)
    elif job == 'num_pdbs_with_affinity':
        df = pd.read_csv(sabdab_tsv, sep='\t', index_col=None, keep_default_na=False)
        count_num_of_pdbs_with_affinity(df, exclude_pdbs)
    

