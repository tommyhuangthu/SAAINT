#!/usr/bin/env python3
import os, sys, argparse, re, warnings, subprocess, gzip, glob, string, random, copy, math
from Bio import PDB, BiopythonWarning, Align
from Bio.Align import substitution_matrices
from utils import read_list, read_fasta, write_fasta, Entity, SelectChains, \
        run_abrsa_pdb, run_unidesign_find_interface_residues, \
        parse_pdb_fasta_header, run_faspr_repack_pdb, run_pulchra, reformat_pulchra_rebuilt_pdb, \
        determine_chain_type_by_pdb_content, reindex_pdb, \
        extract_pdb_chains_to_file, pdb_to_fasta, reindex_pdb_by_list, calculate_mean_radius
from fetch_pdb_web_info import fetch_pdb_web_info
from utils import unidesign, faspr, pulchra, abrsa, abrsa_pdb, abalign, abalign_lib, tmalign, struct_ref_vl, struct_ref_vh
import numpy as np

verbose      = False
asym_id_type = 'auto'

max_len_pep = 50
min_len_ab_chain = 60
min_tmscore_ab_domain = 0.4

min_len_vh = 80
max_len_vh = 150
min_len_vl = 80
max_len_vl = 150
min_len_fab = 180
max_len_fab = 260
min_len_vhvl = 180
max_len_vhvl = 280
max_radius_scfv = 20.0
max_pdb_seq_diff = 60

min_inf_res_num_vhm_vlm = 8
min_inf_res_num_vh_vl = 20
min_inf_res_num_fabh_vl = 35
min_inf_res_num_vh_fabl = 35
min_inf_res_num_fabh_fabl = 35
min_inf_res_num_vhvl_vhvl = 80
max_cdr_inf_res_ratio = 0.7


def string_contains_substr(string, substr_list):
    for substr in substr_list:
        if substr in string:
            return True
    return False


def update_HL_ab_type(all_paired_HLs, pdb_entities, pdbid, title):
    for paired_HLs in all_paired_HLs:
        for paired_HL in paired_HLs:
            H_chain_id, L_chain_id = paired_HL[0], paired_HL[1]
            type_H, type_L = 'N.A.', 'N.A.'
            if H_chain_id != '':
                ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
                filled_seq_len_H, pdb_seq_len_H = ent_H.get_filled_pdb_len(index_H), ent_H.get_real_pdb_len(index_H)
                if filled_seq_len_H - pdb_seq_len_H > max_pdb_seq_diff: filled_seq_len_H = pdb_seq_len_H
                spe_H, abrsa_H, name_H, radius_H = ent_H.get_species(), ent_H.get_abrsa_type(index_H), ent_H.get_name(), ent_H.get_mean_radius(index_H)
                if verbose:
                    print(f'H_chain_id: {H_chain_id}, pdb_seq_len_H: {pdb_seq_len_H}, filled_seq_len_H: {filled_seq_len_H}, name_H: {name_H}, title: {title}, len_all_paired_HLs: {len(all_paired_HLs)}')
            if L_chain_id != '':
                ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
                filled_seq_len_L, pdb_seq_len_L = ent_L.get_filled_pdb_len(index_L), ent_L.get_real_pdb_len(index_L)
                if filled_seq_len_L - pdb_seq_len_L > max_pdb_seq_diff: filled_seq_len_L = pdb_seq_len_L
                spe_L, abrsa_L, name_L, radius_L = ent_L.get_species(), ent_L.get_abrsa_type(index_L), ent_L.get_name(), ent_L.get_mean_radius(index_L)
                if verbose:
                    print(f'L_chain_id: {L_chain_id}, pdb_seq_len_L: {pdb_seq_len_L}, filled_seq_len_L: {filled_seq_len_L}, name_L: {name_L}, title: {title}, len_all_paired_HLs: {len(all_paired_HLs)}')
            # infer the most possible ab type for heavy chain
            if H_chain_id != '':
                if abrsa_H == 'heavy':
                    spe_strs = ['lama', 'vicugna', 'camel']
                    name_strs1 = ['nanobo', 'vhh', 'nano-bo']
                    name_strs2 = ['single domain', 'single-domain', 'single chain', 'single-chain', 'nb', ' nab']
                    name_strs3 = ['nano']
                    title_strs = ['single domain', 'single-domain', 'single chain', 'single-chain', 'nanobo', 'vhh', ' nb', ' nab', 'camel', 'lama', 'vicugna']
                    if string_contains_substr(name_H, name_strs1) \
                            or (string_contains_substr(name_H, name_strs2) and string_contains_substr(spe_H, spe_strs)) \
                            or (string_contains_substr(name_H, name_strs3) and string_contains_substr(spe_H, spe_strs) and (len(all_paired_HLs)==1 or L_chain_id=='')) \
                            or (string_contains_substr(title, title_strs) and string_contains_substr(spe_H, spe_strs) and (len(all_paired_HLs)==1 or L_chain_id=='')):
                        if is_val_within(filled_seq_len_H, min_len_vh, max_len_vh): 
                            type_H = 'VHH'
                        elif is_val_within(filled_seq_len_H, 1, min_len_vh-1): 
                            type_H = 'VHH-'
                        else: 
                            type_H = 'VHH+'
                    else:
                        if is_val_within(filled_seq_len_H, min_len_vh, max_len_vh):
                            type_H = 'VH'
                        elif is_val_within(filled_seq_len_H, min_len_fab, max_len_fab) and \
                                (string_contains_substr(name_H, ['fab', 'heavy', 'mab', 'antibo', 'immunoglo', 'igh']) 
                                        or string_contains_substr(title, ['fab', 'heavy', 'mab', 'antibo', 'immunoglo', 'igh'])):
                            type_H = 'FabH'
                        elif is_val_within(filled_seq_len_H, 1, min_len_vh-1):
                            type_H = 'VH-'
                        else:
                            type_H = 'VH+'
                elif abrsa_H == 'heavy_light':
                    if radius_H < max_radius_scfv:
                        if is_val_within(filled_seq_len_H, min_len_vhvl, max_len_vhvl):
                            type_H = 'scFv'
                        elif is_val_within(filled_seq_len_H, 1, min_len_vhvl-1):
                            type_H = 'scFv-'
                        else:
                            type_H = 'scFv+'
                    else:
                        if is_val_within(filled_seq_len_H, min_len_vhvl, max_len_vhvl):
                            type_H = 'VHVL'
                        elif is_val_within(filled_seq_len_H, 1, min_len_vhvl-1):
                            type_H = 'VHVL-'
                        else:
                            type_H = 'VHVL+'
            # infer the most possible ab type for light chain
            if L_chain_id != '':
                if abrsa_L == 'light':
                    if is_val_within(filled_seq_len_L, min_len_vl, max_len_vl):
                        type_L = 'VL'
                    elif is_val_within(filled_seq_len_L, min_len_fab, max_len_fab) and \
                            (string_contains_substr(name_L, ['fab', 'light', 'mab', 'antibo', 'immunoglo']) 
                                    or string_contains_substr(title, ['fab', 'light', 'mab', 'antibo', 'immunoglo'])):
                        type_L = 'FabL'
                    elif is_val_within(filled_seq_len_L, 1, min_len_vl-1):
                        type_L = 'VL-'
                    else:
                        type_L = 'VL+'
                elif abrsa_L == 'heavy_light':
                    if radius_L < max_radius_scfv:
                        if is_val_within(filled_seq_len_L, min_len_vhvl, max_len_vhvl):
                            type_L = 'scFv'
                        elif is_val_within(filled_seq_len_L, 1, min_len_vhvl-1):
                            type_L = 'scFv-'
                        else:
                            type_L = 'scFv+'
                    else:
                        if is_val_within(filled_seq_len_L, min_len_vhvl, max_len_vhvl):
                            type_L = 'VHVL'
                        elif is_val_within(filled_seq_len_L, 1, min_len_vhvl-1):
                            type_L = 'VHVL-'
                        else:
                            type_L = 'VHVL+'
            # update ab type for paired HL chains
            if type_H == 'VH+' and type_L == 'VL+':
                if is_val_within(filled_seq_len_H, min_len_fab, max_len_fab):
                    type_H = 'FabH'
                if is_val_within(filled_seq_len_L, min_len_fab, max_len_fab):
                    type_L = 'FabL'
            elif type_H == 'FabH' and type_L == 'VL+':
                if is_val_within(filled_seq_len_L, min_len_fab, max_len_fab):
                    type_L = 'FabL'
            elif type_H == 'VH+' and type_L == 'FabL':
                if is_val_within(filled_seq_len_H, min_len_fab, max_len_fab):
                    type_H = 'FabH'
            elif type_H == 'VHH+' and type_L == 'FabL':
                if is_val_within(filled_seq_len_H, min_len_fab, max_len_fab) and spe_H == spe_L:
                    type_H = 'FabH'
                elif (not is_val_within(filled_seq_len_H, min_len_fab, max_len_fab)) and spe_H == spe_L:
                    type_H = 'VH+'
            # add ab type to paired_HL
            if H_chain_id != '' and L_chain_id != '':
                paired_HL.append(f'{type_H}:{type_L}')
            elif H_chain_id != '' and L_chain_id == '':
                paired_HL.append(f'{type_H}')
            elif H_chain_id == '' and L_chain_id != '':
                paired_HL.append(f'{type_L}')
    
    return


def build_ab_chain_pairs_within_entity(pdbid, ent_1):
    for i in range(0, len(ent_1.get_chain_ids())-1):
        chain_id_1 = ent_1.get_chain_id(i)
        for j in range(i+1, len(ent_1.get_chain_ids())):
            chain_id_2 = ent_1.get_chain_id(j)
            if ent_1.get_abrsa_type(i) == 'heavy_light' and ent_1.get_abrsa_type(j) == 'heavy_light':
                HL_pair = [chain_id_1, chain_id_2]
                tmp_pdb = f'{pdbid}_{chain_id_1},{chain_id_2}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_1},{chain_id_2}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_1, chain_id_2, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_2} and {chain_id_1}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')
                
                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_1.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                # VHVL:VHVL
                min_inf_res_num = min_inf_res_num_vhvl_vhvl
                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_1.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])
            elif 'heavy' in ent_1.get_abrsa_type(i) and 'light' in ent_1.get_abrsa_type(j):
                HL_pair = [chain_id_1, chain_id_2]
                tmp_pdb = f'{pdbid}_{chain_id_1},{chain_id_2}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_1},{chain_id_2}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_1, chain_id_2, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_2} and {chain_id_1}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')
    
                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_1.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                len_1, len_2 = ent_1.get_real_pdb_len(i), ent_1.get_real_pdb_len(j)
                # Fab:Fab
                if (is_val_within(len_1, min_len_fab, max_len_fab) and is_val_within(len_2, min_len_fab, max_len_fab)): 
                    min_inf_res_num = min_inf_res_num_fabh_fabl
                # VH: Fab
                elif (is_val_within(len_1, min_len_vh, max_len_vh) and is_val_within(len_2, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_vh_fabl
                # Fab: VL
                elif (is_val_within(len_1, min_len_fab, max_len_fab) and is_val_within(len_2, min_len_vl, max_len_vl)):
                    min_inf_res_num = min_inf_res_num_fabh_vl
                # VH-: VL-
                elif (is_val_within(len_1, 1, min_len_vh-1) and is_val_within(len_2, 1, min_len_vl-1)):
                    min_inf_res_num = min_inf_res_num_vhm_vlm
                else:
                    min_inf_res_num = min_inf_res_num_vh_vl

                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_1.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])
            elif 'heavy' in ent_1.get_abrsa_type(j) and 'light' in ent_1.get_abrsa_type(i):
                HL_pair = [chain_id_2, chain_id_1]
                tmp_pdb = f'{pdbid}_{chain_id_2},{chain_id_1}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_2},{chain_id_1}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_2, chain_id_1, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_2} and {chain_id_1}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')
    
                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_1.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                len_1, len_2 = ent_1.get_real_pdb_len(i), ent_1.get_real_pdb_len(j)
                # Fab:Fab
                if (is_val_within(len_2, min_len_fab, max_len_fab) and is_val_within(len_1, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_fabh_fabl
                # VH:Fab
                elif (is_val_within(len_2, min_len_vh, max_len_vh) and is_val_within(len_1, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_vh_fabl
                # Fab:VL
                elif (is_val_within(len_2, min_len_fab, max_len_fab) and is_val_within(len_1, min_len_vl, max_len_vl)):
                    min_inf_res_num = min_inf_res_num_fabh_vl
                # VH-: VL-
                elif (is_val_within(len_2, 1, min_len_vh-1) and is_val_within(len_1, 1, min_len_vl-1)):
                    min_inf_res_num = min_inf_res_num_vhm_vlm
                else:
                    min_inf_res_num = min_inf_res_num_vh_vl

                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_1.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])

    return


def build_ab_chain_pairs_between_entities(pdbid, ent_1, ent_2):
    for i, chain_id_1 in enumerate(ent_1.get_chain_ids()):
        for j, chain_id_2 in enumerate(ent_2.get_chain_ids()):
            if ent_1.get_abrsa_type(i) == 'heavy_light' and ent_2.get_abrsa_type(j) == 'heavy_light':
                HL_pair = [chain_id_1, chain_id_2]
                tmp_pdb = f'{pdbid}_{chain_id_1},{chain_id_2}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_1},{chain_id_2}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_1, chain_id_2, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_1} and {chain_id_2}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')

                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_2.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                len_1, len_2 = ent_1.get_filled_pdb_len(i), ent_2.get_filled_pdb_len(j)
                # set a very high threshold for VHVL:VHVL dibody
                min_inf_res_num = min_inf_res_num_vhvl_vhvl
                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_2.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])

            elif 'heavy' in ent_1.get_abrsa_type(i) and 'light' in ent_2.get_abrsa_type(j):
                HL_pair = [chain_id_1, chain_id_2]
                tmp_pdb = f'{pdbid}_{chain_id_1},{chain_id_2}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_1},{chain_id_2}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_1, chain_id_2, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_1} and {chain_id_2}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')
                
                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_2.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                len_1, len_2 = ent_1.get_real_pdb_len(i), ent_2.get_real_pdb_len(j)
                # Fab:Fab
                if (is_val_within(len_1, min_len_fab, max_len_fab) and is_val_within(len_2, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_fabh_fabl
                # VH:Fab
                elif (is_val_within(len_1, min_len_vh, max_len_vh) and is_val_within(len_2, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_vh_fabl
                # Fab:VL
                elif (is_val_within(len_1, min_len_fab, max_len_fab) and is_val_within(len_2, min_len_vl, max_len_vl)):
                    min_inf_res_num = min_inf_res_num_fabh_vl
                # VH-: VL-
                elif (is_val_within(len_1, 1, min_len_vh-1) and is_val_within(len_2, 1, min_len_vl-1)):
                    min_inf_res_num = min_inf_res_num_vhm_vlm
                else:
                    min_inf_res_num = min_inf_res_num_vh_vl

                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_2.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])

            elif 'heavy' in ent_2.get_abrsa_type(j) and 'light' in ent_1.get_abrsa_type(i):
                HL_pair = [chain_id_2, chain_id_1]
                tmp_pdb = f'{pdbid}_{chain_id_2},{chain_id_1}.pdb'
                tmp_unidesign_out = f'{pdbid}_{chain_id_2},{chain_id_1}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', HL_pair, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, chain_id_2, chain_id_1, tmp_unidesign_out)
                inf_res_1, inf_res_2, HL_inf_res_num, mean_pair_index = parse_HL_unidesign_out(tmp_unidesign_out, chain_id_1, chain_id_2)
                if verbose:
                    print(f'HL_inf_res_num between {chain_id_2} and {chain_id_1}: {HL_inf_res_num}, mean_pair_index: {mean_pair_index}')
                
                cdr_1 = ent_1.get_cdr_res()
                cdr_2 = ent_2.get_cdr_res()
                comm_1 = list(set(cdr_1[0] + cdr_1[1]) & set(inf_res_1))
                comm_2 = list(set(cdr_2[0] + cdr_2[1]) & set(inf_res_2))
                if inf_res_1: cdr_inf_res_ratio_1 = len(comm_1)/len(inf_res_1)
                else: cdr_inf_res_ratio_1 = 0
                if inf_res_2: cdr_inf_res_ratio_2 = len(comm_2)/len(inf_res_2)
                else: cdr_inf_res_ratio_2 = 0
                if verbose:
                    print(f'cdr_inf_res_ratio_1: {cdr_inf_res_ratio_1}, cdr_inf_res_ratio_2: {cdr_inf_res_ratio_2}')
                if cdr_inf_res_ratio_1 > max_cdr_inf_res_ratio or cdr_inf_res_ratio_2 > max_cdr_inf_res_ratio:
                    continue
                
                len_1, len_2 = ent_1.get_real_pdb_len(i), ent_2.get_real_pdb_len(j)
                # Fab:Fab
                if (is_val_within(len_2, min_len_fab, max_len_fab) and is_val_within(len_1, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_fabh_fabl
                # VH:Fab
                elif (is_val_within(len_2, min_len_vh, max_len_vh) and is_val_within(len_1, min_len_fab, max_len_fab)):
                    min_inf_res_num = min_inf_res_num_vh_fabl
                # Fab:VL
                elif (is_val_within(len_2, min_len_fab, max_len_fab) and is_val_within(len_1, min_len_vl, max_len_vl)):
                    min_inf_res_num = min_inf_res_num_fabh_vl
                # VH-: VL-
                elif (is_val_within(len_2, 1, min_len_vh-1) and is_val_within(len_1, 1, min_len_vl-1)):
                    min_inf_res_num = min_inf_res_num_vhm_vlm
                else:
                    min_inf_res_num = min_inf_res_num_vh_vl
                
                if HL_inf_res_num >= min_inf_res_num:
                    ent_1.paired_chains[i].append([chain_id_2, HL_inf_res_num, mean_pair_index])
                    ent_2.paired_chains[j].append([chain_id_1, HL_inf_res_num, mean_pair_index])
    return


def build_all_ab_chain_pairs(pdbid, unidesign, pdb_entities):
    # find all possible HL pairs
    for index_1 in range(0, len(pdb_entities)):
        ent_1 = pdb_entities[index_1]
        build_ab_chain_pairs_within_entity(pdbid, ent_1)
        for index_2 in range(index_1+1, len(pdb_entities)):
            ent_2 = pdb_entities[index_2]
            build_ab_chain_pairs_between_entities(pdbid, ent_1, ent_2)
    if verbose:
        print(f'candidate HL pairings are as follows:')
        for ent in pdb_entities:
            for i, chain_id in enumerate(ent.get_chain_ids()):
                if ent.get_abrsa_type(i) == 'antigen': continue
                print(f'chain_id: {chain_id}, paired chains: {ent.paired_chains[i]}')
    
    multi_pair_exists = False
    for ent in pdb_entities:
        for i, chain_id in enumerate(ent.get_chain_ids()):
            if ent.get_abrsa_type(i) == 'antigen': continue
            if len(ent.paired_chains[i])>1:
                multi_pair_exists = True
                break
        if multi_pair_exists: break

    if multi_pair_exists:
        # determine best pairings by greedy search
        find_best_pairing_by_greedy_search(pdb_entities)
    else:
        # pair all unpaired ab chains with a dummay chain
        for ent in pdb_entities:
            for i, chain_id in enumerate(ent.get_chain_ids()):
                if ent.get_abrsa_type(i) == 'antigen': continue
                if not ent.paired_chains[i]:
                    ent.paired_chains[i].append(['', 0, 1000])

    # save and reformat best pairings
    all_paired_HLs = []
    for ent in pdb_entities:
        HLs= []
        for i, chain_id in enumerate(ent.get_chain_ids()):
            if 'heavy' in ent.get_abrsa_type(i):
                new_pair = True
                pair_id, inf_res_num, mean_pair_index = ent.paired_chains[i][0]
                HL = [chain_id, pair_id, inf_res_num, mean_pair_index]
                LH = [pair_id, chain_id, inf_res_num, mean_pair_index]
                if pair_id:
                    paired_ent, j = find_ent_by_chain_id(pdb_entities, pair_id)
                    if paired_ent.get_abrsa_type(j) == 'heavy':
                        continue
                for paired_HLs in all_paired_HLs:
                    if HL in paired_HLs or LH in paired_HLs:
                        new_pair = False
                        break
                if new_pair and (HL not in HLs) and (LH not in HLs):
                    HLs.append(HL)
            elif ent.get_abrsa_type(i) == 'light':
                new_pair = True
                pair_id, inf_res_num, mean_pair_index = ent.paired_chains[i][0]
                HL = [pair_id, chain_id, inf_res_num, mean_pair_index]
                for paired_HLs in all_paired_HLs:
                    if HL in paired_HLs:
                        new_pair = False
                        break
                if new_pair and (HL not in HLs):
                    HLs.append(HL)
        # reformat HLs
        if HLs:
            HLs_list = [[] for i in range(3)]
            for HL in HLs:
                id1, id2, num3, index4 = HL
                if id1 and id2:
                    HLs_list[0].append(HL)
                elif id1 and (not id2):
                    HLs_list[1].append(HL)
                else:
                    HLs_list[2].append(HL)
            for i in range(3):
                if HLs_list[i]:
                    all_paired_HLs.append(HLs_list[i])
    return all_paired_HLs


def is_chain_id_paired(chain_id, chain_pairs):
    for pair in chain_pairs:
        if chain_id and chain_id in pair:
            return True
    return False

def find_best_pairing_by_greedy_search(pdb_entities):
    # temporarily pair all unpaired ab chains with a dummay chain
    chain_ids, paired_chains = [], []
    for ent in pdb_entities:
        for i, chain_id in enumerate(ent.get_chain_ids()):
            if ent.get_abrsa_type(i) != 'antigen':
                ent.paired_chains[i].append(['', 0, 1000])
                chain_ids.append(chain_id)
                paired_chains.append(ent.paired_chains[i])
    if verbose:
        print(f'chain_ids: {chain_ids}, paired_chains: {paired_chains}')
    
    # greedy search
    pair_scores = dict()
    for i, chain_id in enumerate(chain_ids):
        for j, paired_chain in enumerate(paired_chains[i]):
            pair_id = paired_chain[0]
            pair_score = 3*paired_chain[1] - paired_chain[2]
            if (chain_id, pair_id) in pair_scores or (pair_id, chain_id) in pair_scores: continue
            pair_scores[(chain_id, pair_id)] = pair_score
    sorted_pair_scores = dict(sorted(pair_scores.items(), key=lambda item: item[1], reverse=True))
    if verbose:
        print(f'sorted_pair_scores: {sorted_pair_scores}')
    sorted_chain_pairs = []
    for key in sorted_pair_scores:
        sorted_chain_pairs.append([key[0], key[1], sorted_pair_scores[key]])
    greedy_pairs = []
    unused_pairs = []
    greedy_score = 0
    for i, pair_i in enumerate(sorted_chain_pairs):
        if is_chain_id_paired(pair_i[0], greedy_pairs) or is_chain_id_paired(pair_i[1], greedy_pairs):
            unused_pairs.append(pair_i)
        else:
            greedy_pairs.append(pair_i)
            greedy_score += pair_i[2]
    if verbose:
        print(f'initial greedy_score: {greedy_score}, greedy_pairs: {greedy_pairs}, unused_pairs: {unused_pairs}')

    best_pairs = copy.deepcopy(greedy_pairs)
    best_unused_pairs = copy.deepcopy(unused_pairs)
    best_score = greedy_score
    
    # heuristic search to improve greedy_score if possible
    while True:
        greedy_pair = greedy_pairs[0]
        insert_from_unused = []
        added_chain_ids = []
        for j, pair_j in enumerate(unused_pairs):
            if greedy_pair[0] and greedy_pair[0] in pair_j and greedy_pair[0] not in added_chain_ids:
                insert_from_unused.append(pair_j)
                if pair_j[0] and pair_j[0] not in added_chain_ids: added_chain_ids.append(pair_j[0])
                if pair_j[1] and pair_j[1] not in added_chain_ids: added_chain_ids.append(pair_j[1])
            if greedy_pair[1] and greedy_pair[1] in pair_j and greedy_pair[1] not in added_chain_ids:
                insert_from_unused.append(pair_j)
                if pair_j[0] and pair_j[0] not in added_chain_ids: added_chain_ids.append(pair_j[0])
                if pair_j[1] and pair_j[1] not in added_chain_ids: added_chain_ids.append(pair_j[1])
        greedy_pairs.pop(0)
    
        updated_greedy_pairs = []
        add_to_unused = []
        for pair_i in greedy_pairs:
            if is_chain_id_paired(pair_i[0], insert_from_unused) or is_chain_id_paired(pair_i[1], insert_from_unused): 
                add_to_unused.append(pair_i)
                continue
            updated_greedy_pairs.append(pair_i)
        greedy_pairs = updated_greedy_pairs
        for i in insert_from_unused:
            unused_pairs.remove(i)
            greedy_pairs.append(i)
        
        unused_pairs.append(greedy_pair)
        for i in add_to_unused:
            unused_pairs.append(i)
        greedy_score = 0
        for pair_i in greedy_pairs:
            greedy_score += pair_i[2]
        
        if verbose:
            print(f'updated greedy_score: {greedy_score}, greedy_pairs: {greedy_pairs}, unused_pairs: {unused_pairs}')
        
        if greedy_score > best_score:
            best_score = greedy_score
            best_pairs = copy.deepcopy(greedy_pairs)
            best_unused_pairs = copy.deepcopy(unused_pairs)
        else:
            break
    
    if verbose:   
        print(f'best_score: {best_score}, best_pairs: {best_pairs}, best_unused_pairs: {best_unused_pairs}')

    best_pairs_id_only = []
    for pair in best_pairs:
        best_pairs_id_only.append([pair[0], pair[1]])


    # update ent.paired_chains
    for ent in pdb_entities:
        for i, chain_id in enumerate(ent.get_chain_ids()):
            if ent.get_abrsa_type(i) == 'antigen': continue
            for paired_chain in ent.paired_chains[i]:
                pair_id = paired_chain[0]
                if [chain_id, pair_id] in best_pairs_id_only or [pair_id, chain_id] in best_pairs_id_only:
                    ent.paired_chains[i] = [paired_chain]
            if verbose:
                print(f'chain_id: {chain_id}, paired_chains: {ent.paired_chains[i]}')
    return


def build_ab_ag_pairs(unidesign, pdbid, paired_HLs, pdb_entities, all_paired_ab_ags):
    paired_ab_ags = []
    for paired_HL in paired_HLs:
        H_chain_id, L_chain_id, HL_inf_res_num, ab_type = paired_HL[0], paired_HL[1], paired_HL[2], paired_HL[3]
        ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
        ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
        paired_ags = []
        # do not pair single FabH or FabL with Ag molecules
        if ab_type == 'FabH' or ab_type == 'FabL': 
            paired_ags.append(['', 0, 0, 0.00])
            paired_ab_ags.append([paired_HL, paired_ags])
            continue
        for ent_ag in pdb_entities:
            # disable the situation where an Ab chain may be an Ag chain because crystal packing may introduce fake Ab-Ag interactions
            is_ag = True
            for abrsa_type in ent_ag.get_abrsa_types():
                if 'antigen' not in abrsa_type:
                    is_ag = False
                    break
            if is_ag == False: continue
            #if ent_ag == ent_H or ent_ag == ent_L: continue
            for ag_chain_id in ent_ag.get_chain_ids():
                if ag_chain_id in [H_chain_id, L_chain_id]: continue
                select_chains = [H_chain_id, L_chain_id, ag_chain_id]
                tmp_pdb = f'{pdbid}_{H_chain_id},{L_chain_id}_{ag_chain_id}.pdb'
                tmp_unidesign_out = f'{pdbid}_{H_chain_id},{L_chain_id}_{ag_chain_id}.unidesign'
                extract_pdb_chains_to_file(f'{pdbid}.pdb', select_chains, tmp_pdb)
                run_unidesign_find_interface_residues(unidesign, tmp_pdb, H_chain_id+L_chain_id, ag_chain_id, tmp_unidesign_out)
                H_inf_res, L_inf_res, ag_inf_res, inf_res_num = parse_ab_ag_unidesign_out(tmp_unidesign_out, H_chain_id, L_chain_id, ag_chain_id)
                if verbose: print(f'pair H, L, and ag chains: {select_chains}, ab-ag_inf_res_num: {inf_res_num}')
                # check if the ab-ag interface is big enough with a threshold of the number of interface residues (default: 10)
                if inf_res_num >= 10:
                    if verbose: print(f'H_inf_res: {H_inf_res}, L_inf_res: {L_inf_res}, ag_inf_res: {ag_inf_res}')
                    hcdr, lcdr = [[],[]], [[], []]
                    ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
                    if ent_H: hcdr = ent_H.get_cdr_res()
                    ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
                    if ent_L: lcdr = ent_L.get_cdr_res()
                    HL_cdr_all = list(set(hcdr[0] + lcdr[0] + hcdr[1] + lcdr[1]))
                    HL_cdr_all = [int(i) for i in HL_cdr_all]
                    HL_cdr_all.sort()
                    if verbose:
                        print(f'HL_cdr_all: {HL_cdr_all}')                        
                    # compute the number of CDR residues on the interface
                    HL_comm = list(set(HL_cdr_all) & set(H_inf_res+L_inf_res))
                    HL_comm.sort()
                    if verbose: print(f'HL_comm: {HL_comm}')
                    cdr_inf_res_num = len(HL_comm)
                    cdr_inf_res_ratio = len(HL_comm)/len(H_inf_res+L_inf_res)
                    if verbose:
                        print(f'cdr_inf_res_num: {cdr_inf_res_num}, cdr_inf_res_ratio: {cdr_inf_res_ratio}')
                    if cdr_inf_res_ratio >= 0.25 and cdr_inf_res_num >= 5: 
                        paired_ags.append([ag_chain_id, inf_res_num, cdr_inf_res_num, round(cdr_inf_res_ratio, 2)])
        if paired_ags:
            paired_ab_ags.append([paired_HL, paired_ags])
        else:
            paired_ags.append(['', 0, 0, 0.00])
            paired_ab_ags.append([paired_HL, paired_ags])
    
    if paired_ab_ags:
        all_paired_ab_ags.append(paired_ab_ags)
    
    return


def parse_HL_unidesign_out(ud_out, H_chain_id, L_chain_id):
    with open(ud_out, 'r') as fp: content = fp.read()
    p1 = re.compile(r'(?<=Interface residues: \n)(.*?)\nTime spent:', re.DOTALL)
    interface = re.findall(p1, content)
    H_inf, L_inf = [], []
    if interface:
        pH = re.compile(r'{}\s+([0-9A-Za-z]+) [A-Z]'.format(H_chain_id), re.DOTALL)
        H_inf = pH.findall(interface[0])
        H_inf2 = [int(i) for i in H_inf]
        pL = re.compile(r'{}\s+([0-9A-Za-z]+) [A-Z]'.format(L_chain_id), re.DOTALL)
        L_inf = pL.findall(interface[0])
        L_inf = [int(i) for i in L_inf]
        return H_inf, L_inf, len(H_inf)+len(L_inf), int((np.mean(np.array(H_inf2))+np.mean(np.array(L_inf)))/2)
    else:
        return [], [], 0, 1000


def parse_ab_ag_unidesign_out(ud_out, H_chain_id, L_chain_id, ag_chain_id):
    with open(ud_out, 'r') as f: 
        content = f.read()
    p1 = re.compile(r'(?<=Interface residues: \n)(.*?)\nTime spent:', re.DOTALL)
    interface = re.findall(p1, content)
    H_inf, L_inf, ag_inf = [], [], []
    if interface:
        if H_chain_id:
            pH = re.compile(r'{}\s+([0-9A-Za-z]+) [A-Z]'.format(H_chain_id), re.DOTALL)
            H_inf = pH.findall(interface[0])
            H_inf = [int(i) for i in H_inf]
        if L_chain_id:
            pL = re.compile(r'{}\s+([0-9A-Za-z]+) [A-Z]'.format(L_chain_id), re.DOTALL)
            L_inf = pL.findall(interface[0])
            L_inf = [int(i) for i in L_inf]
        if ag_chain_id:
            pAg = re.compile(r'{}\s+([0-9A-Za-z]+) [A-Z]'.format(ag_chain_id), re.DOTALL)
            ag_inf = pAg.findall(interface[0])
            ag_inf = [int(i) for i in ag_inf]
        return H_inf, L_inf, ag_inf, len(H_inf)+len(L_inf)+len(ag_inf)
    else:
        return [], [], [], 0


def find_ent_by_chain_id(pdb_entities, chain_id):
    for ent in pdb_entities:
        for i, saved_id in enumerate(ent.get_chain_ids()):
            if chain_id == saved_id:
                return ent, i
    return None, -1


def compute_cdr_residues(ent, pdbid, pdb_model, abrsa_pdb):
    for index, chain_id in enumerate(ent.get_chain_ids()):
        if ent.get_chain_type(index) != 'protein' and ent.get_chain_type(index) != 'peptide' and \
                'heavy' not in ent.get_abrsa_type(index) and 'light' not in ent.get_abrsa_type(index): 
            continue
        tmp_abrsa_pdb_out = f'{pdbid}_{chain_id}.abrsapdb'
        run_abrsa_pdb(abrsa_pdb, f'{pdbid}_prot_{chain_id}.pdb', tmp_abrsa_pdb_out)
        hcdr_res, lcdr_res = parse_abrsapdb_out(tmp_abrsa_pdb_out)
        ent.set_cdr_res(hcdr_res, lcdr_res)
    return


def find_representative_ab_ag_pair_index(paired_ab_ags, pdb_entities):
    rep_ndx = 0
    rep_score = -10000
    for index, paired_ab_ag in enumerate(paired_ab_ags):
        score = 0
        paired_HL, paired_ags = paired_ab_ag[0], paired_ab_ag[1]
        H_chain_id, L_chain_id, HL_inf_res_num = paired_HL[0], paired_HL[1], paired_HL[2]

        # reward H_pdb_seq_len and L_pdb_seq_len
        if H_chain_id:
            ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
            score += ent_H.get_real_pdb_len(index_H)
        if L_chain_id:
            ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
            score += ent_L.get_real_pdb_len(index_L)
        
        # reward HL_inf_res_num
        score += paired_HL[2]
        
        # reward paired ab-ag
        if paired_ags[0][0]:
            score += 100

        # but penalize multi-ag chains
        if len(paired_ags)>1:
            score -= 100*len(paired_ags)
        
        # rewad ab_ag_inf_res_num, cdr_inf_res_num, and cdr_inf_res_ratio
        for paired_ag in paired_ags:
            score += (paired_ag[1] + paired_ag[2] + 10*paired_ag[3])
        
        if score > rep_score:
            rep_ndx = index
            rep_score = score
        
        if verbose:
            print(f'score: {score}, rep_ndx: {rep_ndx}')
        
    return rep_ndx


def is_val_within(val, lower, upper):
    if val >= lower and val <= upper:
        return True
    else:
        return False


def parse_auth_chains(chain_str):
    chain_ids = []
    p1 = re.compile(r'(?<=auth )\w+')
    if chain_str.startswith('Chain '):
        auth_asym_ids = p1.findall(chain_str)
        if auth_asym_ids:
            chain_ids.append(auth_asym_ids[0])
        else:
            p2 = re.compile(r'(?<=Chain )\w+')
            label_asym_ids = p2.findall(chain_str)
            chain_ids.append(label_asym_ids[0])
    elif chain_str.startswith('Chains '):
        chains = chain_str.replace('Chains ', '').split(', ')
        for chain in chains:
            auth_asym_ids = p1.findall(chain)
            if auth_asym_ids:
                chain_ids.append(auth_asym_ids[0])
            else:
                chain_ids.append(chain)
    return chain_ids


def parse_label_chains(chain_str):
    chain_ids = []
    if chain_str.startswith('Chain '):
        p1 = re.compile(r'(?<=Chain )\w+')
        label_asym_ids = p1.findall(chain_str)
        chain_ids.append(label_asym_ids[0])
    elif chain_str.startswith('Chains '):
        p2 = re.compile(r'(\w+)[\[\|]')
        chains = chain_str.replace('Chains ', '').split(', ')
        for chain in chains:
            label_asym_ids = p2.findall(chain)
            if label_asym_ids:
                chain_ids.append(label_asym_ids[0])
            else:
                chain_ids.append(chain)
    return chain_ids


def merge_identical_chain_entities(pdb_entities):
    new_entities = []
    flags = [True for i in range(len(pdb_entities))]
    for index_1 in range(0, len(pdb_entities)):
        seq_1 = pdb_entities[index_1].get_chain_fas_seq()
        for index_2 in range(index_1+1, len(pdb_entities)):
            seq_2 = pdb_entities[index_2].get_chain_fas_seq()
            aligner = Align.PairwiseAligner()
            alignments = aligner.align(seq_1, seq_2)
            aligned_seq1, aligned_seq2 = alignments[0]
            matches = sum(res1 == res2 for res1, res2 in zip(aligned_seq1, aligned_seq2))
            seqid = matches / max(len(seq_1), len(seq_2))
            if(float(seqid) > 0.99):
                pdb_entities[index_1].merge_chain_entities(pdb_entities[index_2])
                flags[index_2] = False
        if flags[index_1] == True:
            new_entities.append(pdb_entities[index_1])
    return new_entities


def get_pdb_entities_from_fasta(pdbid, pdb_fasta_file):
    header_list, seq_dict = read_fasta(pdb_fasta_file)

    global asym_id_type
    if asym_id_type == 'auth':
        asym_id_type = 'auth_asym_id'
    elif asym_id_type == 'label':
        asym_id_type = 'label_asym_id'
    # auto can be either auth or label, depending on the max length of chain id in the official fasta file
    elif asym_id_type == 'auto':
        max_chain_id_len_1 = 1
        max_chain_id_len_2 = 1
        for header in header_list:
            pdbid_ndx, chain_str, chain_name, species = parse_pdb_fasta_header(header)
            chain_ids1 = parse_auth_chains(chain_str)
            for id1 in chain_ids1:
                if max_chain_id_len_1 < len(id1):
                    max_chain_id_len_1 = len(id1)
            chain_ids2 = parse_label_chains(chain_str)
            for id2 in chain_ids2:
                if max_chain_id_len_2 < len(id2):
                    max_chain_id_len_2 = len(id2)
        if max_chain_id_len_1 == 1:
            asym_id_type = 'auth_asym_id'
        elif max_chain_id_len_2 == 1:
            asym_id_type = 'label_asym_id'
        elif max_chain_id_len_1 <= max_chain_id_len_2:
            asym_id_type = 'auth_asym_id'
        elif max_chain_id_len_2 < max_chain_id_len_1:
            asym_id_type = 'label_asym_id'
    print(f'use {asym_id_type} for chain processing')
    
    pdb_entities = []
    for header in header_list:
        pdbid_ndx, chain_str, chain_name, species = parse_pdb_fasta_header(header)
        if asym_id_type == 'auth_asym_id':
            chain_ids = parse_auth_chains(chain_str)
        else:
            chain_ids = parse_label_chains(chain_str)
        fas_seq = seq_dict[header]
        ent = Entity(chain_ids, chain_name, fas_seq, species)
        pdb_entities.append(ent)
    #pdb_entities = merge_identical_chain_entities(pdb_entities)

    # initialize the lists of each ent object
    for ent in pdb_entities:
        ent.init_chain_types()
        ent.init_real_pdb_seqs()
        ent.init_filled_pdb_seqs()
        ent.init_abrsa_types()
        ent.init_paired_chains()
        ent.init_vgene_subgroups()
        ent.init_mean_radii()

    return pdb_entities


def remove_poor_residues_and_chains(model):
    residues_to_remove, chains_to_remove = [], []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ' or residue.get_resname() == 'UNK':
                residues_to_remove.append((chain.id, residue.id))
            else:
                ca_found, p_found, c1p_found = False, False, False
                for atom in residue.get_atoms():
                    if atom.get_name() == 'CA':
                        ca_found = True
                    # atoms P and C1' are for DNA or RNA
                    elif atom.get_name() == 'P':
                        p_found = True
                    elif atom.get_name() == "C1'":
                        c1p_found = True
                if ca_found == False and p_found == False and c1p_found == False:
                    residues_to_remove.append((chain.id, residue.id))
            # some structures have very large values which will violate the PDB format
            # reset these values for the sake of saving pdb models for Pulchra and FASPR processing
            for atom in residue.get_atoms():
                atom.set_occupancy(1.00)
                atom.set_bfactor(0.00)
                atom.set_altloc(' ')
    for rtr in residues_to_remove:
        model[rtr[0]].detach_child(rtr[1])
    for chain in model:
        if len(chain) == 0:
            chains_to_remove.append(chain.id)
    for chain in chains_to_remove:
        model.detach_child(chain)
    return model


def get_pdb_model_from_cif(pdb_cif_path, pdbid):
    auth_chains = True
    if asym_id_type != 'auth_asym_id':
        auth_chains = False
    parser = PDB.MMCIFParser(auth_chains=auth_chains)
    if os.path.exists(f'{pdb_cif_path}/{pdbid}.cif'):
        with open(f'{pdb_cif_path}/{pdbid}.cif', 'r') as f:
            structure = parser.get_structure(pdbid, f)
            model = structure[0]
            model = remove_poor_residues_and_chains(model)
            return model
    elif os.path.exists(f'{pdb_cif_path}/{pdbid[1:3]}/{pdbid}.cif.gz'):
        with gzip.open(f'{pdb_cif_path}/{pdbid[1:3]}/{pdbid}.cif.gz', 'rt') as f:
            structure = parser.get_structure(pdbid, f)
            model = structure[0]
            model = remove_poor_residues_and_chains(model)
            return model


def get_pdb_model_from_pdb(pdbid, pdb_file):
    parser = PDB.PDBParser()
    with open(pdb_file, 'r') as f:
        structure = parser.get_structure(pdbid, f)
        model = structure[0]
        return model


def parse_abrsa_out(abrsa_out):
    with open(abrsa_out, 'r') as f:
        H_chain_found, L_chain_found = False, False
        ab_res = []
        counter = 1
        for line in f.readlines():
            strs = line.rstrip('\n').split()
            if line.startswith('-_EXT'):
                counter += len(strs[2])
            elif line.startswith('H_FR'):
                H_chain_found = True
                for i in range(counter, counter+len(strs[2])): ab_res.append(i)
                counter += len(strs[2])
            elif line.startswith('L_FR'):
                L_chain_found = True
                for i in range(counter, counter+len(strs[2])): ab_res.append(i)
                counter += len(strs[2])
            elif line.startswith('H_CDR'):
                H_chain_found = True
                for i in range(counter, counter+len(strs[1])): ab_res.append(i)
                counter += len(strs[1])
            elif line.startswith('L_CDR'):
                L_chain_found = True
                for i in range(counter, counter+len(strs[1])): ab_res.append(i)
                counter += len(strs[1])
        if H_chain_found == True and L_chain_found == True:
            return 'heavy_light', ab_res
        elif H_chain_found == True and L_chain_found == False:
            return 'heavy', ab_res
        elif H_chain_found == False and L_chain_found == True:
            return 'light', ab_res
        elif H_chain_found == False and L_chain_found == False:
            return 'N.A.', ab_res


def parse_abrsapdb_out(abrsa_out):
    with open(abrsa_out, 'r') as f:
        content = f.read()
        p_hcdr1 = re.compile(r'\n@H_CDR1:,([0-9,]+)')
        p_hcdr2 = re.compile(r'\n@H_CDR2:,([0-9,]+)')
        p_hcdr3 = re.compile(r'\n@H_CDR3:,([0-9,]+)')
        p_lcdr1 = re.compile(r'\n@L_CDR1:,([0-9,]+)')
        p_lcdr2 = re.compile(r'\n@L_CDR2:,([0-9,]+)')
        p_lcdr3 = re.compile(r'\n@L_CDR3:,([0-9,]+)')
        hcdr1_pos, hcdr2_pos, hcdr3_pos, lcdr1_pos, lcdr2_pos, lcdr3_pos = [], [], [], [], [], []
        hcdr1_pos = p_hcdr1.findall(content)
        hcdr1_pos = hcdr1_pos[0].split(',') if hcdr1_pos else []
        hcdr2_pos = p_hcdr2.findall(content)
        hcdr2_pos = hcdr2_pos[0].split(',') if hcdr2_pos else []
        hcdr3_pos = p_hcdr3.findall(content)
        hcdr3_pos = hcdr3_pos[0].split(',') if hcdr3_pos else []
        lcdr1_pos = p_lcdr1.findall(content)
        lcdr1_pos = lcdr1_pos[0].split(',') if lcdr1_pos else []
        lcdr2_pos = p_lcdr2.findall(content)
        lcdr2_pos = lcdr2_pos[0].split(',') if lcdr2_pos else []
        lcdr3_pos = p_lcdr3.findall(content)
        lcdr3_pos = lcdr3_pos[0].split(',') if lcdr3_pos else []
        return hcdr1_pos + hcdr2_pos + hcdr3_pos, lcdr1_pos + lcdr2_pos + lcdr3_pos


def run_abrsa(abrsa, fasta, abrsa_out):
    print('run abrsa, try chothia scheme (default)')
    os.system(f'{abrsa} -c -i {fasta} > {abrsa_out}')
    abrsa_type, ab_res = parse_abrsa_out(f'{abrsa_out}')
    if abrsa_type != 'N.A.':
        return abrsa_type, ab_res
    else:
        print(f'chothia scheme failed, try imgt scheme')
        os.system(f'{abrsa} -g -i {fasta} > {abrsa_out}')
        abrsa_type, ab_res = parse_abrsa_out(f'{abrsa_out}')
        if abrsa_type != 'N.A.':
            return abrsa_type, ab_res
        else:
            print(f'both chothia and imgt schemes failed, try kabat scheme')
            os.system(f'{abrsa} -k -i {fasta} > {abrsa_out}')
            abrsa_type, ab_res = parse_abrsa_out(f'{abrsa_out}')
            return abrsa_type, ab_res


def run_abalign(abalign, fasta, pdbid, chain_id, chain_type = 'heavy'):
    if chain_type == 'heavy' or chain_type == 'light': ctype = chain_type[0]
    #print('run abalign, try chothia scheme')
    os.system(f'{abalign} -g -s -vn 1 -i {fasta} -a{ctype} {pdbid}_{chain_id}.msa -v {pdbid}_{chain_id}.vgene > {pdbid}_{chain_id}.abalign')
    subgroup = parse_abalign_out(f'{pdbid}_{chain_id}.vgene')
    if subgroup != 'N.A.':
        return subgroup
    else:
        #print('chothia scheme failed, try imgt scheme')
        os.system(f'{abalign} -c -s -vn 1 -i {fasta} -a{ctype} {pdbid}_{chain_id}.msa -v {pdbid}_{chain_id}.vgene > {pdbid}_{chain_id}.abalign')
        subgroup = parse_abalign_out(f'{pdbid}_{chain_id}.vgene')
        if subgroup != 'N.A.':
            return subgroup
        else:
            #print('both chothia and imgt schemes failed, try kabat scheme')
            os.system(f'{abalign} -k -s -vn 1 -i {fasta} -a{ctype} {pdbid}_{chain_id}.msa -v {pdbid}_{chain_id}.vgene > {pdbid}_{chain_id}.abalign')
            subgroup = parse_abalign_out(f'{pdbid}_{chain_id}.vgene')
            return subgroup


def run_tmalign(tmalign, query_pdb, ref_pdb, tmalign_out):
    os.system(f'{tmalign} {query_pdb} {ref_pdb} -fast > {tmalign_out}')
    with open(tmalign_out, 'r') as f:
        for line in f.readlines():
            if line.startswith('TM-score') and 'Chain_2' in line:
                strs = line.split()
                return float(strs[1])
    return 0.0


def parse_abalign_out(vgene_file):
    if not os.path.exists(vgene_file):
        return 'N.A.'
    with open(vgene_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('Seq_ID'):
                continue
            elif line:
                strs = line.strip().split(',')
                p = re.compile(r'([A-Z]+[0-9]+)[A-Z-]')
                subgroups = re.findall(p, strs[1])
                if subgroups:
                    subgroup = subgroups[0]
                else:
                    subgroup = 'N.A.'
                return subgroup
    return 'N.A.'


def calc_pairwise_seq_identity(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    aligned_seq1, aligned_seq2 = alignments[0]
    matches = sum(res1 == res2 for res1, res2 in zip(aligned_seq1, aligned_seq2))
    seqid1 = matches / len(seq1)
    seqid2 = matches / len(seq2)
    return seqid1, seqid2

def align_pdb_seq_to_fas_seq(pdb_seq, fas_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.target_left_open_gap_score = 0
    aligner.target_left_extend_gap_score = 0
    aligner.target_internal_open_gap_score = -1
    aligner.target_internal_extend_gap_score = 0
    aligner.target_right_open_gap_score = 0
    aligner.target_right_extend_gap_score = 0
    aligner.query_left_open_gap_score = -2
    aligner.query_left_extend_gap_score = -2
    aligner.query_internal_open_gap_score = -2
    aligner.query_internal_extend_gap_score = -2
    aligner.query_right_open_gap_score = -2
    aligner.query_right_extend_gap_score = -2
    #print(f"aligner mode: {aligner.mode}, mismatch_score: {aligner.mismatch_score}, match_score: {aligner.match_score}, open_gap_score: {aligner.open_gap_score}, extend_gap_score: {aligner.extend_gap_score}, matrix: {aligner.substitution_matrix}")
    alignments = aligner.align(pdb_seq, fas_seq)
    aligned_pdb_seq, aligned_fas_seq = alignments[0]
    res_pos = []
    index = 1
    for res1, res2 in zip(aligned_pdb_seq, aligned_fas_seq):
        if res1 == res2:
            res_pos.append(index)
        index = index+1
    return res_pos



def write_ab_info(f, pdb_entities, chain_mapping, H_chain_id, L_chain_id, HL_inf_res_num, ab_type):
    if H_chain_id != '' and L_chain_id != '':
        ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
        ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
        sgs_H, sgs_L = ent_H.get_vgene_subgroup(index_H), ent_L.get_vgene_subgroup(index_L)
        H_chain_name, L_chain_name = ent_H.get_name(), ent_L.get_name()
        H_species, L_species = ent_H.get_species(), ent_L.get_species()
        H_pdb_seq_len, L_pdb_seq_len = ent_H.get_real_pdb_len(index_H), ent_L.get_real_pdb_len(index_L)
        H_filled_seq_len, L_filled_seq_len = ent_H.get_filled_pdb_len(index_H), ent_L.get_filled_pdb_len(index_L)
        H_filled_pdb_seq, L_filled_pdb_seq = ent_H.get_filled_pdb_seq(index_H), ent_L.get_filled_pdb_seq(index_L)
        H_fas_seq, L_fas_seq = ent_H.get_fas_seq(), ent_L.get_fas_seq()
        H_fas_seq_len, L_fas_seq_len = ent_H.get_fas_seq_len(), ent_L.get_fas_seq_len()
        H_radius, L_radius = ent_H.get_mean_radius(index_H), ent_L.get_mean_radius(index_L)
        f.write('\t'.join([ab_type, sgs_H[0], sgs_L[1],
            '"'+chain_mapping[H_chain_id]+'"', '"'+chain_mapping[L_chain_id]+'"', H_fas_seq, L_fas_seq, H_filled_pdb_seq, L_filled_pdb_seq, str(H_radius), str(L_radius), 
            str(H_fas_seq_len), str(L_fas_seq_len), str(H_pdb_seq_len), str(L_pdb_seq_len), str(H_filled_seq_len), str(L_filled_seq_len), str(HL_inf_res_num), 
            '"'+H_chain_name+'"', '"'+L_chain_name+'"', '"'+H_species+'"', '"'+L_species+'"', '']))
            
    elif H_chain_id != '' and L_chain_id == '':
        ent_H, index_H = find_ent_by_chain_id(pdb_entities, H_chain_id)
        sgs_H = ent_H.get_vgene_subgroup(index_H)
        H_chain_name = ent_H.get_name()
        H_species = ent_H.get_species()
        H_pdb_seq_len = ent_H.get_real_pdb_len(index_H)
        H_filled_seq_len = ent_H.get_filled_pdb_len(index_H)
        H_filled_pdb_seq = ent_H.get_filled_pdb_seq(index_H)
        H_fas_seq = ent_H.get_fas_seq()
        H_fas_seq_len = ent_H.get_fas_seq_len()
        H_radius = ent_H.get_mean_radius(index_H)
        f.write('\t'.join([ab_type, sgs_H[0], sgs_H[1],
            '"'+chain_mapping[H_chain_id]+'"', 'N.A.', H_fas_seq, 'N.A.', H_filled_pdb_seq, 'N.A.', str(H_radius), '0.00', 
            str(H_fas_seq_len), '0', str(H_pdb_seq_len), '0', str(H_filled_seq_len), '0', str(HL_inf_res_num),
            '"'+H_chain_name+'"', 'N.A.', '"'+H_species+'"', 'N.A.', '']))
            
    elif H_chain_id == '' and L_chain_id != '':
        ent_L, index_L = find_ent_by_chain_id(pdb_entities, L_chain_id)
        sgs_L = ent_L.get_vgene_subgroup(index_L)
        L_chain_name = ent_L.get_name()
        L_species = ent_L.get_species()
        L_pdb_seq_len = ent_L.get_real_pdb_len(index_L)
        L_filled_seq_len = ent_L.get_filled_pdb_len(index_L)
        L_filled_pdb_seq = ent_L.get_filled_pdb_seq(index_L)
        L_fas_seq = ent_L.get_fas_seq()
        L_fas_seq_len = ent_L.get_fas_seq_len()
        L_radius = ent_L.get_mean_radius(index_L)
        f.write('\t'.join([ab_type, sgs_L[0], sgs_L[1], 
            'N.A.', '"'+chain_mapping[L_chain_id]+'"', 'N.A.', L_fas_seq, 'N.A.', L_filled_pdb_seq, '0.00', str(L_radius), 
            '0', str(L_fas_seq_len), '0', str(L_pdb_seq_len), '0', str(L_filled_seq_len), str(HL_inf_res_num),
            'N.A.', '"'+L_chain_name+'"', 'N.A.', '"'+L_species+'"', '']))
    return


def write_ag_info(f, pdb_entities, chain_mapping, paired_ags):
    chain_mol_str, chain_spe_str, chain_type_str, chain_id_str = 'N.A.', 'N.A.', 'N.A.', 'N.A.'
    if paired_ags[0][0]:
        chain_id_str = ';'.join([chain_mapping[paired_ags[i][0]] for i in range(0, len(paired_ags))])
        mol_names, ag_species, chain_types = [], [], []
        for i in range(0, len(paired_ags)):
            ent_ag, index_ag = find_ent_by_chain_id(pdb_entities, paired_ags[i][0])
            mol_names.append(ent_ag.get_name())
            ag_species.append(ent_ag.get_species())
            chain_types.append(ent_ag.get_chain_type(index_ag))
        chain_mol_str = ';'.join(mol_names)
        chain_spe_str = ';'.join(ag_species)
        chain_type_str = ';'.join(chain_types)
    
    ab_ag_inf_res_num_str = ';'.join([str(paired_ags[i][1]) for i in range(0, len(paired_ags))])
    cdr_inf_res_num_str = ';'.join([str(paired_ags[i][2]) for i in range(0, len(paired_ags))])
    cdr_inf_res_ratio_str = ';'.join([str(paired_ags[i][3]) for i in range(0, len(paired_ags))])

    if chain_id_str == 'N.A.':
        f.write('\t'.join([chain_id_str, chain_type_str, chain_mol_str, chain_spe_str, 
            ab_ag_inf_res_num_str, cdr_inf_res_num_str, cdr_inf_res_ratio_str+'\n']))
    else:
        f.write('\t'.join(['"'+chain_id_str+'"', '"'+chain_type_str+'"', '"'+chain_mol_str+'"', '"'+chain_spe_str+'"', 
            ab_ag_inf_res_num_str, cdr_inf_res_num_str, cdr_inf_res_ratio_str+'\n']))
    return

def write_paired_ab_ag_ids(f, paired_ab_ags, chain_mapping, rep_ndx):
    for i, paired_ab_ag in enumerate(paired_ab_ags):
        paired_HL, paired_ags = paired_ab_ag[0], paired_ab_ag[1]
        H_chain_id, L_chain_id, HL_inf_res_num = paired_HL[0], paired_HL[1], paired_HL[2]
        H_chain_id_str, L_chain_id_str, ag_chain_id_str = 'N.A.', 'N.A.', 'N.A.'
        if H_chain_id:
            H_chain_id_str = chain_mapping[H_chain_id]
        if L_chain_id:
            L_chain_id_str = chain_mapping[L_chain_id]
        ag_chain_id_str = 'N.A.'
        if paired_ags[0][0]:
            ag_chain_id_str = ';'.join([chain_mapping[paired_ags[i][0]] for i in range(0, len(paired_ags))])
        if i==0:
            if i==rep_ndx:
                f.write(','.join(['rep', H_chain_id_str, L_chain_id_str, ag_chain_id_str]))
            else:
                f.write(','.join(['nonrep', H_chain_id_str, L_chain_id_str, ag_chain_id_str]))
        else:
            if i==rep_ndx:
                f.write(','.join(['\t'+'rep', H_chain_id_str, L_chain_id_str, ag_chain_id_str]))
            else:
                f.write(','.join(['\t'+'nonrep', H_chain_id_str, L_chain_id_str, ag_chain_id_str]))
    f.write('\n')
    return




def write_pdb_info(f, pdbid, mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title, asym_id_type):
    f.write('\t'.join(['"'+pdbid+'"', '"'+title+'"', mut_status, '"'+class_+'"', deposit_date, release_date, 
        '"'+method+'"', resolution, r_free, r_work, pmid, doi, asym_id_type, '']))
    return


def clean_all_temp_files():
    os.system(f'rm {pdbid}*.pdb {pdbid}*.fasta {pdbid}*.abrsa* {pdbid}*.abalign {pdbid}*.msa* {pdbid}*.vgene* {pdbid}*.unidesign {pdbid}*.tmalign 2>/dev/null')


#################################
# Main ab-ag interaction parser
#################################
def parse_ab_ag_interaction(pdbid, pdb_fasta_path, pdb_cif_path):
    print(f'parsing pdb entry {pdbid}')
    l2code = pdbid[1:3]
    pdb_fasta_file = f'{pdb_fasta_path}/{l2code}/{pdbid}.fasta'
    if not os.path.exists(pdb_fasta_file):
        print(f'done parsing (no fasta available for {pdbid})')
        exit(0)

    if not os.path.exists(f'{pdb_cif_path}/{l2code}/{pdbid}.cif.gz'):
        print(f'done parsing (no mmCIF available for {pdbid})')
        exit(0)

    # create pdb_entities from the official pdbid.fasta file
    pdb_entities = get_pdb_entities_from_fasta(pdbid, pdb_fasta_file)
    
    all_chain_abrsa_types = []
    total_chain_num = 0
    for ent in pdb_entities:
        chain_id = ent.get_chain_id(0)
        tmp_fasta = f'{pdbid}_{chain_id}.fasta'
        with open(tmp_fasta, 'w') as f:
            f.write(f'>{chain_id}\n{ent.get_fas_seq()}\n')
        tmp_abrsa_out = f'{pdbid}_{chain_id}.abrsa'
        abrsa_type, ab_res = run_abrsa(abrsa, tmp_fasta, tmp_abrsa_out)
        ent.set_fas_ab_res(ab_res)
        all_chain_abrsa_types.append(abrsa_type)
        total_chain_num += len(ent.get_chain_ids())
        # temporarily set the AbRSA type
        for i, chain_id in enumerate(ent.get_chain_ids()):
            ent.set_abrsa_type(i, abrsa_type)
        if verbose:
            print(f'chain_ids: {ent.get_chain_ids()}, fas_seq: {ent.get_fas_seq()}, fas_seq_len: {len(ent.get_fas_seq())}, abrsa_type (by fasta): {abrsa_type}')
    
    if not ('heavy' in all_chain_abrsa_types or 'heavy_light' in all_chain_abrsa_types or 'light' in all_chain_abrsa_types):
        clean_all_temp_files()
        print(f'done parsing (no ab chain found in {pdbid})')
        exit(0)

    single_letters = list(string.printable)
    chars_to_be_removed = [' ', '.', '`', ':', '<', '>', '"', '/', '\\', '|', '?', '*', '\x0c', '\n', '\t', '\r', '\x0b']
    for char in chars_to_be_removed:
        single_letters.remove(char)
    if total_chain_num >= len(single_letters):
        clean_all_temp_files()
        print(f'skip parsing (too many, {total_chain_num}) chains in the structure')
        exit(0)
    
    # create pdb_model from the official pdbid.cif file
    pdb_model = get_pdb_model_from_cif(pdb_cif_path, pdbid)

    # create one-to-one mapping for chain ids with >=2 letters
    pdb_model_chain_ids = []
    for chain in pdb_model: 
        pdb_model_chain_ids.append(chain.get_id())

    chain_mapping = dict()
    for id_ in pdb_model_chain_ids:
        # for single-letter chain, map it to itself
        if len(id_) == 1:
            chain_mapping[id_] = id_
            single_letters.remove(id_)
    for id_ in pdb_model_chain_ids:
        # for chain with >=2 letters, map it to one of the remaining single letters in the single_letters list
        if len(id_) > 1:
            if len(single_letters) > 0:
                chain_mapping[id_] = single_letters.pop(0)
    if verbose:
        print('chain mapping:', chain_mapping)

    # update pdb_entities and do chain name mapping
    new_entities = []
    for ent in pdb_entities:
        if ent.get_chain_id(0) in pdb_model_chain_ids:
            new_ids = ent.get_chain_ids()
            for i, new_id in enumerate(new_ids):
                new_ids[i] = chain_mapping[new_id]
            ent.set_chain_ids(new_ids)
            new_entities.append(ent)
    pdb_entities = new_entities

    # change chain id in pdb_model correspondingly
    for chain in pdb_model:
        if len(chain.get_id())>1:
            temp_id = chain.get_id()
            chain.id = chain_mapping[temp_id]
            chain_mapping[chain.id] = temp_id

    # run Pulchra and FASPR to rebuild the PDB entity
    io = PDB.PDBIO()
    io.set_structure(pdb_model)
    for ent in pdb_entities:
        for index, chain_id in enumerate(ent.get_chain_ids()):
            select_chains = [chain_id]
            chain_type = determine_chain_type_by_pdb_content(pdb_model, chain_id)
            ent.set_chain_type(index, chain_type)
            if chain_type == 'protein' or chain_type == 'peptide': 
                tmp_pdb = f'{pdbid}_prot_{chain_id}.pdb'
                io.save(tmp_pdb, select=SelectChains(select_chains))
                # temporarily reindex from 1
                reindex_pdb(tmp_pdb, start_index=1)
                run_pulchra(pulchra, tmp_pdb)
                reformat_pulchra_rebuilt_pdb(tmp_pdb, chain_id)
                # extract fasta sequence for each protein chain in the structure
                tmp_fasta = f'{pdbid}_prot_{chain_id}.fasta'
                pdb_to_fasta(tmp_pdb, tmp_fasta)
                headers, seq_dict = read_fasta(tmp_fasta)
                pdb_seq = seq_dict[headers[0]]
                pdb_res_pos = align_pdb_seq_to_fas_seq(pdb_seq, ent.get_fas_seq())
                # reindex according to alignment to fas-seq
                reindex_pdb_by_list(tmp_pdb, pdb_res_pos)
                start, end = pdb_res_pos[0], pdb_res_pos[-1]
                filled_pdb_seq = ent.get_fas_seq()[start-1:end]
                ent.set_real_pdb_seq(index, pdb_seq)
                ent.set_filled_pdb_seq(index, filled_pdb_seq)

                # set short filled_pdb_seq's chain type as peptide
                if len(filled_pdb_seq) < max_len_pep:
                    ent.set_chain_type(index, 'peptide')
                
                # create pdb fasta seq
                tmp_fasta = f'{pdbid}_pdb_{chain_id}.fasta'
                with open(tmp_fasta, 'w') as f:
                    f.write(f'>{chain_id}\n{filled_pdb_seq}\n')
                tmp_abrsa_out = f'{pdbid}_{chain_id}.abrsa'
                # run AbRSA to determine if the pdb seq has VH, VL, or their combination
                abrsa_type, ab_res = run_abrsa(abrsa, tmp_fasta, tmp_abrsa_out)

                # re-determine the ab/ag type for chain entity when possible
                if abrsa_type == 'heavy_light':
                    tmp_H_id, tmp_L_id = chain_id + 'H', chain_id + 'L'
                    ent.set_vgene_subgroup(index, ['N.A.', 'N.A.'])
                    ent.set_abrsa_type(index, 'heavy_light')
                    # calculate mean radius of the chain
                    mean_radius = calculate_mean_radius(tmp_pdb)
                    ent.set_mean_radius(index, mean_radius)
                elif abrsa_type == 'heavy':
                    subgroup_H = run_abalign(abalign, tmp_fasta, pdbid, chain_id, chain_type='heavy')
                    ent.set_vgene_subgroup(index, [subgroup_H, 'N.A.'])
                    ent.set_abrsa_type(index, 'heavy')
                    tmalign_out = f'{pdbid}_{chain_id}.tmalign'
                    tmscore = run_tmalign(tmalign, tmp_pdb, struct_ref_vh, tmalign_out)
                    if tmscore < min_tmscore_ab_domain and ent.get_real_pdb_len(index) > min_len_ab_chain:
                        ent.set_abrsa_type(index, 'antigen')
                    if verbose:
                        print(f'tmscore between chain {chain_id} and struct_ref_vh ({struct_ref_vh}): {tmscore}')
                    # calculate mean radius of the chain
                    mean_radius = calculate_mean_radius(tmp_pdb)
                    ent.set_mean_radius(index, mean_radius)
                elif abrsa_type == 'light':
                    subgroup_L = run_abalign(abalign, tmp_fasta, pdbid, chain_id, chain_type='light')
                    ent.set_vgene_subgroup(index, ['N.A.', subgroup_L])
                    ent.set_abrsa_type(index, 'light')
                    tmalign_out = f'{pdbid}_{chain_id}.tmalign'
                    tmscore = run_tmalign(tmalign, tmp_pdb, struct_ref_vl, tmalign_out)
                    if tmscore < min_tmscore_ab_domain and ent.get_real_pdb_len(index) > min_len_ab_chain:
                        ent.set_abrsa_type(index, 'antigen')
                    if verbose:
                        print(f'tmscore between chain {chain_id} and struct_ref_vl ({struct_ref_vl}): {tmscore}')
                    # calculate mean radius of the chain
                    mean_radius = calculate_mean_radius(tmp_pdb)
                    ent.set_mean_radius(index, mean_radius)
                else:
                    if ent.get_abrsa_type(index) == 'N.A.':
                        ent.set_abrsa_type(index, 'antigen')
                    else: # the sequence was previously determined as an ab chain based on fas-seq
                        pdb_ab_res = list(set(pdb_res_pos) & set(ent.get_fas_ab_res()))
                        pdb_ab_res.sort()
                        pdb_ab_res_num, pdb_ab_res_ratio = len(pdb_ab_res), len(pdb_ab_res)/len(pdb_res_pos)
                        if verbose:
                            print(f'chain_id: {chain_id}, pdb_ab_res_num: {pdb_ab_res_num}, pdb_ab_res_ratio: {pdb_ab_res_ratio}')
                        if pdb_ab_res_num < 5 or pdb_ab_res_ratio < 0.2:
                            ent.set_abrsa_type(index, 'antigen')
                        # calculate mean radius of the chain
                        mean_radius = calculate_mean_radius(tmp_pdb)
                        ent.set_mean_radius(index, mean_radius)

                    ent.set_vgene_subgroup(index, ['N.A.', 'N.A.'])

                if verbose:
                    print(f'protein chain_id: {chain_id}, filled_pdb_seq: {ent.get_filled_pdb_seq(index)}, filled_pdb_len: {ent.get_filled_pdb_len(index)}, real_pdb_len: {ent.get_real_pdb_len(index)}, abrsa_type (by filled_pdb_seq): {abrsa_type}, mean_radius: {ent.get_mean_radius(index)}')

            else:
                tmp_pdb = f'{pdbid}_nonprot_{chain_id}.pdb'
                io.save(tmp_pdb, select=SelectChains(select_chains))
                reindex_pdb(tmp_pdb, start_index=1)
                ent.set_real_pdb_seq(index, ent.get_fas_seq())
                ent.set_filled_pdb_seq(index, ent.get_fas_seq())
                abrsa_type = 'antigen'
                ent.set_abrsa_type(index, 'antigen')

    prot_chains = glob.glob(f'{pdbid}_prot_*.pdb')
    if prot_chains:
        # FASPR can only deal with protein chains, ignoring non-protein chains
        os.system(f'cat {pdbid}_prot_*.pdb > {pdbid}_prot.pdb 2>/dev/null')
        run_faspr_repack_pdb(faspr, f'{pdbid}_prot.pdb')
        
        os.system(f'cat {pdbid}_prot.pdb {pdbid}_nonprot_*.pdb > {pdbid}.pdb 2>/dev/null')
    else:
        clean_all_temp_files()
        print(f'done parsing (no protein chain found in {pdbid})')
        exit(0)

    # compute CDR residues using AbRSA_pdb based on structure
    for ent in pdb_entities:
        compute_cdr_residues(ent, pdbid, pdb_model, abrsa_pdb)

    all_paired_HLs = build_all_ab_chain_pairs(pdbid, unidesign, pdb_entities)

    if all_paired_HLs:
        mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title = fetch_pdb_web_info(pdb_fasta_path, pdbid)
        update_HL_ab_type(all_paired_HLs, pdb_entities, pdbid, title)
    
    if verbose:
        print(f'all_paired_HLs: {all_paired_HLs}')
    
    # pair Ab chains with Ag chains
    all_paired_ab_ags = []
    for paired_HLs in all_paired_HLs:
        build_ab_ag_pairs(unidesign, pdbid, paired_HLs, pdb_entities, all_paired_ab_ags)
    if verbose:
        print('all_paired_ab_ags:', all_paired_ab_ags)

    # output paired Ab-Ag interactions or Ab info
    if all_paired_ab_ags:
        if verbose:
            print(', '.join([f'pdbid: {pdbid}', f'mutation(s): {mut_status}', f'classification: {class_}', f'deposit_date: {deposit_date}', f'release_date: {release_date}',
                f'method: {method}', f'resolution: {resolution}', f'r_free: {r_free}', f'r_work: {r_work}', f'pmid: {pmid}', f'doi: {doi}']))

        print(f'ab-ag interactions found and written to file')
        l2code = pdbid[1:3]
        if not os.path.exists(l2code): os.mkdir(l2code)
        header_line = ['PDB_ID', 'Title', 'Mutation(s)', 'Classification', 'Deposit_date', 'Release_date', 'Method', 'Resolution', 'R_free', 'R_work', 'PMID', 'DOI', 
                'Asym_ID_type', 'Ab_type', 'H_subgroup', 'L_subgroup', 'H_chain_ID', 'L_chain_ID', 'H_fas_seq', 'L_fas_seq', 'H_filled_pdb_seq', 'L_filled_pdb_seq', 'H_mean_radius', 'L_mean_radius', 
                'H_fas_seq_len', 'L_fas_seq_len', 'H_pdb_seq_len', 'L_pdb_seq_len', 'H_filled_seq_len', 'L_filled_seq_len', 'HL_inf_res_num', 'H_mol_name', 'L_mol_name', 'H_species', 'L_species', 
                'Ag_chain_ID(s)', 'Ag_type(s)', 'Ag_mol_name(s)', 'Ag_species', 'Ab_ag_inf_res_num', 'CDR_inf_res_num', 'CDR_inf_res_ratio\n']
        f_all = open(f'{l2code}/{pdbid}_aai_all.tsv', 'w')
        f_rep = open(f'{l2code}/{pdbid}_aai_rep.tsv', 'w')
        f_pairs = open(f'{l2code}/{pdbid}_paired_ab_ag_ids.tsv', 'w')
        f_all.write('\t'.join(header_line))
        f_rep.write('\t'.join(header_line))
        for paired_ab_ags in all_paired_ab_ags:
            # write representative ab-ag to file
            rep_ndx = find_representative_ab_ag_pair_index(paired_ab_ags, pdb_entities)
            rep_HL, rep_ags = paired_ab_ags[rep_ndx][0], paired_ab_ags[rep_ndx][1]
            H_chain_id, L_chain_id, HL_inf_res_num, ab_type = rep_HL[0], rep_HL[1], rep_HL[2], rep_HL[4]
            write_pdb_info(f_rep, pdbid, mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title, asym_id_type)
            write_ab_info(f_rep, pdb_entities, chain_mapping, H_chain_id, L_chain_id, HL_inf_res_num, ab_type)
            write_ag_info(f_rep, pdb_entities, chain_mapping, rep_ags)

            # write all paired ab-ag to file
            for i, paired_ab_ag in enumerate(paired_ab_ags):
                cur_HL, cur_ags = paired_ab_ag[0], paired_ab_ag[1]
                H_chain_id, L_chain_id, HL_inf_res_num, ab_type = cur_HL[0], cur_HL[1], cur_HL[2], cur_HL[4]
                write_pdb_info(f_all, pdbid, mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title, asym_id_type)
                write_ab_info(f_all, pdb_entities, chain_mapping, H_chain_id, L_chain_id, HL_inf_res_num, ab_type)
                write_ag_info(f_all, pdb_entities, chain_mapping, cur_ags)
            
            # write paired ab-ag ids to file
            write_paired_ab_ag_ids(f_pairs, paired_ab_ags, chain_mapping, rep_ndx)
        f_rep.close()
        f_all.close()
        f_pairs.close()
    else:
        print('no ab-ag interaction or ab found')
    
    # clean all temporary files
    clean_all_temp_files()

    print('done parsing (normal exit)')
    return


#  MAIN 
if __name__ == '__main__':
    warnings.simplefilter('ignore', BiopythonWarning)
    usage = f'Usage: python {sys.argv[0]} [options] pdbid'
    
    # add arguments and options
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('pdbid', type=str, help='Input a legal PDB ID')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', default=False, 
            help='Do not print status messages')
    parser.add_argument('-v', '--verbose', action='store_true', default=True, 
            help='Print status messages (useful for debugging)')
    parser.add_argument('-a', '--asym_id_type', type=str, choices=['auto', 'auth', 'label'], default='auto', 
            help='Choose an asym_id_type for parsing chains in structure')
    args = parser.parse_args()
    pdbid, verbose, asym_id_type = '', False, 'auto'
    
    # parse arguments and options
    if args.pdbid: 
        pdbid = args.pdbid
        pdbid = pdbid.lower()
    else:
        quit()
    
    verbose = args.verbose
    asym_id_type = args.asym_id_type
    
    # set PDB mmCIF and FASTA paths respectively
    abs_path     = os.path.dirname(os.path.realpath(__file__))
    pdb_fasta_path = f'{abs_path}/../database/fasta_divided'
    pdb_cif_path   = f'{abs_path}/../database/mmCIF_divided'
    
    # copy Abalign/lib because Abalign will always search from lib/ from the current path
    if not os.path.exists('lib'):
        os.system(f'cp -r {abalign_lib} .')
    
    # run a set of programs to parse ab-ag interactions (AAIs)
    parse_ab_ag_interaction(pdbid, pdb_fasta_path, pdb_cif_path)
