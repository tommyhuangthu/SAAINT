#!/usr/bin/env python3
import os, sys

def get_rep_ab_ag_pair_chain_ids(pdbid, aff_H_chain_id, aff_L_chain_id):
    l2code = pdbid[1:3]
    file1 = 'database/saaint_divided/{}/{}_paired_ab_ag_ids.tsv'.format(l2code, pdbid)
    with open(file1, 'r') as f:
        for line in f.readlines():
            string = ',{},{},'.format(aff_H_chain_id, aff_L_chain_id)
            if string in line:
                pair_strs = line.rstrip().split('\t')
                for pair_str in pair_strs:
                    strs = pair_str.split(',')
                    flag = strs[0]
                    if flag == 'rep':
                        return strs[1], strs[2]
    return 'N.A.', 'N.A.'


def merge_affinity_data_and_rep_chain_ids(old_affinity_file, saaint_file, rep_affinity_file):
    lines_affinity = open(old_affinity_file, 'r').readlines()
    lists1 = []
    for i, line1 in enumerate(lines_affinity):
        strs1 = line1.rstrip('\n').split('\t')
        lists1.append(strs1)

    lines_saaint = open(saaint_file, 'r').readlines()
    lists2 = []
    for i, line1 in enumerate(lines_saaint):
        strs1 = line1.rstrip('\n').split('\t')
        lists2.append(strs1)

    fo = open(rep_affinity_file, 'w')
    fo.write('\t'.join(['PDB_ID', 'PMID', 'DOI', 'Asym_ID_type', 'H_chain_ID', 'L_chain_ID', 'Ag_chain_ID(s)', 'Ag_type(s)', 'Affinity_KD(nM)', 'Affinity_method', 'Affinity_temp(K)', 'Affinity_notes\n']))
    for i, list1 in enumerate(lists1):
        pdbid = list1[0]
        if pdbid.startswith('PDB_ID'): continue
        aff_H_chain_id, aff_L_chain_id = list1[1], list1[2]
        rep_H_chain_id, rep_L_chain_id = get_rep_ab_ag_pair_chain_ids(pdbid, aff_H_chain_id, aff_L_chain_id)
        if rep_H_chain_id != 'N.A.' or rep_L_chain_id != 'N.A.':
            found = False
            for j, list2 in enumerate(lists2):
                pdbid2, pmid, doi, asym_id, H_chain_id, L_chain_id, ag_chain_id, ag_type = list2[0], list2[10], list2[11], list2[12], list2[16], list2[17], list2[35], list2[36]
                aff_kd, aff_method, aff_temp, aff_remark = list1[4], list1[5], list1[6], list1[7]
                if pdbid2.startswith('PDB_ID'): continue
                if pdbid == pdbid2 and rep_H_chain_id == H_chain_id and rep_L_chain_id == L_chain_id:
                    found = True
                    break
            if found == False:
                print(f'did not find saaint data for {list1}, further check needed')
            else:
                fo.write('\t'.join([pdbid, pmid, doi, asym_id, rep_H_chain_id, rep_L_chain_id, ag_chain_id, ag_type, aff_kd, aff_method, aff_temp, aff_remark+'\n']))
    fo.close()
    return


if __name__ == '__main__':
    old_aff_file = 'manual_ab_ag_affinity.tsv'
    saaint_file = 'saaintdb_2025011208_all.tsv'
    rep_aff_file=  'saaintdb_affinity_rep.tsv'
    merge_affinity_data_and_rep_chain_ids(old_aff_file, saaint_file, rep_aff_file)
