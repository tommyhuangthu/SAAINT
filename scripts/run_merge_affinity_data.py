#!/usr/bin/env python3
import os, sys

def get_rep_ab_ag_pair_chain_ids(pdbid):
    l2code = pdbid[1:3]
    file1 = 'database/saaint_divided/{}/{}_paired_ab_ag_ids.tsv'.format(l2code, pdbid)
    model_HL_dict = dict()
    with open(file1, 'r') as f:
        for line in f.readlines():
            HL_pairs = []
            cols = line.rstrip().split('\t')
            for col in cols[1:]:
                strs = col.split(',')
                HL_pairs.append([strs[1], strs[2]])
            key = int(cols[0])
            if key not in model_HL_dict:
                model_HL_dict[key] = [HL_pairs]
            else:
                model_HL_dict[key].append(HL_pairs)
    if pdbid == '9axf':
        print(model_HL_dict)
    return model_HL_dict


def merge_affinity_data_and_rep_chain_ids(old_affinity_file, saaint_file, rep_affinity_file):
    # read old_affinity_file contents into lists1
    lines_affinity = open(old_affinity_file, 'r').readlines()
    lists1 = []
    for i, line1 in enumerate(lines_affinity):
        strs1 = line1.rstrip('\n').split('\t')
        lists1.append(strs1)
    
    # read saaintdb_all.tsv file into lists2
    lines_saaint = open(saaint_file, 'r').readlines()
    lists2 = []
    for i, line1 in enumerate(lines_saaint):
        strs1 = line1.rstrip('\n').split('\t')
        lists2.append(strs1)

    #
    fo = open(rep_affinity_file, 'w')
    fo.write('\t'.join(['PDB_ID', 'PMID', 'DOI', 'Model_index', 'Asym_ID_type', 'H_chain_ID', 'L_chain_ID', 'Ag_chain_ID(s)', 'Ag_type(s)', 'Affinity_KD(nM)', 'Affinity_method', 'Affinity_temp(K)', 'Affinity_notes\n']))
    for i, list1 in enumerate(lists1):
        pdbid = list1[0]
        if pdbid.startswith('PDB_ID'): 
            continue
        aff_H_chain_id, aff_L_chain_id = list1[1], list1[2]
        model_HL_dict = get_rep_ab_ag_pair_chain_ids(pdbid)
        for key in model_HL_dict:
            model_HLs = model_HL_dict[key]
            for model_HL in model_HLs:
                for j, list2 in enumerate(lists2):
                    pdbid2, pmid, doi, model_index, asym_id, H_chain_id, L_chain_id, ag_chain_id, ag_type = list2[0], list2[10], list2[11], list2[12], list2[13], list2[17], list2[18], list2[36], list2[37]
                    aff_kd, aff_method, aff_temp, aff_remark = list1[4], list1[5], list1[6], list1[7]
                    if pdbid2.startswith('PDB_ID'): 
                        continue
                    if pdbid == pdbid2 and [aff_H_chain_id, aff_L_chain_id] in model_HL and [H_chain_id, L_chain_id] in model_HL and ag_chain_id != 'N.A.':
                        fo.write('\t'.join([pdbid, pmid, doi, model_index, asym_id, H_chain_id, L_chain_id, ag_chain_id, ag_type, aff_kd, aff_method, aff_temp, aff_remark+'\n']))
    fo.close()
    return


if __name__ == '__main__':
    old_aff_file = 'saaintdb/manual_ab_ag_affinity.tsv'
    saaint_file = 'saaintdb/saaintdb_20250523_all.tsv'
    new_aff_file=  'saaintdb/saaintdb_affinity_all.tsv'
    merge_affinity_data_and_rep_chain_ids(old_aff_file, saaint_file, new_aff_file)
