#!/usr/bin/env python3

import os
from utils import read_list

# get the absolute path of current python script file
abs_path = os.path.dirname(os.path.realpath(__file__))

def get_mmcif_update_date(update_time_file):
    with open(update_time_file, 'r') as f:
        update_date = f.readline().rstrip().split()[0]
        return update_date


if __name__ == '__main__':
    dts = get_mmcif_update_date(f'{abs_path}/../database/mmCIF_update_time')
    
    # move rsync models to processed models
    if not os.path.exists(f'{abs_path}/../database/models/processed_pdb_models_{dts}'):
        os.system(f'mkdir -p {abs_path}/../database/models/processed_pdb_models_{dts}')
    os.system(f'rsync -az {abs_path}/../database/saaint_divided/*.pdb {abs_path}/../database/models/processed_pdb_models_{dts}/')
    os.system(f'tar -czf {abs_path}/../database/models/processed_pdb_models_{dts}.tar.gz {abs_path}/../database/models/processed_pdb_models_{dts}')
    
    # copy mmcif.gz to unprocessed_mmcifs
    if not os.path.exists(f'{abs_path}/../database/models/unprocessed_mmcifs_{dts}'):
        os.system(f'mkdir {abs_path}/../database/models/unprocessed_mmcifs_{dts}')
    # create saaintdb pdb list
    os.system(f'cut -f1 {abs_path}/../saaintdb/saaintdb_{dts}_all.tsv | grep -v PDB_ID | sort -u > {abs_path}/../list_saaintdb_pdbs_{dts}.txt')
    ents = read_list(f'{abs_path}/../list_saaintdb_pdbs_{dts}.txt')
    for ent in ents:
        l2code = ent[1:3]
        os.system(f'cp {abs_path}/../database/mmCIF_divided/{l2code}/{ent}.cif.gz {abs_path}/../database/models/unprocessed_mmcifs_{dts}')
    os.system(f'gunzip -r {abs_path}/../database/models/unprocessed_mmcifs_{dts}')
    os.system(f'tar -czf {abs_path}/../database/models/unprocessed_mmcifs_{dts}.tar.gz {abs_path}/../database/models/unprocessed_mmcifs_{dts}')

    
