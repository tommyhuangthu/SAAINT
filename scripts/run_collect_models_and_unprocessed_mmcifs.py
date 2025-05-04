#!/usr/bin/env python3

import os
from utils import read_list

# get the absolute path of current python script file
abs_path = os.path.dirname(os.path.realpath(__file__))


if __name__ == '__main__':
    # move pdb models to processed models
    if not os.path.exists(f'{abs_path}/../database/processed_models'):
        os.system(f'mkdir {abs_path}/../database/processed_models')
    os.system(f'mv mkdir {abs_path}/../database/saaint_divided/*.pdb {abs_path}/../database/processed_models')
    os.system(f'tar -cvf - {abs_path}/../database/processed_models | split -b 90M - {abs_path}/../database/processed_models.tar.part_')
    
    # copy mmcif.gz to unprocessed_mmcifs
    if not os.path.exists(f'{abs_path}/../database/unprocessed_mmcifs'):
        os.system(f'mkdir {abs_path}/../database/unprocessed_mmcifs')
    ents = read_list('list_saaintdb_pdbs.txt')
    for ent in ents:
        l2code = ent[1:3]
        os.system(f'cp {abs_path}/../database/mmCIF_divided/{l2code}/{ent}.cif.gz {abs_path}/../database/unprocessed_mmcifs')
    
