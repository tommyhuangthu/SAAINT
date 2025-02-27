#!/usr/bin/env python3
import glob, os, sys
from utils import get_l2codes, get_pdbids_by_l2code, is_fasta_legal


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: {} mmCIF_folder fasta_folder'.format(sys.argv[0]))
        exit(-1)

    cif_folder = sys.argv[1]
    fas_folder = sys.argv[2]

    db_path = f'/home/xiaoqiah/turbo/work/SAAINT/database'
    os.chdir(db_path)

    cif_path = f'{db_path}/{cif_folder}'
    fas_path = f'{db_path}/{fas_folder}'
    
    cif_l2codes = get_l2codes(cif_path)
    print(cif_l2codes)
    for l2code in cif_l2codes:
        cif_pdbids = get_pdbids_by_l2code(cif_path, l2code)
        for pdbid in cif_pdbids:
            fasta = f'{fas_path}/{l2code}/{pdbid}.fasta'
            if not os.path.exists(fasta):
                print(f'{fasta} does not exist, skip')
            else:
                legal = is_fasta_legal(fasta)
                if not legal:
                    print(f'{fasta} is not a legal fasta file')

