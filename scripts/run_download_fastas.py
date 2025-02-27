#!/usr/bin/env python3
import glob, os, sys, datetime
from utils import get_l2codes, get_pdbids_by_l2code, is_fasta_legal


def get_undownloaded_fasta_pdbids(fasta_path, l2code, cif_pdbid_lst):
    if not os.path.exists(os.path.join(fasta_path, l2code)):
        os.mkdir(os.path.join(fasta_path, l2code))
        return cif_pdbid_lst
    lst = []
    for pdbid in cif_pdbid_lst:
        if not os.path.exists(os.path.join(fasta_path, l2code, f'{pdbid}.fasta')):
            lst.append(pdbid)
    return lst

def download_pdb_fastas_in_bulk(undownloaded_pdbid_lst, fasta_path, l2code):
    # create a text file to record all entries with fastas to be downloaded
    file_name = f'fastas_for_download.txt'
    with open(file_name, 'w') as f:
        for entcode in undownloaded_pdbid_lst:
            f.write(f'https://www.rcsb.org/fasta/entry/{entcode}\n')
    # use wget to download fastas with the text file as input
    os.system(f'wget -i {file_name}')
    
    # deal with the downloaded fastas
    # and specifically check if the file is in legal fasta format
    for entcode in undownloaded_pdbid_lst:
        if is_fasta_legal(f'{entcode}'):
            os.system(f'mv {entcode} {fasta_path}/{l2code}/{entcode}.fasta')
        else:
            os.remove(f'{entcode}')
    os.remove(file_name)
    return

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
    for l2code in cif_l2codes:
        cif_pdbids = get_pdbids_by_l2code(cif_path, l2code)
        undownloaded_pdbids = get_undownloaded_fasta_pdbids(fas_path, l2code, cif_pdbids)
        if undownloaded_pdbids:
            download_pdb_fastas_in_bulk(undownloaded_pdbids, fas_path, l2code)

    # add timestamp
    datetime_str = datetime.datetime.now().strftime('%Y-%m-%d|%H:%M:%S')
    history = f'{db_path}/fasta_update_time'
    with open(history, 'w') as f:
        f.write(f'{datetime_str}\n')
