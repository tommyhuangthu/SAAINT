#!/usr/bin/env python3
import glob, os, datetime
from utils import is_fasta_legal, parse_rsync_mmcif_out

def update_pdb_fastas(fasta_path, update_dict):
    for l2code in update_dict:
        # create a text file to record all entries with fastas to be downloaded
        file_name = f'fastas_for_download.txt'
        with open(file_name, 'w') as f:
            for entcode in update_dict[l2code]:
                f.write(f'https://www.rcsb.org/fasta/entry/{entcode}\n')
        os.system(f'wget -i {file_name}')
        
        # check if the file is in legal fasta format
        for entcode in update_dict[l2code]:
            if is_fasta_legal(f'{entcode}'):
                os.system(f'mv {entcode} {fasta_path}/{l2code}/{entcode}.fasta')
            else:
                os.remove(f'{entcode}')
        os.remove(file_name)
    return


if __name__ == '__main__':
    abs_path = os.path.dirname(os.path.realpath(__file__))
    obsolete_list, update_dict = parse_rsync_mmcif_out(f'{abs_path}/../record/rsync_mmCIF.out')

    db_path = f'{abs_path}/../database'

    with open(f'{db_path}/list_obsolete_cifs.txt', 'w') as f:
        for ent in obsolete_list:
            f.write(f'{ent}\n')
    print(f'number of obsolete mmcif entries: {len(obsolete_list)}')
    
    num_update_entries = 0
    with open(f'{db_path}/list_update_cifs.txt', 'w') as f:
        for l2code in update_dict:
            num_update_entries += len(update_dict[l2code])
            for ent in update_dict[l2code]:
                f.write(f'{ent}\n')
    print(f'number of update mmcif entries: {num_update_entries}')
    
    # download fasta files for updated mmcif entries
    fasta_path = f'{abs_path}/../database/fasta_divided'
    update_pdb_fastas(fasta_path, update_dict)

    # add timestamp
    datetime_str = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
    history = f'{db_path}/fasta_update_time'
    with open(history, 'w') as f:
        f.write(f'{datetime_str}\n')
