#!/usr/bin/env python3
import os, sys, datetime

# the absolute path (abs_path) of this python script
abs_path = os.path.dirname(os.path.realpath(__file__))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} assemblies_folder'.format(sys.argv[0]))
        exit(-1)

    cif_folder = sys.argv[1]
    # path to local assemblies folder
    db_path = f'{abs_path}/../database'
    mmcif_path = f'{db_path}/{cif_folder}'
    
    if not os.path.exists(mmcif_path):
        os.mkdir(mmcif_path)
    
    # rsync local mmCIF folder
    os.system(f'rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/assemblies/mmCIF/divided/ {mmcif_path}')

    # add timestamp
    datetime_str = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
    history = f'{db_path}/assemblies_mmCIF_update_time'
    with open(history, 'w') as f:
        f.write(f'{datetime_str}\n')
