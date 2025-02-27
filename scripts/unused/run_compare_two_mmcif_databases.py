#!/usr/bin/env python3
import glob, os

if __name__ == '__main__':
    abs_path = os.path.dirname(os.path.realpath(__file__))
    old_dir = f'{abs_path}/../database/mmCIF_divided'
    new_dir = f'{abs_path}/../database/mmCIF_divided_beta'

    new_entries = []
    obs_entries = []
    
    # find new entries in current mmCIF repository relative to a previous repo.
    for l2full in glob.glob(f'{new_dir}/*'):
        l2code = l2full[-2::1]
        for ent in glob.glob(f'{l2full}/*'):
            entcode = ent[-11:-7]
            # exists in new but not old -> new entry
            if not os.path.exists(f'{old_dir}/{l2code}/{entcode}.cif.gz'):
                new_entries.append(entcode)

    # find obsolete entries in current mmCIF repo. relative to a previous repo.
    for l2full in glob.glob(f'{old_dir}/*'):
        l2code = l2full[-2::1]
        for ent in glob.glob(f'{l2full}/*'):
            entcode = ent[-11:-7]
            # exists in old but not new -> obsolete entry
            if not os.path.exists(f'{new_dir}/{l2code}/{entcode}.cif.gz'):
                obs_entries.append(entcode)

    #print('new PDB entries are:\n', new_entries)
    #print('obsolete PDB entries are:\n', obs_entries)

    print('number of new mmcif entries:', len(new_entries))
    with open(f'{abs_path}/../database/list_new_cifs.txt', 'w') as f:
        for ent in new_entries:
            f.write(f'{ent}\n')

    print('number of obsolete mmcif entries:', len(obs_entries))
    with open(f'{abs_path}/../database/list_obsolete_cifs.txt', 'w') as f:
        for ent in obs_entries:
            f.write(f'{ent}\n')
