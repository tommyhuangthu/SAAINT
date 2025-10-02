#!/usr/bin/env python3
import os, pandas, subprocess

abs_path = os.path.dirname(os.path.realpath(__file__))
db_path = f'{abs_path}/../database'


def get_mmcif_update_date(update_time_file):
    with open(update_time_file, 'r') as f:
        update_date = f.readline().rstrip().split()[0]
        return update_date


if __name__ == '__main__':
    dts = get_mmcif_update_date(f'{db_path}/mmCIF_update_time')
    for tag in ['all']:
        tmp = f'saaintdb/saaintdb_{dts}_{tag}.tmp'
        xlsx = tmp.replace('.tmp', '.xlsx')
        tsv = tmp.replace('.tmp', '.tsv')
        
        folder_name='saaint_divided'
        subprocess.run([f'cat {db_path}/{folder_name}/[0-9a-z][0-9a-z]/*_{tag}.tsv | head -n 1 > {tmp}'], 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subprocess.run([f'cat {db_path}/{folder_name}/[0-9a-z][0-9a-z]/*_{tag}.tsv | grep -v PDB_ID >> {tmp}'], 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        df = pandas.read_csv(tmp, sep='\t', keep_default_na=False)
        df = df.sort_values(by=['Deposit_date', 'Release_date', 'PDB_ID', 'Ag_chain_ID(s)', 'H_chain_ID', 'L_chain_ID'], ascending=[True, True, True, True, True, True])
        print(df)
        
        df.to_excel(f'{xlsx}', sheet_name=f'{tag}', na_rep='', index=False)
        df.to_csv(f'{tsv}', sep='\t', na_rep='', index=False)
        os.system(f'rm {tmp}')
