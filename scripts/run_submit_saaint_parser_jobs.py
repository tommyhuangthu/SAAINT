#!/usr/bin/env python3
import os, sys, subprocess, glob
from utils import read_list, get_l2codes, get_pdbids_by_l2code, abalign_lib, parse_rsync_mmcif_out

# get the absolute path of current python script file
abs_path = os.path.dirname(os.path.realpath(__file__))

def delete_saaint_results(saaint_path, obsolete_list, update_dict):
    for pdbid in obsolete_list:
        if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv'):
            print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv')
            os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv')
        if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv'):
            print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv')
            os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv')
        if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv'):
            print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv')
            os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv')
        for file in glob.glob(f'{saaint_path}/{pdbid}_model_*.pdb'):
            print(f'deleting {file}')
            os.remove(f'{file}')
    for l2code in update_dict:
        for pdbid in update_dict[l2code]:
            if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv'):
                print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv')
                os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_all.tsv')
            if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv'):
                print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv')
                os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_aai_rep.tsv')
            if os.path.exists(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv'):
                print(f'deleting {saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv')
                os.remove(f'{saaint_path}/{pdbid[1:3]}/{pdbid}_paired_ab_ag_ids.tsv')
            for file in glob.glob(f'{saaint_path}/{pdbid}_model_*.pdb'):
                print(f'deleting {file}')
                os.remove(f'{file}')
    return


def allocate_entries_by_cif_path(cif_path, n_cpu=300):
    list_of_entries = []
    for i in range(n_cpu):
        list_of_entries.append([])
    l2codes = get_l2codes(cif_path)
    
    index = 0
    for l2code in l2codes:
        pdbids = get_pdbids_by_l2code(cif_path, l2code)
        for pdbid in pdbids:
            list_of_entries[index].append(pdbid)
            index = index+1
            if index >= n_cpu: index = index - n_cpu

    return l2codes, list_of_entries


def allocate_entries_by_cif_list(list_file, n_cpu=2):
    entries = read_list(list_file)
    list_of_entries = []
    for i in range(n_cpu):
        list_of_entries.append([])

    index = 0
    l2codes = []
    for entry in entries:
        #if os.path.exists(f'/home/xiaoqiah/turbo/work/SAAINT/database/saaint_divided/{entry[1:3]}/{entry}_aai_all.tsv'):
        #    continue
        if not entry[1:3] in l2codes:
            l2codes.append(entry[1:3])
        list_of_entries[index].append(entry)
        index = index + 1
        if index >= n_cpu: index = index - n_cpu

    return l2codes, list_of_entries


def submit_sbatch_jobs(l2codes, list_of_entries, work_dir='./', n_cpu=300):
    if not os.path.exists(f'{work_dir}'):
        os.mkdir(f'{work_dir}')
    os.chdir(f'{work_dir}')

    # copy Abalign lib/ folder to working directory
    if not os.path.exists('lib/'):
        os.system(f'cp -r {abalign_lib} .')

    # create l2code folders to save 
    for l2code in l2codes:
        if not os.path.exists(f'{l2code}'):
            os.mkdir(f'{l2code}')

    # create record/ directory to save sbatch script output and logs
    if not os.path.exists('record/'):
        os.mkdir('record/')
    
    for index in range(len(list_of_entries)):
        if list_of_entries[index]:
            # create sbatch script
            job = f'record/run_saaint.{index}'
            with open(job, 'w') as f:
                f.write('''#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --error=%s.err
#SBATCH --output=%s.out
#SBATCH --account=%s
#SBATCH --partition=standard
#SBATCH --nodes=1                 ## how many computers do you need
#SBATCH --ntasks-per-node=1       ## how many processors do you need on each computer
#SBATCH --mem=5G                 ## how much memory do you need on each computer
#SBATCH --time=24:00:00           ## how long does this need to run (DD-HH:MM:ss)
'''%(job, job, job, 'jiex99'))
                for ent in list_of_entries[index]:
                    f.write(f'{abs_path}/run_saaint_parser.py {ent}\n')

            # submit sbatch script
            subprocess.Popen(['sbatch {}'.format(job)], shell=True)
    return


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: {} -path <mmCIF_path> <work_dir> <n_cpus>\n or\nUsage: {} -list <mmCIF_list> <work_dir> <n_cpus>'.format(sys.argv[0], sys.argv[0]))
        exit(-1)
    
    type_ = sys.argv[1]
    
    if type_ == '-path':
        cif_path = sys.argv[2]
        saaint_path = sys.argv[3]
        n_cpu = int(sys.argv[4])
        l2codes, list_of_entries = allocate_entries_by_cif_path(cif_path, n_cpu=n_cpu)
        submit_sbatch_jobs(l2codes, list_of_entries, work_dir=saaint_path, n_cpu=n_cpu)
    elif type_ == '-list':
        cif_list = sys.argv[2]
        saaint_path = sys.argv[3]
        n_cpu = int(sys.argv[4])
        

        # delete previously calculated saaint results for obsolete entries
        obsolete_list, update_dict = parse_rsync_mmcif_out(f'{abs_path}/../record/rsync_mmCIFs.out')
        #saaint_path = f'{abs_path}/../database/saaint_divided'
        delete_saaint_results(saaint_path, obsolete_list, update_dict)

        l2codes, list_of_entries = allocate_entries_by_cif_list(cif_list, n_cpu=n_cpu)
        submit_sbatch_jobs(l2codes, list_of_entries, work_dir=saaint_path, n_cpu=n_cpu)

