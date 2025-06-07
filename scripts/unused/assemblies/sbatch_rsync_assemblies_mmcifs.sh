#!/bin/bash
#SBATCH --job-name=rsync_assemblies_mmCIF
#SBATCH --error=/home/xiaoqiah/turbo/work/SAAINT/record/rsync_assemblies_mmCIF.err
#SBATCH --output=/home/xiaoqiah/turbo/work/SAAINT/record/rsync_assemblies_mmCIF.out
#SBATCH --account=jiex99
#SBATCH --partition=standard
#SBATCH --nodes=1                 ## how many computers do you need
#SBATCH --ntasks-per-node=3       ## how many processors do you need on each computer
#SBATCH --mem=5G                  ## how many cpus or processors do you need on each computer
#SBATCH --time=72:00:00           ## how long does this need to run (DD-HH:MM:ss)
/home/xiaoqiah/turbo/work/SAAINT/scripts/run_rsync_assemblies_mmcifs.py  assemblies_mmCIF_divided
