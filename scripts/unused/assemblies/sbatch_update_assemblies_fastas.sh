#!/bin/bash
#SBATCH --job-name=/home/xiaoqiah/turbo/work/SAAINT/record/update_assemblies_fastas
#SBATCH --error=/home/xiaoqiah/turbo/work/SAAINT/record/update_assemblies_fastas.err
#SBATCH --output=/home/xiaoqiah/turbo/work/SAAINT/record/update_assemblies_fastas.out
#SBATCH --account=jiex99
#SBATCH --partition=standard
#SBATCH --nodes=1                 ## how many computers do you need
#SBATCH --ntasks-per-node=1       ## how many processors do you need on each computer
#SBATCH --mem=5G                  ## how many cpus or processors do you need on each computer
#SBATCH --time=48:00:00           ## how long does this need to run (DD-HH:MM:ss)
/home/xiaoqiah/turbo/work/SAAINT/scripts/run_update_assemblies_fastas.py
