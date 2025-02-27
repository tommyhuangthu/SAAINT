#!/bin/bash
#SBATCH --job-name=rsync_gzip
#SBATCH --error=rsync_gzip.log
#SBATCH --output=rsync_gzip.out
#SBATCH --account=jiex0
#SBATCH --partition=standard
#SBATCH --nodes=1                 ## how many computers do you need
#SBATCH --ntasks-per-node=8       ## how many processors do you need on each computer
#SBATCH --mem=5G                  ## how many cpus or processors do you need on each computer
#SBATCH --time=72:00:00           ## how long does this need to run (DD-HH:MM:ss)
mkdir -p /scratch/jiex_root/jiex0/xiaoqiah/rsync_gzip/
mv /home/xiaoqiah/turbo/work/Databases/mmCIF_divided-20230519 /scratch/jiex_root/jiex0/xiaoqiah/rsync_gzip/
pigz -p 8 -1 -r /scratch/jiex_root/jiex0/xiaoqiah/rsync_gzip/mmCIF_divided-20230519
mv /scratch/jiex_root/jiex0/xiaoqiah/rsync_gzip/mmCIF_divided-20230519 /home/xiaoqiah/turbo/work/Databases
