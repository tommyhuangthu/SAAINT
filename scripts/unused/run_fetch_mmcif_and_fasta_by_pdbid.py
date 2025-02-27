#!/usr/bin/env python3
import os, sys

cif_path   = 'fetched_cifs'
fasta_path = 'fetched_fastas'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} pdbid'.format(sys.argv[0]))
        print('Example: {} 5zxv'.format(sys.argv[0]))
        exit(-1)

    pdbid = sys.argv[1]

    os.system(f'wget https://files.rcsb.org/download/{pdbid}.cif ; mv {pdbid}.cif {cif_path}/{pdbid}.cif')
    os.system(f'wget https://www.rcsb.org/fasta/entry/{pdbid} ; mv {pdbid} {fasta_path}/{pdbid}.fasta')

