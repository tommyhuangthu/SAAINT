# SAAINT
This package provides the source code for SAAINT-parser and SAAINT-DB.

## Installation and running SAAINT-parser
1. Clone this repository
   ```bash
    git clone https://github.com/tommyhuangthu/SAAINT.git
    ```
1. Rsync mmCIF files from the Protein Data Bank (PDB) and download mmCIF-associated FASTA files
   
   ```bash
    ./sbatch_rsync_mmcifs.sh
    ```
   The sbatch_rsync_mmcifs.sh bash script will ```call scripts/run_rsync_mmcifs.py``` to download mmCIF files from PDBe.
   In fact, users can choose to download mmCIFs from RCSB PDB, PDBe, or PDBj (see instructions at
   https://www.wwpdb.org/ftp/pdb-ftp-sites).

   * After rsyncing mmCIF files, users can download the mmCIF-associated FASTA files (if available) using the following command
   ```bash
    rsync -rlpt -v -z --delete ftp.pdbj.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
    ```

1. Run and test SAAINT-parser

1. Run SAAINT-parser to process all mmCIF files

1. Build SAAINT-DB

## Analyzing SAAINT-DB
1. Analyze SAAINT-DB
   ```bash
    python3 scripts/run_saaint_analyzer.py saaintdb_2025012908_all.xlsx -j <task_name>
    ```
