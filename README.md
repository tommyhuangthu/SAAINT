# SAAINT
This package provides the source code for SAAINT-parser and SAAINT-DB.

## Installation and running SAAINT-parser
1. Clone this repository
   ```bash
    git clone https://github.com/tommyhuangthu/SAAINT.git
    ```
1. Rsync mmCIF files from the Protein Data Bank (PDB)
   Users can choose to rsync mmCIF files from one of the three sites: RCSB PDB, PDBe, or PDBj, by
   running one of the following commands 
   ```bash
    rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
    ```
   ```bash
    rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ ./mmCIF
    ```
   ```bash
    rsync -rlpt -v -z --delete ftp.pdbj.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
    ```

1. Download mmCIF-associated FASTA files

1. Run/Test SAAINT-parser on a single PDB entry

1. Run SAAINT-parser to process all mmCIF files

1. Build SAAINT-DB

## Analyzing SAAINT-DB
1. Analyze SAAINT-DB
   ```bash
    python scripts/run_saaint_analyzer.py saaintdb_2025012908_all.xlsx -j <task_name>
    ```
