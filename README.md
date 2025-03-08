# SAAINT
This package provides the source code for SAAINT-parser and SAAINT-DB.

## Installation and running SAAINT-parser
1. Clone this repository
   ```bash
    git clone https://github.com/tommyhuangthu/SAAINT.git
    ```
1. Rsync mmCIF files from the Protein Data Bank (PDB)


1. Download mmCIF-associated FASTA files

1. Run/Test SAAINT-parser on a single PDB entry

1. Run SAAINT-parser to process all mmCIF files

1. Build SAAINT-DB

## Analyzing SAAINT-DB
1. Build SAAINT-DB
   ```bash
    python scripts/run_saaint_analyzer.py saaintdb_2025012908_all.xlsx -j <task_name>
    ```
