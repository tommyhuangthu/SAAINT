# SAAINT
This package provides the source code for SAAINT-parser and SAAINT-DB.

## Installation and running SAAINT-parser
1. Clone this repository
   ```bash
    git clone https://github.com/tommyhuangthu/SAAINT.git
    ```
1. Rsync mmCIF files from the Protein Data Bank (PDB) and download mmCIF-associated FASTA files
   * rsync PDB mmCIF files
   ```bash
    ./sbatch_rsync_mmcifs.sh
    ```
   The provided bash script ```sbatch_rsync_mmcifs.sh``` will call ```scripts/run_rsync_mmcifs.py``` to download mmCIF files from PDBe.
   However, users can modify ```scripts/run_rsync_mmcifs.py``` to download mmCIFs from RCSB PDB, PDBe, or PDBj (follow the instructions at
   https://www.wwpdb.org/ftp/pdb-ftp-sites).

   * Download associated FASTA files
   ```bash
    ./sbatch_update_fastas.sh
    ```
   The provided bash script ```sbatch_update_fastas.sh``` will call ```scripts/run_update_fastas.py``` to download or update FASTA files
   from the RCSB PDB website. For example, the script can automatically create ```database/fastas/5zxv.fasta``` by retrieving the FASTA content
   from ```https://www.rcsb.org/fasta/entry/5zxv```, where 5zxv is a valid PDB entry.

1. Run and test SAAINT-parser
   ```bash
    python3 ./scripts/run_saaint_parser.py 5zxv
    ```
   or:
   ```bash
    python3 ./scripts/run_saaint_parser.py -v 5zxv
    ```
   This command will print out the calculation details, which can help debugging.

1. Run SAAINT-parser to process all mmCIF files
   ```bash
    python3 ./scripts/run_submit_saaint_parser_jobs.py -path <mmCIF_path> <work_dir> <n_cpus>
    ```
   or:
   ```bash
    python3 ./scripts/run_submit_saaint_parser_jobs.py -list <mmCIF_list> <work_dir> <n_cpus>
    ```

1. Build SAAINT-DB
   ```bash
    python3 ./scripts/run_saaintdb_builder.py
    ```

## Analyzing SAAINT-DB
1. Analyze SAAINT-DB
   ```bash
    python3 scripts/run_saaint_analyzer.py saaintdb_2025012908_all.xlsx -j <task_name>
    ```
