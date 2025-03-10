# SAAINT
SAAINT stands for <u>S</u>tructural <u>A</u>ntibody and <u>A</u>ntibody-antigen <u>INT</u>eraction. This package provides the source code for implementing SAAINT-parser, and building, analyzing, and updating the resulting SAAINT database, SAAINT-DB.

## Installation and running SAAINT-parser
1. Configure the environment
   ```bash
   conda create -n saaint python=3.10 biopython=1.84 numpy=1.26.3 matplotlib=3.10.0 seaborn=0.13.2 urllib3=2.2.1 pandas=2.2.3
   conda activate saaint
   ```
   
1. Clone this repository
   ```bash
   git clone https://github.com/tommyhuangthu/SAAINT.git
   cd SAAINT/
   mkdir -p database/mmCIF_divided database/fasta_divided database/saaint_divided
   ```
   
1. Rsync mmCIF files from the Protein Data Bank (PDB) and download mmCIF-associated FASTA files
   * Download mmCIF files from the PDB
   ```bash
   ./scripts/sbatch_rsync_mmcifs.sh
   ```
   The SLURM script ```scripts/sbatch_rsync_mmcifs.sh``` will call ```scripts/run_rsync_mmcifs.py``` to download mmCIF files from the PDBe.
   The downloaded mmCIF files are saved into the database/mmCIF_divided directory, and specifically, a set of folders named with the two middle letters of PDB entries.
   Users can also modify ```scripts/run_rsync_mmcifs.py``` to download mmCIFs from RCSB PDB, PDBe, or PDBj (see the instructions at
   https://www.wwpdb.org/ftp/pdb-ftp-sites), depending on their locations.

   * Download the associated FASTA files
   ```bash
   ./scripts/sbatch_update_fastas.sh
   ```
   The SLURM script ```scripts/sbatch_update_fastas.sh``` will call ```scripts/run_update_fastas.py``` to download or update FASTA files
   from the RCSB PDB website. For example, for PDB entry ```5zxv```, the script can automatically retrieve its FASTA content from ```[https://www.rcsb.org/fasta/entry/5zxv](https://www.rcsb.org/fasta/entry/5ZXV/display)```
   and save it into ```database/fasta_divided/zx/5zxv.fasta```.

1. Run and test SAAINT-parser
   For example, to extract the antibodies and antibo
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
