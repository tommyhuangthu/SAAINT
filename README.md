# SAAINT

SAAINT stands for <strong><ins>S</ins>tructures of <ins>A</ins>ntibodies and Antibody-<ins>A</ins>ntigen <ins>INT</ins>eractions</strong>.

This package provides: 
* Implementation of the ```SAAINT-parser``` workflow, designed for fast and accurate extraction and annotation of antibody (Ab) structures and antibody-antigen interactions (AAIs) from the Protein Data Bank (PDB). This results in a comprehensive and up-to-date structural antibody database, ```SAAINT-DB```.
* Source code for building, analyzing, and updating ```SAAINT-DB```.

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
   mkdir -p record database/mmCIF_divided database/fasta_divided database/saaint_divided
   ```

2. Clone the UniDesign, FASPR, and Pulchra repositories (used by SAAINT-parser)
   
   ```bash
   git clone https://github.com/tommyhuangthu/UniDesign.git
   git clone https://github.com/tommyhuangthu/FASPR.git
   git clone https://github.com/euplotes/pulchra.git
   ```
   Then ```cd``` into each repository and build executibles of UniDesign, FASPR, and pulchra, following their instructions. \
   For clarity, we assume these repositories are saved to the user's home directory (any directory is fine as long as the correct directories are specified for these executibles; see below).
   
1. Modify related directories and paths within the source code before running
   
   Some directories or paths are hard coded within the source code, and thus a direct run of the code may cause errors.
   Assume the user home directory is ```/home/username``` and the SAAINT package is saved as ```/home/username/SAAINT```, the following paths or directories need changes:\
   In ```scripts/sbatch_rsync_mmcifs.sh```: Changes ```/home/xiaoqiah/turbo/work/SAAINT``` to ```/home/username/SAAINT```. \
   In ```scripts/sbatch_update_fastas.sh```: Changes ```/home/xiaoqiah/turbo/work/SAAINT``` to ```/home/username/SAAINT```. \
   In ```scripts/utils.py```: Changes the paths for UniDesign, FASPR, and Pulchra.
   
   
1. Rsync mmCIF files from the Protein Data Bank (PDB) and download mmCIF-associated FASTA files
   
   * Download mmCIF files from the PDB
   ```bash
   sbatch scripts/sbatch_rsync_mmcifs.sh
   ```
   The SLURM script ```scripts/sbatch_rsync_mmcifs.sh``` executes ```scripts/run_rsync_mmcifs.py``` to download mmCIF files from PDBe.
   The downloaded files are stored in the ```database/mmCIF_divided``` directory, organized into subfolders named after the two middle letters of the corresponding PDB entries.
   For example, for PDB entry ```5zxv```, the script automatically retrieves its mmCIF file and saves it as ```database/mmCIF_divided/zx/5zxv.cif.gz```.
   Users can modify ```scripts/run_rsync_mmcifs.py``` to download mmCIF files more efficiently from RCSB PDB, PDBe, or PDBj, depending on their location. For more details, refer to the instructions at [wwPDB](https://www.wwpdb.org/ftp/pdb-ftp-sites).

   * Download the associated FASTA files
   ```bash
   sbatch scripts/sbatch_update_fastas.sh
   ```
   The SLURM script ```scripts/sbatch_update_fastas.sh``` executes ```scripts/run_update_fastas.py``` to download or update FASTA files from the RCSB PDB website (See below for updating).
   The script ```scripts/run_update_fastas.py``` analyzes the output of the SLURM job ```scripts/sbatch_rsync_mmcifs.sh``` to identify the PDB entries that need updating or have become obsolete, generating two mmCIF list files: ```database/list_update_cifs.txt``` and ```database/list_obsolete_cifs.txt```.
   The downloaded FASTA files are stored in the ```database/fasta_divided``` directory, organized into subfolders named after the two middle letters of the corresponding PDB entries.
   For example, for PDB entry ```5zxv```, the script automatically retrieves its FASTA content from [RCSB PDB: 5zxv](https://www.rcsb.org/fasta/entry/5zxv/display) and saves it as ```database/fasta_divided/zx/5zxv.fasta```.

1. Run and test SAAINT-parser

   Users can run SAAINT-parser with a PDB entry as the only input. For example:
   ```bash
   python3 ./scripts/run_saaint_parser.py 5zxv
   ```
   To enable verbose output for debugging, use the ```--verbose``` (or ```-v```) option:
   ```bash
   python3 ./scripts/run_saaint_parser.py -v 5zxv
   ```
   Successfully identified antibodies and antibody-antigen interactions are saved in a folder named after the two middle letters of the PDB entry. For example, for ```5zxv```, the results are stored in the ```zx``` folder:\
   ```zx/5zxv_aai_all.tsv```: contains all identified antibodies or AAIs.\
   ```zx/5zxv_aai_rep.tsv```: records representative antibodies or AAIs (see Reference for details).\
   ```zx/5zxv_paired_ab_ag_ids.tsv```: lists paired antibody-antigen chain IDs, and labels representative and nonrepresentative pairs.

## Building SAAINT-DB

1. Run SAAINT-parser to process all mmCIF files
   
   To process all mmCIF files using SAAINT-parser, run:
   ```bash
   python3 ./scripts/run_submit_saaint_parser_jobs.py -path <mmCIF_path> <work_dir> <n_cpus>
   ```
   ```<mmCIF_path>```: Path to the directory containing the downloaded mmCIF files (e.g., database/mmCIF_divided). \
   ```<work_dir>```: Directory for saving intermediate outputs and final results generated by SAAINT-parser, e.g., ```database/saaint_divided```. \
   ```<n_cpus>```: Number of CPU cores allocated for processing, e.g., ```300``` (A large number of CPUs is recommended for initial construction). \
   We recommend using this command for the initial construction of the SAAINT-DB database from scratch.

   Alternatively, to process a specific set of mmCIF files, run:
   ```bash
   python3 ./scripts/run_submit_saaint_parser_jobs.py -list <mmCIF_list> <work_dir> <n_cpus>
   ```
   ```<mmCIF_list>```: List of mmCIF files to be processed, e.g., ```database/list_update_cifs.txt```. \
   ```<work_dir>``` and ```<n_cpus>```: Same as described above.\
   We recommend using the command with ```-list``` option for updating SAAINT-DB efficiently.
      
1. Collect SAAINT-parser results to build the SAAINT-DB datasets

   To build the SAAINT-DB databases, run:
   ```bash
   python3 ./scripts/run_saaintdb_builder.py
   ```
  This command merges all ```*_aai_all.tsv``` and ```*_aai_rep.tsv``` files to generate the SAAINT-DB full and representative datasets. \
  The SAAINT-DB full dataset is saved as ```saaintdb_{timestamp}_all.xlsx``` and ```saaintdb_{timestamp}_all.tsv```. \
  The representative dataset is saved as ```saaintdb_{timestamp}_rep.xlsx``` and ```saaintdb_{timestamp}_rep.tsv```. Here, ```{timestamp}``` (e.g., ```2025012908```) represents the date and time of the database build. \
  For more details, refer to ```saaintdb/saaintdb_2025012908_all.xlsx```, ```saaintdb/saaintdb_2025012908_all.tsv```, ```saaintdb/saaintdb_2025012908_rep.xlsx```, and ```saaintdb/saaintdb_2025012908_all.tsv```.

## Analyzing SAAINT-DB and antibody-antigen binding affinity data

1. Analyze SAAINT-DB
   
   To make statistical analyses of the SAAINT-DB data entries, run:
   ```bash
   python3 scripts/run_saaintdb_analyzer.py saaintdb/saaintdb_2025012908_all.xlsx -j <jobname>
   ```
   Here, ```<jobname>``` takes the following options:\
   ```date```: Analyzes the deposition and release dates of PDB entries.\
   ```classification```: Examines the classification of PDB entries.\
   ```method```: Analyzes the experimental methods used to determine PDB structures.\
   ```resolution```: analyzes the resolution of X-ray and EM structures\
   ```publication```: Analyzes PDB-associated publication details, including PMID and DOI.\
   ```asym_id```: Analyzes the types of asym_id used to parse PDB entries.\
   ```plot_pdb_num```: Plots the number of PDB entries over released years.\
   ```count_pdb_num```: Counts the number of PDB entries per release year.\
   ```count_data_num```: Counts the total number of data entries.\
   ```ab_spe```: Analyzes the top species for antibody heavy and light chains.\
   ```ab_type```: Analyzes antibody types annotated by the SAAINT parser.\
   ```HL_inf_res_num```: Analyzes and plots the number of heavy chain-light chain interface residues per antibody types.\
   ```HL_chain_len```: Analyzes and plots the length distribution of heavy and light chains.\
   ```radius```: Analyzes and plots the distribution of mean radius for scFv and VHVL types.\
   ```ag_spe```: Analyzes the top antigen species.\
   ```ag_type```: Analyzes antigen types.\
   ```ab_ag_inf_res_num```: Plots a histogram of antibody-antigen interface residue counts.\
   ```cdr_inf_res_num```: Plots a histogram of interface CDR residue counts.\
   ```cdr_inf_res_ratio```: Plots a histogram of the ratio of interface CDR residues to total interface antibody residues.\
   ```ag_chain_num```: Analyzes the number of antigens with varying chain counts.\
   ```num_pdbs_with_paired_vhvl```: Counts the PDB entries containing paired VH/VL chains.\
   ```num_pdbs_with_ag```: Counts the PDB entries that include an antigen.

1. Analyze antibody-antigen binding affinity data

   To make statistical analyses of the SAAINT-DB data entries, run:
   ```bash
   python3 scripts/run_saaintdb_analyzer.py saaintdb/saaintdb_affinity_rep.tsv -j <jobname>
   ```
   Here, ```<jobname>``` takes the following options:\
   ```affinity```: Plots a histogram of pKd values.\
   ```num_pdbs_with_affinity```: Counts the number of PDB entries containing antibody-antigen binding affinity data.

## Updating SAAINT-DB

1. Update mmCIF and FASTA files
   
   * Update mmCIF files
   ```bash
   sbatch scripts/sbatch_rsync_mmcifs.sh
   ```
   The SLURM script ```scripts/sbatch_rsync_mmcifs.sh``` runs ```scripts/run_rsync_mmcifs.py``` to update mmCIF files from PDBe. The output file ```rsync_mmCIF.out``` logs newly added, updated, and deleted PDB entries.

   * Update FASTA files
   ```bash
   sbatch scripts/sbatch_update_fastas.sh
   ```
   The SLURM script ```scripts/sbatch_update_fastas.sh``` runs ```scripts/run_update_fastas.py``` to update FASTA files. Specifically, it downloads FASTA files for newly added and updated PDB entries and remove FASTA files for deleted entries. The script generates two output files: ```database/list_obsolete_cifs.txt``` that records deleted PDB entries;  ```database/list_update_cifs.txt``` that records newly added and updated PDB entries. Additionally, the script deletes SAAINT-parser results for deleted PDB entries.

1. Run SAAINT-parser to process updated PDB entries
   
   To process updated mmCIF files, run:
   ```bash
   python3 ./scripts/run_submit_saaint_parser_jobs.py -list <mmCIF_list> <work_dir> <n_cpus>
   ```
   ```<mmCIF_list>```: List of mmCIF files to be processed, e.g., ```database/list_update_cifs.txt```. \
   ```<work_dir>```: Directory for saving intermediate outputs and final results generated by SAAINT-parser, e.g., ```database/saaint_divided```. \
   ```<n_cpus>```: Number of CPU cores allocated for processing, e.g., ```5``` (A small number of CPUs is recommended for updating).

1. Collect SAAINT-parser results to update SAAINT-DB

   To update the SAAINT-DB databases, run:
   ```bash
   python3 ./scripts/run_saaintdb_builder.py
   ```

## SAAINT-DB vs. SAbDab

   Below we provide a head-to-head statistical comparison between SAAINT-DB (updated on May 1, 2025) and SAbDab (updated on May 2, 2025) 

Item for comparison                            | SAAINT-DB    | SAbDab
--------------------------------------------:  | -----------: | -------:
#{PDB entries}                                 |  9,757       |  9,521
#{data entries}                                | 19,128       | 18,744
#{PDB entries with ≥1 paired VH/VL}            |  7,822       |  7,622
#{data entries with ≥1 paired VH/VL}           | 14,747       | 14,420
#{PDB entries with Ag}                         |  7,489       |  7,788
#{data entries with Ag}                        | 14,316       | 14,869
#{PDB entries with Ab-Ag binding affinity data}|  1,331       |    736
#{nonredundant Ab-Ag binding affinity data}    |  1,444       |    736

## Reference

[Huang X, Zhou J, Chen S, Xia X, Chen YE, Xu J. Fast, accurate parsing of antibody structures and antibody-antigen interactions enables a comprehensive structural antibody database. bioRxiv. doi: 10.1101/2025.02.25.640196](https://www.biorxiv.org/content/10.1101/2025.02.25.640196v1.full)

## Contact

Please report any bugs or suggestions to [Xiaoqiang Huang](xiaoqiah@umich.edu)
