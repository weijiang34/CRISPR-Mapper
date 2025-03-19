# CRISPR Spacer Prediction and Host Assignment Tool

This tool predicts CRISPR spacers from bacterial MAGs (Metagenome-Assembled Genomes) and maps them to viral contigs using BLASTn. It then assigns viral contigs to bacterial hosts based on the CRISPR spacer matches. CRISPR-Spacer prediction is implemented by using CRT (CRISPR Recognition Tool).

---


## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Input Files](#input-files)
5. [Output Files](#output-files)
6. [Parameters](#parameters)
7. [Dependencies](#dependencies)
8. [License](#license)
9. [Contact](#contact)

---


## Overview

This tool performs the following steps:
1. **CRISPR Spacer Prediction**: Predicts CRISPR spacers from bacterial MAGs using the CRT: [http://www.room220.com/crt/](http://www.room220.com/crt/).
2. **BLASTn Mapping**: Maps the predicted CRISPR spacers to viral contigs using BLASTn.
3. **Host Assignment**: Assigns viral contigs to bacterial hosts based on the BLASTn results.

The tool is designed to handle large datasets efficiently by utilizing multiprocessing for parallel execution.

---


## Installation

### Prerequisites
- Python 3.7 or higher
- Biopython (`pip install biopython`)
- BLAST+ (install from [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download))
- CRT (CRISPR Recognition Tool) - Ensure it is available in the `dependency` directory.

### Installation Steps
1. Clone the repository:
   ```bash
   git clone https://github.com/weijiang34/CRISPR-Mapper.git
   cd CRISPR-Mapper
   ```
2. Install Python dependencies:
   ```bash
   # mamba/micromamba: (suggested)
   mamba create -f environment.yml
   # conda/miniconda:
   conda create -f environment.yml
   ```
3. Ensure BLAST and CRT are installed and accessible in your system or the ```dependency``` directory.
4. Test installation with: 
   ```bash
   python crispr-mapper.py --help
   ```


## Usage  
### Command-Line Interface 
```bash
python crispr_host_assignment.py \
    -m <mags_dir> \
    -v <virus_path> \
    -o <out_dir> \
    -t <threads> \
    -d <dependency_dir> \
    -i <identity_threshold> \
    -c <coverage_threshold> \
    -e <evalue_threshold>
```
### Example
```bash
python crispr_host_assignment.py \
    -m /path/to/mags \
    -v /path/to/virus_contigs.fa \
    -o /path/to/output \
    -t 8 \
    -d /path/to/dependency \
    -i 95 \
    -c 0.90 \
    -e 1e-5
```


## Input Files
### Required Inputs
1. **MAGs Directory (--mags_dir):**
   - A directory containing bacterial MAGs in FASTA format (```.fa``` or ```.fasta```).
   - Example:
   ```bash
   mags_dir/
    ├── mag1.fa
    ├── mag2.fa
    └── mag3.fa
   ```
2. **Viral Contigs File (--virus):**
   - A FASTA file containing viral contigs.
   - Default: ```dependency/```


## Output Files
The tool generates the following output files in the specified output directory (```--out_dir```):
1. **CRISPR:**
   - ```CRISPR_Spacer.fa```: Merged FASTA file of all predicted CRISPR spacers.
   - ```CRISPR_Spacer_shortname.fa```: FASTA file with renamed sequence headers for BLAST compatibility.
   - ```CRISPR_Spacer_name_mapping.txt```: Mapping table for original and renamed sequence headers.
2. **BLASTn Results:**
   - ```Virus_vs_CRISPR_Spacer.shortname.blast```: BLASTn output with renamed sequence headers.
   - ```Virus_vs_CRISPR_Spacer.blast```: BLASTn output with original sequence headers.
3. **Host Assignment Results:**
   - Host assignment results are saved in the output directory with details of viral contigs assigned to bacterial hosts.


## Parameters

| Parameter              | Description                                                                 | Default Value |
|------------------------|-----------------------------------------------------------------------------|---------------|
| `-m`, `--mags_dir`     | Directory containing bacterial MAGs in FASTA format.                        | Required      |
| `-v`, `--virus`        | Path to viral contigs FASTA file.                                           | Required      |
| `-o`, `--out_dir`      | Output directory for results.                                               | Required      |
| `-t`, `--threads`      | Number of threads for parallel processing.                                  | 1             |
| `-d`, `--dependency`   | Directory containing CRT and other dependencies.                            | `dependency/` |
| `-i`, `--ident`        | Identity threshold for filtering BLAST results (%).                        | 95            |
| `-c`, `--coverage`     | Coverage threshold (matched length / spacer length) for filtering BLAST results. | 0.90         |
| `-e`, `--eval`         | E-value threshold for filtering BLAST results.                              | 1e-5          |


## Dependencies

### Python Libraries
- `argparse`
- `os`
- `sys`
- `subprocess`
- `multiprocessing`
- `logging`
- `Biopython` (for BLAST utilities)

### External Tools
- **BLAST+**: Required for `makeblastdb` and `blastn`. Install from [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
- **CRT (CRISPR Recognition Tool)**: Ensure it is available in the `dependency` directory.


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/weijiang34/CRISPR-Mapper) or contact the author at [wjiang34-c@my.cityu.edu.hk](wjiang34-c@my.cityu.edu.hk).
