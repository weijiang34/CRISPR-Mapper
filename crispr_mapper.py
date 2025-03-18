import os 
import sys
import argparse
import subprocess
import multiprocessing
import logging

from utils.crispr_utils import predict_crispr
from utils.rename_sequences_utils import rename_sequences, name_mapping_blastn_output
from utils.blast_utils import makeblastdb, blastn
from utils.assign_utils import host_assign

logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)  

def main(mags_dir, virus_path, out_dir, threads, dependency_dir, ident_threshold, coverage_threshold, eval_threshold):
    '''
    Main function to predict CRISPR spacers from bacterial MAGs and map them to viral contigs
    '''
    # Convert all paths to absolute paths
    mags_dir = os.path.abspath(mags_dir)
    virus_path = os.path.abspath(virus_path)
    out_dir = os.path.abspath(out_dir)
    dependency_dir = os.path.abspath(dependency_dir)

    # Local varibales
    blastdb_dir = os.path.join(out_dir, 'CRISPR_Spacer_db')
    all_spacer_file = os.path.join(out_dir, 'CRISPR_Spacer.fa')
    all_spacer_with_shortname_path = os.path.join(out_dir, 'CRISPR_Spacer_shortname.fa')
    spacer_name_mapping_table_path = os.path.join(out_dir, 'CRISPR_Spacer_name_mapping.txt')
    molecule_tpye = 'nucl'
    blast_out_path = os.path.join(out_dir, 'Virus_vs_CRISPR_Spacer.shortname.blast')
    blast_out_renamed_path = os.path.join(out_dir, 'Virus_vs_CRISPR_Spacer.blast')

    # Create new folder for the output
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Create a new folder for the temporary files
    tmp_dir = os.path.join(out_dir, 'crispr_tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # PART A: CRISPR spacers related
    '''
    Predict CRISPR spacers from bacterial MAGs and merge them into a single file
    '''
    # Parse MAGs directory
    mags_dir = os.path.abspath(mags_dir)
    mag_files = [os.path.basename(f) for f in os.listdir(mags_dir) if f.endswith('.fa')]
    # Predict CRISPR spacers for each MAG in parallel
    with multiprocessing.Pool(threads) as pool:
        pool.starmap(predict_crispr, [(mags_dir, mag_file, tmp_dir, dependency_dir) for mag_file in mag_files])
    # Merge all CRISPR spacers into a single file
    crispr_files = [f for f in os.listdir(tmp_dir) if f.endswith('_CRISPR_Spacer.fa')]
    try: 
        with open(all_spacer_file, 'w') as out_file:
            for crispr_file in crispr_files:
                with open(os.path.join(tmp_dir, crispr_file), 'r') as in_file:
                    out_file.write(in_file.read())
        logging.info(f"Merged CRISPR spacers to {all_spacer_file}")
    except:
        logging.error(f"Failed to merge CRISPR spacers to {all_spacer_file}")
                
    # Rename the sequence headers to avoid BLAST errors
    logging.info(f"Renaming MAGs' names to avoid makeblastdb errors")
    try:
        rename_sequences(all_spacer_file, all_spacer_with_shortname_path, spacer_name_mapping_table_path, max_length=32, num_workers=threads, chunk_size=1000)
    except: 
        logging.error(f"Failed to rename the MAGs' names in {all_spacer_file}")
        
    # PART B: BLASTn related
    '''
    Create a BLAST database from the CRISPR spacers and run BLASTn to map viral contigs to CRISPR spacers.
    '''
    try:
        logging.info(f"Creating BLAST database from {all_spacer_with_shortname_path}")
        makeblastdb(all_spacer_with_shortname_path, molecule_tpye, blastdb_dir)
        logging.info(f"Created BLAST database: {blastdb_dir}")
        # os.remove(all_spacer_with_shortname_path)
    except:
        logging.error(f"Failed to create BLAST database from {all_spacer_with_shortname_path}")
    # Run BLASTn
    try:
        logging.info(f"Running BLASTn: {virus_path} -> {blastdb_dir}")
        blastn(virus_path, blastdb_dir, blast_out_path, threads)
        logging.info(f"BLASTn completed: {virus_path} -> {blastdb_dir}")
    except:
        logging.error(f"Failed to run BLASTn: {virus_path} -> {blastdb_dir}")
    
    # Rename the sequence headers back to the original names in the BLAST output
    try:
        logging.info(f"Renaming MAGs' names back to the original names in the BLAST output")
        name_mapping_blastn_output(blast_out_path, blast_out_renamed_path, spacer_name_mapping_table_path, num_workers=threads, chunk_size=1000)
        # os.remove(blast_out_path)
    except:
        logging.error(f"Failed to rename the MAGs' names back to the original names in the BLAST output")
        
    # PART C: Assign the viral contigs to the CRISPR spacers
    '''
    Assign hosts to the viral contigs based on the BLASTn results
    '''
    try:
        logging.info(f"Assigning hosts to the viral contigs based on the BLASTn results")
        host_assign(blast_out_renamed_path, out_dir, ident_threshold, coverage_threshold, eval_threshold)
    except:
        logging.error(f"Failed to assign hosts to the viral contigs based on the BLASTn results")

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict CRISPR spacers from bacterial MAGs and map them to viral contigs')
    parser.add_argument('-m', '--mags_dir', help='Directory containing bacterial MAGs', required=True)
    parser.add_argument('-v', '--virus', help='Path to viral contigs fasta file', required=True, default='')
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('-t', '--threads', type=int, help='Number of threads', default=1)
    parser.add_argument('-d', '--dependency', help='Dependency directory where CRT exists', default='dependency/')
    parser.add_argument('-i', '--ident', help='Identity threshold for filtering BLAST results', type=float, default=95)
    parser.add_argument('-c', '--coverage', help='Coverage (matched length over spacer length) threshold for filtering BLAST results', type=float, default=0.90)
    parser.add_argument('-e', '--eval', help='E-value threshold for filtering BLAST results', type=float, default=1e-5)
    args = parser.parse_args()
    
    main(args.mags_dir, args.virus, args.out_dir, args.threads, args.dependency, args.ident, args.coverage, args.eval)
