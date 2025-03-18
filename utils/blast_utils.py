import os 
import logging
import subprocess

logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def makeblastdb(fasta_file, dbtype, out_db):
    '''
    Create a BLAST database from a fasta file
    '''
    # Local variables
    db_name = os.path.basename(out_db)
    # Create a BLAST database from the fasta file
    makeblastdb_command = f"makeblastdb -in {fasta_file} -dbtype {dbtype} -out {os.path.join(out_db, db_name)}"
    try:
        os.makedirs(out_db, exist_ok=True)
        subprocess.run(makeblastdb_command, check=True, text=True, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create BLAST database from {fasta_file}: {e}")
    return

def blastn(query, db_dir, out_file, threads):
    '''
    Run BLASTn to map viral contigs to CRISPR spacers
    '''
    db_name = os.path.basename(db_dir)
    db = os.path.join(db_dir, db_name)
    # Run BLASTn
    blastn_command = f"blastn -query {query} -db {db} -out {out_file} -outfmt \"6 qseqid sseqid evalue pident length slen\" -num_threads {threads}"
    
    try:
        subprocess.run(blastn_command, check=True, text=True, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to run BLASTn: {query} -> {db}: {e}")
    return

if __name__ == '__main__':
    pass