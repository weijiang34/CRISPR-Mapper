import os
import logging
import hashlib
from collections import defaultdict
import concurrent.futures

logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def read_chunks(file_obj, chunk_size):
    '''
    Read the file and chunk it
    '''
    chunk = []
    for line in file_obj:
        chunk.append(line)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
    
def rename_sequences(input_file, output_file, mapping_file, max_length=50, num_workers=1, chunk_size=1000):
    """
    Preprocess the CRISPR-Spacer sequence file, and rename all names to short names to avoid makeblastdb errors. 
    Assign tasks in chunks to ensure order consistency.
    """
    def rename_chunk(chunk_lines, max_length=50):
        """
        Process a chunk of lines and return the processed lines and mapping.
        """
        def generate_short_name(long_name, length=50):
            """
            Generate a short name for a sequence using its hash value.
            """
            hash_object = hashlib.sha256(long_name.encode())
            hash_hex = hash_object.hexdigest()
            return hash_hex[:length]
        
        # Initialize lcoal variables
        name_mapping = {}
        used_names = set()
        name_counters = defaultdict(int)
        # Output lines
        processed_lines = []
        chunk_mapping_lines = []
        for line in chunk_lines:
            line = line.strip()
            if line.startswith(">"):  # sequence name
                long_name = line[1:]  # drop the '>'
                if long_name not in name_mapping:
                    # generate a unique short name
                    short_name_base = generate_short_name(long_name, max_length)
                    counter = name_counters[short_name_base]
                    unique_short_name = f"hash_{short_name_base}" if counter == 0 else f"hash_{short_name_base}_{counter}"
                    while unique_short_name in used_names:
                        counter += 1
                        unique_short_name = f"hash_{short_name_base}_{counter}"
                    used_names.add(unique_short_name)
                    name_mapping[long_name] = unique_short_name
                    name_counters[short_name_base] = counter + 1
                # add the new name to the processed lines
                processed_lines.append(f">{name_mapping[long_name]}")
                chunk_mapping_lines.append(f"{long_name}\t{name_mapping[long_name]}")
            else:
                processed_lines.append(line)
        return processed_lines, chunk_mapping_lines
    
    # initialize shared resources
    mapping_lines = []    
    if num_workers > 1:  # if the file is larger than 10MB
        logging.info(f"Processing file in parallel with {num_workers} workers...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(mapping_file, 'w') as mapfile:
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                futures = []
                for chunk in read_chunks(infile, chunk_size):
                    futures.append(executor.submit(rename_chunk, chunk, max_length))
                # get the results in order
                for future in futures:
                    processed_lines, chunk_mapping_lines = future.result()
                    outfile.write("\n".join(processed_lines) + "\n")
                    mapping_lines.extend(chunk_mapping_lines)
    else:
        logging.info("Processing file in single thread...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(mapping_file, 'w') as mapfile:
            for chunk in read_chunks(infile, chunk_size):
                processed_lines, chunk_mapping_lines = rename_chunk(chunk, max_length)
                outfile.write("\n".join(processed_lines) + "\n")
                mapping_lines.extend(chunk_mapping_lines)

    # write the mapping file
    with open(mapping_file, 'w') as mapfile:
        mapfile.write("\n".join(mapping_lines) + "\n")
        
def name_mapping_blastn_output(input_file, output_file, mapping_file, num_workers=1, chunk_size=1000):
    """
    Process a BLASTn output file, assign tasks in chunks to ensure order consistency.
    """
    def load_mapping(mapping_file):
        """
        Load the mapping file into a dictionary.
        """
        mapping = {}
        with open(mapping_file, 'r') as f:
            for line in f:
                long_name, short_name = line.strip().split('\t')
                mapping[short_name] = long_name
        return mapping

    def name_mapping_chunk(chunk_lines, mapping):
        """
        Mapping names of a chunk of lines and return the processed lines.
        """
        processed_lines = []
        for line in chunk_lines:
            line_list = line.strip().split("\t")
            try:
                line_list.insert(2, mapping[line_list[1]])
            except KeyError:
                line_list.insert(2, '')
                logging.warning(f"Failed to find the mapping for {line_list[1]}")
            processed_lines.append("\t".join(line_list))
            
        return processed_lines
    
    # load the mapping
    mapping = load_mapping(mapping_file)
    if num_workers > 1: 
        logging.info(f"Processing file in parallel with {num_workers} workers...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                futures = []
                for chunk in read_chunks(infile, chunk_size):
                    futures.append(executor.submit(name_mapping_chunk, chunk, mapping))
                # get the results in order
                for future in futures:
                    processed_lines = future.result()
                    outfile.write("\n".join(processed_lines) + "\n")
    else:
        logging.info("Processing file in single thread...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for chunk in read_chunks(infile, chunk_size):
                processed_lines = name_mapping_chunk(chunk, mapping)
                outfile.write("\n".join(processed_lines) + "\n")
    os.system(f"sed -i '1i\qseqid\tsseqid_short\tsseqid\tevalue\tpident\tlength\tslen' {output_file}")
        
if __name__=='__main__':
    pass