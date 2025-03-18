import os 
import subprocess
import logging

logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def predict_crispr(mags_dir, mag_file, tmp_dir, dependency_dir):
    '''
    Predict CRISPR spacers from a MAG using CRT
    '''
    output_file = mag_file.split('.fa')[0] + '.crispr'
    java_command = f'java -cp {dependency_dir}/CRT1.2-CLI.jar crt {os.path.join(mags_dir, mag_file)} {tmp_dir}/{output_file}'
    try:
        subprocess.run(java_command, check=True, text=True, shell=True)
        logging.info(f"Processed {mag_file} -> {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process {mag_file}: {e}")

    parse_CRT(output_file, tmp_dir)

def parse_CRT(crispr_file, tmp_dir):
    '''
    Parse the CRT output file to extract CRISPR spacers
    '''
    # Extract CRISPR spacers from the CRT output file
    bin_header = os.path.basename(crispr_file).split('.')[0]
    with open(os.path.join(tmp_dir, crispr_file), 'r') as in_file:
        lines = in_file.readlines()
        spacers = []
        in_crispr_section = False
        in_sequence_section = False
        for i, line in enumerate(lines):
            if line.startswith('ORGANISM:'):
                cnt = 0
            else:
                # Check if we are entering CRISPR section
                if line.startswith('CRISPR '):
                    # print("break0")
                    in_crispr_section = True
                    continue
                if in_crispr_section:
                    line_list = line.split("\t")
                    # print("break1: ", line_list)
                    # Check if we are in sequence section or exiting CRISPR section
                    if line_list[0]==("--------"):
                        in_sequence_section = not in_sequence_section
                        if not in_sequence_section:
                            in_crispr_section = False
                            # print("break2")
                            continue
                    if in_sequence_section:
                        try:
                            spacer = line_list[3]
                            if spacer=='\n':
                                continue
                            spacers.append((cnt, spacer))
                            cnt += 1
                        except:
                            continue

        with open(os.path.join(tmp_dir, f"{bin_header}_CRISPR_Spacer.fa"), 'w') as out_file:
            for spacer in spacers:
                # Create spacer ID and sequence
                spacer_id = f"{bin_header}_spacer_{spacer[0]}"
                spacer_seq = spacer[1]
                out_file.write(f">{spacer_id}\n{spacer_seq}\n")
            logging.info(f"""Parsed {len(spacers)} CRISPR spacers from {crispr_file}""")
                    
    return 

if __name__ == '__main__':
    pass