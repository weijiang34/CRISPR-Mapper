import pandas as pd
import os 

def host_assign(basltn_out, out_dir, ident_threshold, coverage_threshold, eval_threshold):
    '''
    Assign hosts to the viral contigs based on the BLASTn results.
    '''
    # Read the BLASTn output
    blastn_df = pd.read_csv(basltn_out, sep='\t', header=0)
    blastn_df['coverage'] = blastn_df['length'] / blastn_df['slen']
    # Filter the results based on the identity and coverage thresholds
    blastn_df = blastn_df[(blastn_df['pident'] >= ident_threshold) & (blastn_df['coverage'] >= coverage_threshold) & (blastn_df['evalue'] <= eval_threshold)]
    blastn_df['MAG'] = blastn_df['sseqid'].str.split('_spacer_').str[0]
    blastn_df['CRISPR_Spacer'] = blastn_df['sseqid'].str.split('_spacer_').str[1]
    blastn_df = blastn_df.drop(columns=['sseqid_short', 'sseqid', 'length', 'slen']).rename({'qseqid': 'v_contig'}, axis=1)
    blastn_df = blastn_df.drop_duplicates()
    blastn_df = blastn_df.loc[:,['v_contig', 'MAG', 'CRISPR_Spacer', 'evalue', 'pident', 'coverage']]
    # Save the results
    out_file = os.path.join(out_dir, 'host_assignments.csv')
    blastn_df.to_csv(out_file, index=False)
    
    return

if __name__=='__main__':
    pass
