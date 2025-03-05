'''This script takes IgBLAST annotations that have been appended with B cell
lineage (cluster) assignments and creates a FASTA file usable as a search 
database for mapping Ig-seq peptide sequences to their respective lineages.

Use python `make_searchable.py --help` for instructions
'''

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
from datetime import datetime
import random
import string

def setup_logging(tsv_path):
    """Setup log file."""
    base_dir = os.path.dirname(tsv_path)
    out_dir = os.path.join(base_dir, "make_searchable")
    os.makedirs(out_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    logfile = os.path.join(out_dir, f'log_searchDB_{timestamp}.txt')
    logging.basicConfig(filename=logfile, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s')
    return logfile

def read_extension_sequences(file_path):
    """Read extension sequences from a file into a list."""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def collapse_and_deduplicate(df):
    """Collapse counts and remove duplicates."""
    original_count = len(df)
    df['Collapsed'] = df.groupby('sequence_aa')['nt_seq_count'].transform('sum')
    deduplicated_df = df.drop_duplicates(subset='sequence_aa')
    deduplicated_count = len(deduplicated_df)
    logging.info(f"Rows before deduplication: {original_count}, after deduplication: {deduplicated_count}")
    return deduplicated_df

def generate_unique_string(length=8):
    """Generate a unique string of characters."""
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def create_bcrseq_id(row, chain_type):
    """Create a unique bcrseq_id for each sequence."""
    unique_str = generate_unique_string()
    if chain_type == 'heavy':
        bcrseq_id = f"BCRseq_{unique_str}_{row['ClusterID']}_{row['Collapsed']}_{row['cdr3_aa']}_{row['cdr2_aa']}_{row['cdr1_aa']}_{row['v_call']}_{row['j_call']}"
    else:
        bcrseq_id = f"BCRseq_{unique_str}_{chain_type}_{row['Collapsed']}_{row['cdr3_aa']}_{row['cdr2_aa']}_{row['cdr1_aa']}_{row['v_call']}_{row['j_call']}"
    return bcrseq_id

def determine_chain_type(df):
    """
    Determine the chain type ('heavy', 'lambda', or 'kappa') by tallying gene 
    call column data.
    """
    heavy_tally = 0
    kappa_tally = 0
    lambda_tally = 0
    gene_cols = ['v_call', 'd_call', 'j_call', 'c_call']

    for col in gene_cols:
        # Process each gene call: take the first comma-separated element
        calls = df[col].dropna().apply(lambda x: x.split(',')[0])
        tally_heavy = calls.apply(lambda s: 1 if len(s) >= 3 and s[2] == 'H' else 0).sum()
        tally_kappa = calls.apply(lambda s: 1 if len(s) >= 3 and s[2] == 'K' else 0).sum()
        tally_lambda = calls.apply(lambda s: 1 if len(s) >= 3 and s[2] == 'L' else 0).sum()

        logging.info(f"{col} tally - Heavy: {tally_heavy}, Kappa: {tally_kappa}, Lambda: {tally_lambda}")
        heavy_tally += tally_heavy
        kappa_tally += tally_kappa
        lambda_tally += tally_lambda

    logging.info(f"Total tally - Heavy: {heavy_tally}, Kappa: {kappa_tally}, Lambda: {lambda_tally}")

    # Determine which chain type has the highest tally.
    if heavy_tally >= kappa_tally and heavy_tally >= lambda_tally:
        return "heavy"
    elif kappa_tally >= lambda_tally:
        return "kappa"
    else:
        return "lambda"

def write_unaltered_sequence(row, output_handle):
    """Write an unaltered sequence to the FASTA file."""
    header = row['bcrseq_id']
    record = SeqRecord(Seq(row['sequence_aa']), id=header, description="")
    SeqIO.write(record, output_handle, "fasta")

def extend_sequence(row, extensions, output_handle):
    """Extend the sequence and write each version to the FASTA file."""
    sequence_aa = row['sequence_aa']
    fwr4_index = sequence_aa.find(row['fwr4_aa']) if pd.notna(row['fwr4_aa']) else -1
    if fwr4_index != -1:
        insert_index = fwr4_index + len(row['fwr4_aa'])
        for ext in extensions:
            mod_seq = sequence_aa[:insert_index] + ext
            header = f"{row['bcrseq_id']}_{ext}"
            record = SeqRecord(Seq(mod_seq), id=header, description="")
            SeqIO.write(record, output_handle, "fasta")
    else:
        # Handle case where fwr4_aa is present but not found in sequence_aa
        logging.warning(f"fwr4_aa not found in sequence_aa for row with bcrseq_id {row['bcrseq_id']}")

def process_tsv_to_fasta(tsv_path, extensions, logfile):
    df = pd.read_csv(tsv_path, sep='\t')
    unique_sequences = collapse_and_deduplicate(df)
    
    # Generate bcrseq_id for each unique sequence
    chain_type = determine_chain_type(df)
    unique_sequences = unique_sequences.copy() # Avoid SettingWithCopyWarning
    unique_sequences['bcrseq_id'] = unique_sequences.apply(lambda row: create_bcrseq_id(row, chain_type), axis=1)
        
    # FASTA file remains in the base directory of the TSV file
    fasta_file = os.path.splitext(tsv_path)[0] + '_searchDB.fasta'
    
    # Modified TSV file will be written in the 'make_searchable' directory
    base_dir = os.path.dirname(tsv_path)
    out_dir = os.path.join(base_dir, "make_searchable")
    modified_tsv_filename = os.path.basename(os.path.splitext(tsv_path)[0] + '_BCRseqID.tsv')
    modified_tsv_path = os.path.join(out_dir, modified_tsv_filename)
    
    extended = 0
    unextended = 0

    with open(fasta_file, 'w') as output_handle:
        for _, row in unique_sequences.iterrows():
            if pd.isna(row['fwr4_aa']):
                write_unaltered_sequence(row, output_handle)
                unextended += 1
            else:
                extend_sequence(row, extensions, output_handle)
                extended += 1

    # Write the modified DataFrame to a new TSV file in the output directory
    unique_sequences.to_csv(modified_tsv_path, sep='\t', index=False)
    logging.info(f"Modified TSV file with 'Collapsed' and 'bcrseq_id' columns created at: {modified_tsv_path}")
    logging.info(f"Number of extended unique peptide sequences: {extended}")
    logging.info(f"Number of unextended unique peptide sequences: {unextended}")
    logging.info(f"FASTA file created at: {fasta_file}")

def main():
    parser = argparse.ArgumentParser(description='Create a FASTA file with sequences extended by specified sequences.')
    parser.add_argument('tsv_path', type=str, help='Path to the input TSV file. This is an IgBLAST/MiXCR-annotated, B-cell-lineage-clustered file.')
    parser.add_argument('extensions_txt', type=str, help='Path to the file with '
                        'extension sequences. A simple txt file with a single unlabeled column of '
                        'amino acid sequences to be added after transcripts that have an FR4. These '
                        'are meant to be whatever downstream amino acid sequence(s) that is/are '
                        'required for proteomic digestion but not captured by the BCR transcript '
                        'amplification primers.')
    args = parser.parse_args()

    logfile = setup_logging(args.tsv_path)
    extensions = read_extension_sequences(args.extensions_txt)
    process_tsv_to_fasta(args.tsv_path, extensions, logfile)
    logging.info(f"Processing complete. Log file created at {logfile}")

if __name__ == "__main__":
    main()
