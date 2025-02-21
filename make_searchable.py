'''This script takes IgBLAST annotations that have been appended with B cell
lineage (cluster) assignments and creates a FASTA file usable as a search 
database for mapping Ig-seq peptide sequences to their respective lineages.

USAGE:
python make_searchable.py clustered_data.tsv extension_seqs.txt
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
    directory = os.path.dirname(tsv_path)
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    logfile = os.path.join(directory, f'log_searchDB_{timestamp}.txt')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
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

def create_bcrseq_id(row):
    """Create a unique bcrseq_id for each sequence."""
    unique_str = generate_unique_string()
    bcrseq_id = f"BCRseq_{unique_str}_{row['ClusterID']}_{row['Collapsed']}_{row['cdr3_aa']}_{row['cdr2_aa']}_{row['cdr1_aa']}_{row['v_call']}_{row['j_call']}"
    return bcrseq_id

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

def process_tsv_to_fasta(tsv_path, extensions):
    df = pd.read_csv(tsv_path, sep='\t')
    unique_sequences = collapse_and_deduplicate(df)
    
    # Generate bcrseq_id for each unique sequence
    unique_sequences['bcrseq_id'] = unique_sequences.apply(create_bcrseq_id, axis=1)
    
    fasta_file = os.path.splitext(tsv_path)[0] + '_bcrseqID_searchDB.fasta'
    modified_tsv_path = os.path.splitext(tsv_path)[0] + '_with_bcrseq_id.tsv'
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

    # Write the modified DataFrame to a new TSV file
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
    'amino acid sequences to be added after transcripts that have an FR4. These'
    ' are meant to be whatever downstream amino acid sequence(s) that is/are '
    'required for proteomic digestion but not captured by the BCR transcript '
    'amplification primers.')
    args = parser.parse_args()

    logfile = setup_logging(args.tsv_path)
    extensions = read_extension_sequences(args.extensions_txt)
    process_tsv_to_fasta(args.tsv_path, extensions)
    logging.info(f"Processing complete. Log file created at {logfile}")

if __name__ == "__main__":
    main()
