#!/usr/bin/env python3
'''
python merge_igseq_bcrseq.py --help

TO DO's:
    -those peptides that match a set of BCR seq translated seqs: check if the BCR seq
        translated seqs are subsequences of one another
    //-annotate where the matches occur: CDRH3? other?
    -implement PSM, AvgM, and PGA thresholds? is this redundant with PD2.5 filters?
    //-allow Rs & Ks to be interpreted as the same
    //-limit filter to 1 match per ClusterID rather than 1 match period
    //-create summary text file
    -Figure out if we should remove duplicates based on "PSMs Peptide ID"
'''

import argparse
import pandas as pd
import csv
import os
import re

def load_data(psm_files, bcrseq_files):
    """
    Loads and concatenates multiple PSM and BCRseq files.
    
    For PSM files:
        - Each file is read with tab-separated values, quoting enabled, and all columns are stripped of surrounding quotes.
        - The resulting dataframes are concatenated into one.
    
    For BCRseq files:
        - Each file is read with tab-separated values and columns are stripped of surrounding quotes.
        - Dataframes are concatenated using union of columns (to account for files with differing column sets).
    """
    print("Loading data...")
    # Load and concatenate all PSM files
    psm_list = []
    for file in psm_files:
        df = pd.read_csv(file, sep='\t', quoting=csv.QUOTE_ALL, dtype=str)
        df.columns = df.columns.str.strip('"')
        psm_list.append(df)
    pep_data = pd.concat(psm_list, ignore_index=True)
    print(f"Combined PSM file row count: {pep_data.shape[0]}")
    
    # Load and concatenate all BCRseq files (with potential differences in columns)
    bcrseq_list = []
    for file in bcrseq_files:
        df = pd.read_csv(file, sep='\t', dtype=str)
        df.columns = df.columns.str.strip('"')
        bcrseq_list.append(df)
    bcrseq_data = pd.concat(bcrseq_list, ignore_index=True, sort=False)
    print(f"Combined BCRseq file row count: {bcrseq_data.shape[0]}")
    
    return pep_data, bcrseq_data

def filter_igseq_data(pep_data):
    # Report initial count
    initial_count = pep_data.shape[0]
    
    # Ensure 'Protein Accessions' is string
    pep_data['Protein Accessions'] = pep_data['Protein Accessions'].astype(str)
    
    # Split 'Protein Accessions', strip whitespace and quotes, remove empty strings
    pep_data['accession_list'] = pep_data['Protein Accessions'].str.split(';').apply(
        lambda x: [acc.strip().strip('"') for acc in x if acc.strip().strip('"')]
    )

    # Filter out rows where any accession does not start with 'BCRseq_'
    def starts_with_bcrseq(accessions):
        return all(acc.startswith('BCRseq_') for acc in accessions)
    pep_data['valid_accessions'] = pep_data['accession_list'].apply(starts_with_bcrseq)
    pep_data = pep_data[pep_data['valid_accessions']]
    after_bcrseq_filter_count = pep_data.shape[0]
    print(f"IgSeq data count after filtering out accessions not starting with 'BCRseq_': {after_bcrseq_filter_count}")
    
    ## Why did I remove this step? Seems important to identify "informative" peptides
    ## Specifically, consider filtering out peptides that map to >1 CLUSTER (not just >1 BCRseq sequence), specifically in the CDR3 
    # **Commented out the filter that excludes rows with multiple accessions**
    # This allows us to retain peptides that map to multiple BCRseq sequences
    # pep_data = pep_data[pep_data['accession_list'].apply(len) == 1]
    # after_multival_filter_count = pep_data.shape[0]
    # print(f"IgSeq data count after filtering out rows with multiple accessions: {after_multival_filter_count}")
    
    # Expand the data so each row corresponds to one accession
    pep_data = pep_data.explode('accession_list').reset_index(drop=True)
    pep_data['accession'] = pep_data['accession_list']
    pep_data['accession'] = pep_data['accession'].apply(lambda x: x.split('_')[1] if '_' in x else x)
    
    return pep_data

def parse_bcrseq_data(bcrseq_data, suffixes):
    # Ensure all region sequences are strings
    region_names = ['fwr1_aa', 'cdr1_aa', 'fwr2_aa', 'cdr2_aa', 'fwr3_aa', 'cdr3_aa', 'fwr4_aa']
    for name in region_names:
        bcrseq_data[name] = bcrseq_data[name].fillna('')
    
    longest_suffix = max(suffixes, key=len)
    
    # Function to compute full_bcr_seq and region positions
    def compute_full_seq_and_positions(row):
        sequences = [row[name].upper() for name in region_names]
        lengths = [len(seq) for seq in sequences]
        # Compute cumulative positions
        positions = []
        start = 0
        for name, seq_len in zip(region_names, lengths):
            end = start + seq_len
            positions.append({'region': name, 'start': start, 'end': end})
            start = end
        # Add suffix
        suffix_len = len(longest_suffix)
        positions.append({'region': 'suffix', 'start': start, 'end': start + suffix_len})
        # Build full sequence
        full_seq = ''.join(sequences) + longest_suffix
        return pd.Series({
            'full_bcr_seq': full_seq,
            'region_positions': positions
        })
    
    bcrseq_data[['full_bcr_seq', 'region_positions']] = bcrseq_data.apply(compute_full_seq_and_positions, axis=1)
    
    print(f"Adding downstream peptide artificial extensions a.k.a. suffixes: {suffixes}")
    print(f"BCRseq data count after adding the longest downstream peptide artificial extension: {bcrseq_data.shape[0]}")
    
    # Create normalized full_bcr_seq where 'R' and 'L' are interchangeable
    bcrseq_data['normalized_full_bcr_seq'] = bcrseq_data['full_bcr_seq'].str.replace('[RL]', 'X', regex=True)
    
    # Process 'bcrseq_id' to get part after '_'
    bcrseq_data['bcrseq_id_processed'] = bcrseq_data['bcrseq_id'].apply(lambda x: x.split('_')[1] if '_' in x else x)
    
    return bcrseq_data

def merge_datasets(pep_data, bcrseq_data):
    print("Merging data on processed 'accession' and 'bcrseq_id'...")
    merged_data = pep_data.merge(bcrseq_data, left_on='accession', right_on='bcrseq_id_processed')
    print(f"Merged data count: {merged_data.shape[0]}")
    return merged_data

def find_matched_regions(merged_data):
    # Process 'Annotated Sequence' to extract peptide
    merged_data['peptide_sequence'] = merged_data['Annotated Sequence'].apply(lambda x: x.split('.')[1].upper() if '.' in x else x.upper())
    
    # Optionally, normalize sequences by replacing 'R' and 'L' with 'X'
    merged_data['normalized_peptide_sequence'] = merged_data['peptide_sequence'].str.replace('[RL]', 'X', regex=True)
    merged_data['normalized_full_bcr_seq'] = merged_data['full_bcr_seq'].str.replace('[RL]', 'X', regex=True)
    
    # For each row, find where the peptide occurs in 'full_bcr_seq', and determine the regions it overlaps
    def find_regions(row):
        normalized_full_bcr_seq = row['normalized_full_bcr_seq']
        normalized_peptide_sequence = row['normalized_peptide_sequence']
        region_positions = row['region_positions']
        # Find all occurrences of the peptide in the full sequence
        matches = [m.start() for m in re.finditer('(?={})'.format(re.escape(normalized_peptide_sequence)), normalized_full_bcr_seq)]
        if not matches:
            return pd.Series({'matched_regions': '', 'region_overlap_percentages': '', 'cdr3_coverage_percentage': 0, 'match_found': False})
        # For each match, determine regions and overlaps
        match_info_list = []
        for match_start in matches:
            match_end = match_start + len(normalized_peptide_sequence)
            peptide_length = match_end - match_start
            matched_regions = []
            overlap_percentages = {}
            cdr3_overlap_length = 0
            cdr3_length = 0
            cdr3_present = False
            for comp in region_positions:
                comp_name = comp['region']
                comp_start = comp['start']
                comp_end = comp['end']
                # Calculate overlap
                overlap_start = max(match_start, comp_start)
                overlap_end = min(match_end, comp_end)
                if overlap_start < overlap_end:  # There is an overlap
                    matched_regions.append(comp_name)
                    overlap_length = overlap_end - overlap_start
                    percentage = (overlap_length / peptide_length) * 100
                    overlap_percentages[comp_name] = percentage
                    if comp_name == 'cdr3_aa':
                        cdr3_overlap_length = overlap_length
                        cdr3_length = comp_end - comp_start
                        cdr3_present = True
            # Format the overlap percentages for output
            region_overlap_percentages = ','.join([f"{comp}:{percentage:.1f}%" for comp, percentage in overlap_percentages.items()])
            # Calculate percentage of cdr3_aa covered by the peptide
            if cdr3_present and cdr3_length > 0:
                cdr3_coverage = (cdr3_overlap_length / cdr3_length) * 100
            else:
                cdr3_coverage = 0
            match_info_list.append({
                'matched_regions': ','.join(matched_regions),
                'region_overlap_percentages': region_overlap_percentages,
                'match_start': match_start,
                'match_end': match_end,
                'cdr3_coverage_percentage': cdr3_coverage
            })
        return pd.Series({
            'match_info_list': match_info_list,
            'match_found': True
        })
    
    # Apply the function to each row
    matched_info = merged_data.apply(find_regions, axis=1)
    merged_data = pd.concat([merged_data, matched_info], axis=1)
    
    # Explode the 'match_info_list' to separate rows for each match
    matched_data = merged_data.explode('match_info_list')
    matched_data = matched_data[matched_data['match_found'] == True]
    
    # Extract the match information
    matched_data['matched_regions'] = matched_data['match_info_list'].apply(lambda x: x['matched_regions'])
    matched_data['region_overlap_percentages'] = matched_data['match_info_list'].apply(lambda x: x['region_overlap_percentages'])
    matched_data['match_start'] = matched_data['match_info_list'].apply(lambda x: x['match_start'])
    matched_data['match_end'] = matched_data['match_info_list'].apply(lambda x: x['match_end'])
    matched_data['cdr3_coverage_percentage'] = matched_data['match_info_list'].apply(lambda x: x['cdr3_coverage_percentage'])
    
    matched_data = matched_data.drop(columns=['match_info_list'])
    
    # Report matches
    print(f"Number of rows where peptide sequence was found in full_bcr_seq: {matched_data.shape[0]}")
    
    return matched_data

def filter_matches(matched_data):
    initial_match_count = matched_data.shape[0]
    
    # Identify matches that occur only in FWR or suffix regions
    def is_only_fwr_or_suffix(matched_regions):
        regions = matched_regions.split(',')
        # Exclude 'cdr' regions
        return all(not comp.startswith('cdr') for comp in regions)
    
    matched_data['is_only_fwr_or_suffix'] = matched_data['matched_regions'].apply(is_only_fwr_or_suffix)
    
    # Remove matches that occur only in FWR or suffix regions
    filtered_data = matched_data[~matched_data['is_only_fwr_or_suffix']]
    after_filter_count = filtered_data.shape[0]
    removed_count = initial_match_count - after_filter_count
    print(f"Number of matches removed that occur only in FWR or suffix regions: {removed_count}")
    print(f"Number of matches remaining after removing FWR/suffix-only matches: {after_filter_count}")
    
    return filtered_data

def finalize_data(filtered_data):
    # Now, for each peptide, determine how many unique BCRseq sequences it matches to
    peptide_match_counts = filtered_data.groupby('Annotated Sequence')['bcrseq_id'].nunique().reset_index()
    peptide_match_counts.columns = ['Annotated Sequence', 'bcrseq_match_count']
    filtered_data = filtered_data.merge(peptide_match_counts, on='Annotated Sequence', how='left')
    
    # Provide statistical distribution metrics before filtering
    print("Statistical distribution of peptide-BCRseq match counts before filtering:")
    print(peptide_match_counts['bcrseq_match_count'].describe())
    
    # Now, filter to keep peptides that match to only one BCRseq sequence
    final_data = filtered_data[filtered_data['bcrseq_match_count'] == 1]
    final_count = final_data.shape[0]
    print(f"Number of matches remaining after filtering peptides matching to only one BCRseq sequence: {final_count}")
    
    # Provide statistical distribution metrics after filtering
    print("Statistical distribution of peptide-BCRseq match counts after filtering:")
    print(final_data.groupby('Annotated Sequence')['bcrseq_id'].nunique().describe())
    
    return final_data

def main(psm_files, bcrseq_files, suffix_file=None):
    pep_data, bcrseq_data = load_data(psm_files, bcrseq_files)

    # Compute counts before filtering
    initial_pep_data_count = pep_data.shape[0]
    print(f"Initial IgSeq data count: {initial_pep_data_count}")

    # Filter IgSeq data and report counts
    pep_data = filter_igseq_data(pep_data)

    # Prepare BCRseq data
    suffixes = ['STTA', 'STTAP'] if suffix_file is None else [line.strip() for line in open(suffix_file)]
    bcrseq_data = parse_bcrseq_data(bcrseq_data, suffixes)

    # Merge the data
    merged_data = merge_datasets(pep_data, bcrseq_data)

    # Find matched regions and calculate overlap percentages
    matched_data = find_matched_regions(merged_data)

    # **New steps: Filter matches and provide statistics**
    filtered_data = filter_matches(matched_data)
    final_data = finalize_data(filtered_data)

    # Compute summary statistics
    num_filtered_matches = final_data.shape[0]
    num_total_full_bcr_seq_matched = final_data['full_bcr_seq'].shape[0]
    num_unique_full_bcr_seq_matched = final_data['full_bcr_seq'].nunique()
    num_unique_bcrseq_id_matched = final_data['bcrseq_id'].nunique()

    # Counts of matches within 'cdr3_aa' only
    matches_in_cdr3_only = final_data[final_data['matched_regions'] == 'cdr3_aa'].shape[0]

    # Counts of matches within 'cdr3_aa' plus other regions
    matches_in_cdr3_plus = final_data[final_data['matched_regions'].str.contains('cdr3_aa') & (final_data['matched_regions'] != 'cdr3_aa')].shape[0]

    # Use the first psm file's name to derive output filenames
    base_name = os.path.splitext(psm_files[0])[0]
    output_filename = base_name + '_merged.tsv'
    final_data.to_csv(output_filename, sep='\t', index=False)
    print(f"Data saved to {output_filename}")

    # Create summary file
    summary_filename = base_name + '_match_summary.txt'
    with open(summary_filename, 'w') as f:
        f.write(f"Initial IgSeq data count: {initial_pep_data_count}\n")
        f.write(f"IgSeq data count after filtering out accessions not starting with 'BCRseq_': {pep_data.shape[0]}\n")
        # f.write(f"IgSeq data count after filtering out rows with multiple accessions: {pep_data.shape[0]}\n")
        f.write(f"Merged data count: {merged_data.shape[0]}\n")
        f.write(f"Number of matches before removing FWR/suffix-only matches: {matched_data.shape[0]}\n")
        f.write(f"Number of matches after removing FWR/suffix-only matches: {filtered_data.shape[0]}\n")
        f.write(f"Number of matches remaining after filtering peptides matching to only one BCRseq sequence: {final_data.shape[0]}\n")
        f.write(f"Number of total full_bcr_seq values in final matched data: {num_total_full_bcr_seq_matched}\n")
        f.write(f"Number of unique full_bcr_seq values in final matched data: {num_unique_full_bcr_seq_matched}\n")
        f.write(f"Number of unique bcrseq_id values in final matched data: {num_unique_bcrseq_id_matched}\n")
        f.write(f"Number of matches within cdr3_aa only: {matches_in_cdr3_only}\n")
        f.write(f"Number of matches within cdr3_aa plus other regions: {matches_in_cdr3_plus}\n")
    print(f"Summary saved to {summary_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process IgSeq and BCRseq data files.")
    parser.add_argument("-psm_files", nargs='+', help="File path(s) for the peptide input TXT file(s) (output from Proteome Discoverer 2.5)", required=True)
    parser.add_argument("-bcrseq_files", nargs='+', help="File path(s) for the BCRseq TSV file(s) (already filtered for QC, IgBLAST-ed, clustered)", required=True)
    parser.add_argument("--suffix_file", help="Optional file path for downstream peptide add-ons to simulate e.g. part of the CH1 in BCRseq data (one AA seq per line). Default is the ferret CH1 needed to allow proalanine digestion.", default=None)
    args = parser.parse_args()

    main(args.psm_files, args.bcrseq_files, args.suffix_file)
