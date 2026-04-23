import argparse
import pandas as pd
from collections import defaultdict, Counter
import sys
import os

def main():
    parser = argparse.ArgumentParser(description='Generate gene usage TSV files from BCRseq data.')
    parser.add_argument('files', nargs='+', help='Input BCRseq data TSV files')
    parser.add_argument('-n', '--names', nargs='+', required=True, help='Names for each input dataset')
    parser.add_argument('-filter', dest='lengths_file', help='TSV file with seq_id and seq_length')
    parser.add_argument('-chain', required=True, choices=['heavy', 'kappa', 'lambda'], help='Chain type: heavy, kappa, or lambda')
    parser.add_argument('-raw_reads', action='store_true', help='Use raw read counts without normalization')
    parser.add_argument('-min', dest='min_usage', type=float, default=0, help='Minimum average usage percentage cutoff for genes')
    parser.add_argument('-keep_genes', nargs='+', help='List of gene names to keep even if they are non-Functional/F')

    args = parser.parse_args()

    if len(args.files) != len(args.names):
        parser.error('Number of input files and dataset names must be equal.')

    # Define chain-specific substrings
    chain_gene_substrings = {
        'heavy': ['IGHV', 'IGHD', 'IGHJ'],
        'kappa': ['IGKV', 'IGKJ'],
        'lambda': ['IGLV', 'IGLJ']
    }

    chain = args.chain.lower()
    gene_substrings = chain_gene_substrings[chain]

    # Read lengths file if provided
    if args.lengths_file:
        lengths_df = pd.read_csv(args.lengths_file, sep='\t')
        seq_length_dict = lengths_df.set_index('seq_id')['seq_length'].to_dict()
    else:
        seq_length_dict = None

    # Initialize dictionaries
    counts_raw = {}
    gene_functionality = defaultdict(list)
    per_gene_data = {'v': defaultdict(list), 'd': defaultdict(list), 'j': defaultdict(list)}  # For per-gene data

    # Keep track of filtering information
    filtering_info = []
    overall_stats = []  # For detailed statistics

    gene_types = ['v', 'd', 'j']

    for file, name in zip(args.files, args.names):
        df = pd.read_csv(file, sep='\t')

        # Ensure required columns are present
        required_columns = ['v_call', 'd_call', 'j_call',
                            'v_sequence_alignment', 'd_sequence_alignment', 'j_sequence_alignment',
                            'v_identity', 'd_identity', 'j_identity',
                            'stop_codon',
                            'v_alignment_start', 'v_alignment_end',
                            'd_alignment_start', 'd_alignment_end',
                            'j_alignment_start', 'j_alignment_end']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing columns in {file}: {missing_columns}")
            sys.exit(1)

        # Initial dataset row count
        total_initial_rows = len(df)
        filtering_info.append(f"{name}: Starting with {total_initial_rows} total rows.")

        # Filter rows where 'stop_codon' == 'F'
        df = df[df['stop_codon'] == 'F']
        rows_after_stop_codon_filter = len(df)
        filtered_stop_codon_count = total_initial_rows - rows_after_stop_codon_filter
        filtering_info.append(f"{name}: {filtered_stop_codon_count} rows filtered out due to stop codon not equal to 'F'.")

        # Process gene calls and functionalities
        for gene_type in gene_types:
            df[f'{gene_type}_gene'] = df[f'{gene_type}_call'].str.split(',').str[0]
            df[f'{gene_type}_gene_name'] = df[f'{gene_type}_gene']
            df[f'{gene_type}_functionality'] = df[f'{gene_type}_gene'].str.split('|').str[-1]

        # Now, filter the DataFrame based on sequence lengths and identity
        if seq_length_dict is not None:
            # Function to safely get sequence length
            def get_seq_length(seq):
                if isinstance(seq, str):
                    return len(seq)
                else:
                    return -1  # Return -1 for non-string or NaN values

            # Initialize masks with the same index as df
            pass_length_mask = pd.Series([False]*len(df), index=df.index)
            pass_identity_mask = pd.Series([False]*len(df), index=df.index)
            combined_gene_mask = pd.Series([False]*len(df), index=df.index)
            gene_type_masks = {}  # Store masks for each gene type

            for gene_type in gene_types:
                length_mask = df.apply(lambda row: get_seq_length(row[f'{gene_type}_sequence_alignment']) == seq_length_dict.get(row[f'{gene_type}_gene_name'], -1), axis=1)
                identity_mask = df[f'{gene_type}_identity'] == 100
                combined_mask = length_mask & identity_mask

                # Ensure the masks have the correct index
                length_mask = length_mask.astype(bool)
                length_mask.index = df.index

                identity_mask = identity_mask.astype(bool)
                identity_mask.index = df.index

                combined_mask = combined_mask.astype(bool)
                combined_mask.index = df.index

                gene_type_masks[gene_type] = combined_mask  # Store the mask for this gene type

                # Update the overall masks
                pass_length_mask |= length_mask
                pass_identity_mask |= identity_mask
                combined_gene_mask |= combined_mask

            fail_length_mask = ~pass_length_mask
            fail_identity_mask = ~pass_identity_mask
            fail_both_mask = fail_length_mask & fail_identity_mask

            num_total_rows = len(df)
            num_rows_kept = combined_gene_mask.sum()

            # Calculate number filtered out due to length mismatch only
            num_filtered_length_only = (fail_length_mask & pass_identity_mask).sum()
            # Number filtered out due to identity mismatch only
            num_filtered_identity_only = (pass_length_mask & fail_identity_mask).sum()
            # Number filtered out due to both
            num_filtered_both = (fail_length_mask & fail_identity_mask).sum()

            filtering_info.append(f"{name}: {num_rows_kept}/{num_total_rows} rows kept after filtering based on sequence lengths and identity.")
            filtering_info.append(f"    {num_filtered_length_only} rows filtered out due to length mismatch.")
            filtering_info.append(f"    {num_filtered_identity_only} rows filtered out due to identity not equal to 100.")
            filtering_info.append(f"    {num_filtered_both} rows filtered out due to both length mismatch and identity not equal to 100.")

            df_filtered = df[combined_gene_mask].copy()
            # Create a dictionary of filtered DataFrames per gene type
            df_filtered_by_gene_type = {}
            for gene_type in gene_types:
                df_filtered_by_gene_type[gene_type] = df[gene_type_masks[gene_type]].copy()
        else:
            df_filtered = df
            df_filtered_by_gene_type = {gene_type: df.copy() for gene_type in gene_types}

        # Prepare to collect stats per gene type
        gene_type_stats = {}

        # Process each gene type
        for gene_type in gene_types:
            df_gene_type_filtered = df_filtered_by_gene_type[gene_type].copy()
            initial_gene_type_rows = len(df_gene_type_filtered)

            # Calculate alignment length
            alignment_start_col = f'{gene_type}_alignment_start'
            alignment_end_col = f'{gene_type}_alignment_end'
            df_gene_type_filtered[f'{gene_type}_alignment_length'] = (df_gene_type_filtered[alignment_end_col] - df_gene_type_filtered[alignment_start_col]).abs()

            # Apply the threshold
            if gene_type == 'v':
                threshold = 280
                df_gene_type_filtered = df_gene_type_filtered[df_gene_type_filtered[f'{gene_type}_alignment_length'] > threshold]
            elif gene_type == 'd':
                threshold = 6
                df_gene_type_filtered = df_gene_type_filtered[df_gene_type_filtered[f'{gene_type}_alignment_length'] >= threshold]
            elif gene_type == 'j':
                threshold = 20
                df_gene_type_filtered = df_gene_type_filtered[df_gene_type_filtered[f'{gene_type}_alignment_length'] >= threshold]

            rows_after_alignment_filter = len(df_gene_type_filtered)
            filtered_alignment_length_count = initial_gene_type_rows - rows_after_alignment_filter
            filtering_info.append(f"{name}: {filtered_alignment_length_count} rows filtered out due to {gene_type.upper()} alignment length not meeting threshold.")

            # Collect stats
            gene_type_stats[gene_type] = {
                'initial_rows': initial_gene_type_rows,
                'rows_after_alignment_filter': rows_after_alignment_filter
            }

            # Filter genes based on chain type
            gene_substring_matches = df_gene_type_filtered[f'{gene_type}_gene_name'].apply(
                lambda x: isinstance(x, str) and any(sub in x for sub in gene_substrings)
            )

            df_gene = df_gene_type_filtered[gene_substring_matches]

            # Store functionality status
            for gene, func in zip(df_gene[f'{gene_type}_gene_name'], df_gene[f'{gene_type}_functionality']):
                gene_functionality[gene].append(func)

            # Check if 'nt_seq_count' column exists
            if 'nt_seq_count' in df_gene.columns:
                gene_counts = df_gene.groupby(f'{gene_type}_gene_name')['nt_seq_count'].sum()
            else:
                gene_counts = df_gene[f'{gene_type}_gene_name'].value_counts()

            # Add counts to the counts_raw dictionary
            counts_raw.setdefault(name, {})[gene_type] = gene_counts

            # Collect per-gene data
            for gene in df_gene[f'{gene_type}_gene_name'].unique():
                gene_rows = df_gene[df_gene[f'{gene_type}_gene_name'] == gene].copy()
                # Add the dataset name as a new column
                gene_rows.insert(0, 'dataset', name)
                per_gene_data[gene_type][gene].append(gene_rows)

        # Add gene type stats to overall stats
        overall_stats.append({
            'dataset': name,
            'total_initial_rows': total_initial_rows,
            'rows_after_stop_codon_filter': rows_after_stop_codon_filter,
            'gene_type_stats': gene_type_stats
        })

    # Prepare to collect final gene counts per gene type
    final_gene_counts = {}

    # Process each gene type
    for gene_type in gene_types:
        # Collect counts across datasets
        gene_counts_dict = {}
        for name in args.names:
            if name in counts_raw:
                gene_counts = counts_raw[name].get(gene_type, pd.Series())
                gene_counts_dict[name] = gene_counts
        if not gene_counts_dict:
            print(f"No data available for {gene_type.upper()} genes after filtering.")
            continue

        # Create DataFrame for raw counts
        df_raw = pd.DataFrame(gene_counts_dict).fillna(0)

        # Normalize counts to percentages if -raw_reads is not passed
        if not args.raw_reads:
            df_norm = df_raw.div(df_raw.sum(axis=0), axis=1) * 100
        else:
            # If raw reads are requested, normalized counts are the same as raw counts
            df_norm = df_raw.copy()

        # Combine raw and normalized counts into a single DataFrame with multi-level columns
        df_raw.columns = pd.MultiIndex.from_product([df_raw.columns, ['raw']])
        df_norm.columns = pd.MultiIndex.from_product([df_norm.columns, ['normalized']])
        df_combined = pd.concat([df_raw, df_norm], axis=1)
        df_combined.sort_index(axis=1, level=0, inplace=True)  # Sort columns by dataset name

        # For plotting, we will use the normalized counts
        df_plot = df_norm

        # Apply minimum usage cutoff if specified
        if not args.raw_reads and args.min_usage > 0:
            # Calculate average usage across datasets
            df_plot['average_usage'] = df_plot.mean(axis=1)
            # Filter genes based on min_usage
            df_plot = df_plot[df_plot['average_usage'] >= args.min_usage]
            df_combined = df_combined.loc[df_plot.index]
            # Remove 'average_usage' column
            df_plot = df_plot.drop(columns=['average_usage'])

        if df_plot.empty:
            print(f"No data available for {gene_type.upper()} genes after applying minimum usage cutoff.")
            continue

        # Determine the most common functionality status for each gene
        gene_func_status = {}
        for gene in df_plot.index:
            funcs = gene_functionality.get(gene, [])
            if funcs:
                func_counts = Counter(funcs)
                most_common_func = func_counts.most_common(1)[0][0]
            else:
                most_common_func = 'Unknown'
            gene_func_status[gene] = most_common_func

        # Add gene functionality as a column in df_combined
        df_combined['functionality'] = pd.Series(gene_func_status)

        # Determine the genes to keep
        genes_to_keep = set()
        # Functional/F genes
        for gene, func in gene_func_status.items():
            if func in ['Functional', 'F']:
                genes_to_keep.add(gene)

        # Genes specified in args.keep_genes
        if args.keep_genes:
            genes_to_keep.update(args.keep_genes)

        # Filter df_plot to include only genes in genes_to_keep
        df_plot = df_plot.loc[df_plot.index.intersection(genes_to_keep)]

        if df_plot.empty:
            print(f"No data available for {gene_type.upper()} genes after filtering non-Functional genes.")
            continue

        # Also update df_combined
        df_combined = df_combined.loc[df_plot.index]

        # Report genes kept due to args.keep_genes
        if args.keep_genes:
            genes_kept_due_to_args = set(args.keep_genes).intersection(set(df_plot.index))
            if genes_kept_due_to_args:
                filtering_info.append(f"The following non-Functional genes were kept due to user specification: {', '.join(genes_kept_due_to_args)}")

        # Sort the DataFrame
        def sort_genes(df):
            # Calculate average usage across datasets
            df['average_usage'] = df.mean(axis=1)
            # Add functionality status
            df['functionality'] = df.index.map(lambda gene: gene_func_status.get(gene, 'Unknown'))
            # Map functionality to a sorting key
            func_sort_order = {'Functional': 0, 'F': 0, 'ORF': 1, 'ORF_P': 1, 'Pseudogene': 2, 'P': 2, 'Unknown': 3}
            df['func_sort_key'] = df['functionality'].map(lambda x: func_sort_order.get(x, 3))
            # Sort by functionality and then by average usage
            df_sorted = df.sort_values(by=['func_sort_key', 'average_usage'], ascending=[True, False])
            # Keep the sorted index
            sorted_index = df_sorted.index
            # Drop helper columns
            df_sorted = df_sorted.drop(columns=['average_usage', 'functionality', 'func_sort_key'])
            return df_sorted, sorted_index

        # Sort the DataFrame
        df_sorted, sorted_index = sort_genes(df_plot.copy())

        # Also sort the combined DataFrame and save to TSV
        df_combined_sorted = df_combined.loc[sorted_index]
        output_tsv = f'{gene_type}_gene_usage.tsv'
        df_combined_sorted.to_csv(output_tsv, sep='\t')

        # Record the number of genes included in the final TSV
        num_final_genes = len(df_combined_sorted)
        filtering_info.append(f"Final number of genes in {gene_type.upper()} usage TSV: {num_final_genes}")

        # Collect final gene counts per gene type
        final_gene_counts[gene_type] = num_final_genes

        # Create directory for per-gene TSV files
        output_dir = f'{gene_type}_gene_usage'
        os.makedirs(output_dir, exist_ok=True)

        # For each gene in the final list, save per-gene data
        for gene in df_sorted.index:
            # Collect data from per_gene_data
            gene_data_frames = per_gene_data[gene_type].get(gene, [])
            if gene_data_frames:
                gene_df = pd.concat(gene_data_frames, ignore_index=True)
                # Replace problematic characters in gene name for file naming
                safe_gene_name = gene.replace('|', '_').replace('*', '_')
                gene_output_path = os.path.join(output_dir, f'{safe_gene_name}.tsv')
                gene_df.to_csv(gene_output_path, sep='\t', index=False)

    # Add overall stats to filtering_info
    for stats in overall_stats:
        name = stats['dataset']
        filtering_info.append(f"\nDataset: {name}")
        filtering_info.append(f"  Total initial rows: {stats['total_initial_rows']}")
        filtering_info.append(f"  Rows after stop codon filter: {stats['rows_after_stop_codon_filter']}")
        for gene_type in gene_types:
            if gene_type in stats['gene_type_stats']:
                gene_stats = stats['gene_type_stats'][gene_type]
                filtering_info.append(f"    {gene_type.upper()} gene:")
                filtering_info.append(f"      Initial rows: {gene_stats['initial_rows']}")
                filtering_info.append(f"      Rows after alignment length filter: {gene_stats['rows_after_alignment_filter']}")

    # Add final gene counts per gene type
    filtering_info.append("\nFinal gene counts per gene type:")
    for gene_type, count in final_gene_counts.items():
        filtering_info.append(f"  {gene_type.upper()}: {count} genes included in final usage TSV.")

    # Save the filtering info at the end
    if filtering_info:
        with open('filtering_info.txt', 'w') as f:
            for line in filtering_info:
                f.write(line + '\n')

if __name__ == '__main__':
    main()
