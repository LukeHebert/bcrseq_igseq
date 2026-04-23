import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils import ProtParamData
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Generate smooth histograms of CDR3 AA hydrophobicity from TSV files.")
    parser.add_argument('files', nargs='+', help='Input TSV file paths')
    parser.add_argument('-color', type=str, default='blue', help='Optional: Color for the histograms (name or hex)')
    parser.add_argument('-x_max', type=float, help='Optional: Maximum x-axis value for the plots')
    return parser.parse_args()

valid_aas = set(ProtParamData.kd.keys())

def compute_hydrophobicity(sequence):
    try:
        hydrophobicity = sum(ProtParamData.kd[aa] for aa in sequence) / len(sequence)
        return hydrophobicity
    except KeyError:
        return None

def read_data(file):
    data = pd.read_csv(file, sep='\t')
    cdr3_aa_sequences = data['cdr3_aa'].dropna()
    total_cdr3s = len(cdr3_aa_sequences)
    unique_cdr3_aa = cdr3_aa_sequences.drop_duplicates()
    unique_cdr3s = len(unique_cdr3_aa)
    # Filter sequences that contain only valid amino acids
    unique_cdr3_aa = unique_cdr3_aa[unique_cdr3_aa.apply(lambda seq: set(seq).issubset(valid_aas))]
    hydrophobicity_values = unique_cdr3_aa.apply(compute_hydrophobicity).dropna()
    return hydrophobicity_values, total_cdr3s, unique_cdr3s

def calculate_stats(values):
    stats = {
        'total': len(values),
        'mean': values.mean(),
        'median': values.median(),
        'range': (values.min(), values.max()),
        'std_dev': values.std()
    }
    return stats

def plot_histogram(data, color, alpha, label=None, x_max=None):
    # Ensure data is a DataFrame with the appropriate structure
    if isinstance(data, pd.Series):
        data = pd.DataFrame({'hydrophobicity': data})
    elif isinstance(data, list):
        data = pd.DataFrame({'hydrophobicity': data})

    sns.kdeplot(
        data=data, x='hydrophobicity', color=color, fill=True,
        alpha=alpha, label=label
    )
    plt.axvline(data['hydrophobicity'].mean(), color=color, linestyle='dashed', linewidth=1)

    # Set x-axis limits if x_max is specified
    if x_max is not None:
        plt.xlim(right=x_max)

    # Set y-axis labels to percentage
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))

def save_histogram(file_name):
    plt.xlabel('Average Hydrophobicity')
    plt.ylabel('Percentage')
    plt.legend()
    plt.savefig(file_name, dpi=800)
    plt.close()

def main(args):
    aggregate_hydrophobicity = []
    stats = []
    total_cdr3s_aggregate = 0
    unique_cdr3s_aggregate = 0

    plt.figure(figsize=(10, 6))  # Prepare figure for overlaying histograms

    # Overlay histogram for each file
    for file in args.files:
        hydrophobicity_values, total_cdr3s, unique_cdr3s = read_data(file)
        aggregate_hydrophobicity.extend(hydrophobicity_values.tolist())  # ensure list format
        plot_histogram(hydrophobicity_values, args.color, alpha=0.3, x_max=args.x_max)
        file_stats = calculate_stats(hydrophobicity_values)
        file_stats['filename'] = file
        file_stats['total_cdr3s'] = total_cdr3s
        file_stats['unique_cdr3s'] = unique_cdr3s
        stats.append(file_stats)
        total_cdr3s_aggregate += total_cdr3s
        unique_cdr3s_aggregate += unique_cdr3s

    plt.title('Overlay CDR3 AA Hydrophobicity Distribution for All Files')
    save_histogram('cdr3_hydro_overlay_histogram.png')

    # Aggregate histogram
    plt.figure(figsize=(10, 6))
    aggregate_hydrophobicity = pd.Series(aggregate_hydrophobicity)  # convert to Series before plotting
    plot_histogram(aggregate_hydrophobicity, args.color, alpha=0.5, label='Aggregate', x_max=args.x_max)

    plt.title('Aggregate CDR3 AA Hydrophobicity Distribution')
    save_histogram('cdr3_hydro_aggregate_histogram.png')

    # Calculate aggregate statistics
    aggregate_stats = calculate_stats(aggregate_hydrophobicity)
    aggregate_stats['filename'] = 'Aggregate'
    aggregate_stats['total_cdr3s'] = total_cdr3s_aggregate
    aggregate_stats['unique_cdr3s'] = unique_cdr3s_aggregate

    # Save statistics
    with open('cdr3_hydro_stats.txt', 'w') as f:
        # Save command/arguments at the top
        f.write('Command: ' + ' '.join(sys.argv) + '\n\n')

        # Save statistics for each file
        for stat in stats:
            stat_info = (
                f"{stat['filename']}: Total CDR3s = {stat['total_cdr3s']}, Unique CDR3s = {stat['unique_cdr3s']}, "
                f"Total Used = {stat['total']}, Mean = {stat['mean']:.2f}, "
                f"Median = {stat['median']:.2f}, Range = ({stat['range'][0]:.2f}, {stat['range'][1]:.2f}), "
                f"Standard Deviation = {stat['std_dev']:.2f}"
            )
            f.write(stat_info + '\n')

        # Save statistics for aggregate data
        f.write('\nAggregate Data Statistics:\n')
        aggregate_stat_info = (
            f"Aggregate: Total CDR3s = {aggregate_stats['total_cdr3s']}, Unique CDR3s = {aggregate_stats['unique_cdr3s']}, "
            f"Total Used = {aggregate_stats['total']}, Mean = {aggregate_stats['mean']:.2f}, "
            f"Median = {aggregate_stats['median']:.2f}, Range = ({aggregate_stats['range'][0]:.2f}, {aggregate_stats['range'][1]:.2f}), "
            f"Standard Deviation = {aggregate_stats['std_dev']:.2f}"
        )
        f.write(aggregate_stat_info + '\n')

if __name__ == "__main__":
    args = parse_args()
    main(args)
