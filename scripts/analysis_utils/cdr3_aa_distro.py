import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Generate amino acid usage bar plots from TSV files.")
    parser.add_argument('files', nargs='+', help='Input TSV file paths')
    parser.add_argument('-color', type=str, default='blue', help='Optional: Color for the bar plots (name or hex)')
    return parser.parse_args()

def read_data(file):
    data = pd.read_csv(file, sep='\t', low_memory=False)
    # Keep only unique cdr3_aa values after dropping NaNs
    cdr3_aa = data['cdr3_aa'].dropna()
    initial_count = len(cdr3_aa)
    unique_cdr3_aa = cdr3_aa.drop_duplicates()
    unique_count = len(unique_cdr3_aa)
    return unique_cdr3_aa.tolist(), initial_count, unique_count

def count_amino_acids(sequences):
    total_counts = {}
    for seq in sequences:
        analysis = ProteinAnalysis(seq)
        aa_counts = analysis.count_amino_acids()
        for aa, count in aa_counts.items():
            total_counts[aa] = total_counts.get(aa, 0) + count
    total_aa = sum(total_counts.values())
    return total_counts, total_aa

def normalize_counts(counts, total_aa):
    percentages = {aa: (count / total_aa) * 100 for aa, count in counts.items()}
    return percentages

def plot_amino_acid_usage(percentages_list, labels, title, filename, colors=None, show_legend=True):
    plt.figure(figsize=(12, 8))
    amino_acids = sorted(list(set().union(*[p.keys() for p in percentages_list])))
    x = np.arange(len(amino_acids))
    width = 0.8 / len(percentages_list)  # Adjust bar width based on the number of datasets

    for i, percentages in enumerate(percentages_list):
        usage = [percentages.get(aa, 0) for aa in amino_acids]
        if colors:
            color = colors[i]
        else:
            color = None
        plt.bar(x + i * width, usage, width=width, label=labels[i], color=color)

    plt.xticks(x + width * (len(percentages_list) - 1) / 2, amino_acids)
    plt.xlabel('Amino Acids')
    plt.ylabel('Percentage Usage')
    plt.title(title)
    if show_legend:
        plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=800)
    plt.close()

def plot_aggregate_amino_acid_usage(mean_percentages, std_devs, title, filename, color):
    plt.figure(figsize=(12, 8))
    amino_acids = sorted(mean_percentages.keys())
    x = np.arange(len(amino_acids))

    usage = [mean_percentages[aa] for aa in amino_acids]
    errors = [std_devs[aa] for aa in amino_acids]

    plt.bar(x, usage, yerr=errors, capsize=5, color=color)

    plt.xticks(x, amino_acids)
    plt.xlabel('Amino Acids')
    plt.ylabel('Percentage Usage')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=800)
    plt.close()

def save_percentages(percentages, filename, initial_count=None, unique_count=None, std_devs=None, command_args=None):
    with open(filename, 'w') as f:
        if command_args:
            f.write(f'Command: {" ".join(command_args)}\n\n')
        if initial_count is not None and unique_count is not None:
            f.write(f'Total initial CDR3s: {initial_count}\n')
            f.write(f'Total unique CDR3s: {unique_count}\n\n')
        if std_devs:
            f.write('Amino Acid\tPercentage\tStdDev\n')
            for aa in sorted(percentages.keys()):
                std_dev = std_devs.get(aa, 0)
                f.write(f'{aa}\t{percentages[aa]:.2f}\t{std_dev:.2f}\n')
        else:
            f.write('Amino Acid\tPercentage\n')
            for aa in sorted(percentages.keys()):
                f.write(f'{aa}\t{percentages[aa]:.2f}\n')

def main(args):
    aggregate_sequences = []
    percentages_list = []
    labels = []
    total_initial_cdr3s = 0
    total_unique_cdr3s = 0

    # Process each file individually
    for file in args.files:
        sequences, initial_count, unique_count = read_data(file)
        aggregate_sequences.extend(sequences)
        total_initial_cdr3s += initial_count
        total_unique_cdr3s += unique_count
        counts, total_aa = count_amino_acids(sequences)
        percentages = normalize_counts(counts, total_aa)
        percentages_list.append(percentages)

        # Use the basename for labels
        filename = os.path.basename(file)
        labels.append(filename)

        # Save individual amino acid percentages to a file in the current directory
        save_percentages(
            percentages,
            f'aa_distro_{filename.replace(".tsv","")}.txt',
            initial_count=initial_count,
            unique_count=unique_count
        )

    # Plot overlay of amino acid usage for all files
    plot_amino_acid_usage(
        percentages_list,
        labels,
        'Amino Acid Usage in CDR3 Sequences (Overlay)',
        'aa_distro_overlay.png'
    )

    # Calculate mean and std dev of percentages across files
    amino_acids = sorted(list(set().union(*[p.keys() for p in percentages_list])))
    mean_percentages = {}
    std_percentages = {}

    for aa in amino_acids:
        aa_percentages = [p.get(aa, 0) for p in percentages_list]
        mean_percentages[aa] = np.mean(aa_percentages)
        if len(aa_percentages) > 1:
            std_percentages[aa] = np.std(aa_percentages, ddof=1)
        else:
            std_percentages[aa] = 0.0  # If only one sample, std dev is zero

    # Save aggregate percentages along with total initial and unique CDR3 counts
    save_percentages(
        mean_percentages,
        'aa_distro_aggregate.txt',
        initial_count=total_initial_cdr3s,
        unique_count=len(set(aggregate_sequences)),
        std_devs=std_percentages,
        command_args=sys.argv  # Include command-line arguments
    )

    # Plot aggregate amino acid usage with error bars
    plot_aggregate_amino_acid_usage(
        mean_percentages,
        std_percentages,
        'Amino Acid Usage in CDR3 Sequences (Aggregate)',
        'aa_distro_aggregate.png',
        args.color
    )

if __name__ == "__main__":
    args = parse_args()
    main(args)
