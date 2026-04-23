import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Generate histograms of CDR3 AA lengths from TSV files.")
    parser.add_argument('files', nargs='+', help='Input TSV file paths')
    parser.add_argument('-color', type=str, default='blue', help='Optional: Color for the histograms (name or hex)')
    parser.add_argument('-y_mode', type=str, default='counts', choices=['counts', 'percentage'], 
                        help='Y-axis mode: "counts" for total lengths or "percentage" for percentage of lengths')
    parser.add_argument('-x_max', type=int, help='Optional: Maximum x-axis value for the plots')
    return parser.parse_args()

def read_data(file):
    data = pd.read_csv(file, sep='\t')
    unique_cdr3_aa = data['cdr3_aa'].dropna().drop_duplicates()
    valid_lengths = unique_cdr3_aa.apply(len)
    return valid_lengths


def calculate_stats(lengths):
    stats = {
        'total': len(lengths),
        'mean': lengths.mean(),
        'median': lengths.median(),
        'range': (lengths.min(), lengths.max()),
        'std_dev': lengths.std()
    }
    return stats

def plot_histogram(data, color, alpha, label=None, y_mode='counts', x_max=None):
    # Ensure data is a DataFrame with the appropriate structure
    if isinstance(data, pd.Series):
        data = pd.DataFrame({'lengths': data})
    elif isinstance(data, list):
        data = pd.DataFrame({'lengths': data})

    if y_mode == 'percentage':
        data['weights'] = 100*(np.ones_like(data['lengths']) / len(data['lengths']))
        weight_col = data['weights']
    else:
        weight_col = None

    # Determine bins based on x_max if provided
    max_length = data['lengths'].max()
    if x_max is not None and x_max < max_length:
        max_length = x_max

    bins = np.arange(0, max_length + 2)

    sns.histplot(
        data=data, x='lengths', bins=bins, color=color, element='step', fill=True,
        alpha=alpha, weights=weight_col, label=label
    )
    plt.axvline(data['lengths'].mean(), color=color, linestyle='dashed', linewidth=1)

    # Set x-axis limits if x_max is specified
    if x_max is not None:
        plt.xlim(left=0, right=x_max)

def save_histogram(file_name, y_mode):
    plt.xlabel('CDR3 AA Length')
    plt.ylabel('Frequency' if y_mode == 'counts' else 'Percentage')
    if y_mode == 'percentage':
        plt.legend()
    plt.savefig(file_name, dpi=800)
    plt.close()

def main(args):
    aggregate_lengths = []
    stats = []
    
    plt.figure(figsize=(10, 6))  # Prepare figure for overlaying histograms

    # Overlay histogram for each file
    for file in args.files:
        cdr3_aa_lengths = read_data(file)
        aggregate_lengths.extend(cdr3_aa_lengths.tolist())  # ensure list format
        plot_histogram(cdr3_aa_lengths, args.color, alpha=0.3, y_mode=args.y_mode, x_max=args.x_max)
        file_stats = calculate_stats(cdr3_aa_lengths)
        file_stats['filename'] = file
        stats.append(file_stats)

    plt.title('Overlay CDR3 AA Length Distribution for All Files')
    save_histogram('cdr3len_overlay_histogram.png', args.y_mode)

    # Aggregate histogram
    plt.figure(figsize=(10, 6))
    aggregate_lengths = pd.Series(aggregate_lengths)  # convert to Series before plotting
    plot_histogram(aggregate_lengths, args.color, alpha=0.5, label='Aggregate', y_mode=args.y_mode, x_max=args.x_max)

    plt.title('Aggregate CDR3 AA Length Distribution')
    save_histogram('cdr3len_aggregate_histogram.png', args.y_mode)

    # Calculate and include aggregate statistics
    aggregate_stats = calculate_stats(aggregate_lengths)
    aggregate_stats['filename'] = 'Aggregate'
    stats.append(aggregate_stats)

    # Save statistics
    with open('cdr3len_stats.txt', 'w') as f:
        # Write the command line arguments at the top
        command_line = ' '.join(sys.argv)
        f.write('Command Line: ' + command_line + '\n\n')

        # Write statistics for each file including the aggregate data
        for stat in stats:
            stat_info = (
                f"{stat['filename']}: Total = {stat['total']}, Mean = {stat['mean']:.2f}, "
                f"Median = {stat['median']}, Range = {stat['range']}, "
                f"Standard Deviation = {stat['std_dev']:.2f}"
            )
            f.write(stat_info + '\n')

if __name__ == "__main__":
    args = parse_args()
    main(args)
