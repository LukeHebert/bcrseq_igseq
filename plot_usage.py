#!/usr/bin/env python3
"""
Takes a .tsv file from MiXCR's exportAlignments function and plots the gene
usage distribution of either V, D, or J segments as an ordered bar plot. It 
also overlays somatic hypermutation plots on top of each bar.

Example:
python plot_usage.py myseqs.tsv --segment V --tool igblast "my V gene usage title"
"""

__author__ = "Luke S Hebert"
__license__ = "MIT"

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(args):
    """ Main entry point of the app """

    '''Plot a bar graph of germline-mapped gene usage (one plot for Vs, one for
    Js, one for Ds)'''

    #user arguments
    input_file = getattr(args, 'input')
    x = getattr(args, 'segment').lower()
    tool_used = getattr(args, 'tool')
    count_max = getattr(args, 'read_max')
    shm_max = getattr(args, 'SHM_max')
    plot_title = getattr(args, 'plot_title')
    max_v = getattr(args, 'max_genes')
    
    
    #load the TSV as a dataframe
    print('\nLoading data and creating a new column...')
    if tool_used == 'MiXCR':
        cols = [f'best{x.upper()}Gene', f'{x.lower()}BestIdentityPercent']
    else: # tool_used == 'IgBLAST'
        cols = [f'{x}_call', f'{x}_identity']
    df = pd.read_table(input_file, sep='\t', usecols=cols)
    #create SHM column (the inverse percentage of % identity columns)
    if tool_used == 'MiXCR':
        df[f'{x}_shm'] = 100-(df[cols[1]]*100)
    else: # tool_used == 'IgBLAST'
        df[f'{x}_shm'] = 100-df[cols[1]]

    # Count the occurrences of each gene
    gene_counts = df[cols[0]].value_counts()
    # Add gene count to DataFrame
    df['gene_count'] = df[cols[0]].map(gene_counts)
    # Sort DataFrame by gene_count in descending order
    df = df.sort_values(by='gene_count', ascending=False)

    #get a list of ordered gene names for ordering the plot X axis
    gene_counts = df[cols[0]].value_counts().to_dict()
    #order the keys and values by most to least frequent
    gene_ordered = [[g,c] for g,c in 
        zip(list(gene_counts.keys()), list(gene_counts.values()))]
    gene_ordered = sorted(gene_ordered, key=lambda x: x[1], reverse=True)
    #below, [:bar_num+1] gets rid of the low-count bars on the tail of the plot
    x_order = [x[0] for x in gene_ordered][:max_v+1]
    #make a list of xtick label colors so that newly discovered genes stand out
    xtickcolors = ['red' if x.count('_') == 2 else 'black' for x in x_order]


    #plot the bar graph of gene usage + violin plot of SHM and save as a PNG
    sns.set_style("darkgrid")
    #create twin Y axes
    fig, ax1 = plt.subplots() # initializes figure and plots
    ax2 = ax1.twinx() # applies twinx to ax2 which is the second y axis.
    print('\nPlotting usage (gene call counts)...')
    if tool_used == 'MiXCR':
        gene_col = f'best{x.upper()}Gene'
        identity_col = f'{x.lower()}BestIdentityPercent'
    else: # tool_used == 'IgBLAST'
        gene_col = f'{x}_call'
        identity_col = f'{x}_identity'
    sns.countplot(
        x=gene_col, 
        data=df,
        order=x_order,
        color='tomato',
        ax=ax1)
    print('\nPlotting somatic hypermutation distributions...')
    sns.violinplot(
        x=gene_col, 
        y=f'{x}_shm', 
        data=df,
        order=x_order,
        cut=0, 
        color='maroon', 
        alpha=0.1, 
        linewidth=0,
        ax=ax2)
    #configure axes
    ax2.grid(None)
    ax1.set_ylabel(f'Reads Mapped to {x.upper()} Germline', color='tomato')
    ax1.set_xlabel('') #get rid of X axis label; it's implied
    ax2.set_ylabel(f'Read-Germline % Dissimilarity', color='maroon')
    ax1.set_ylim(bottom=0, top=count_max) #read-count-per-V y axis maximum
    ax2.set_ylim(bottom=0, top=shm_max) #SHM y axis maximum
    if thousands(ax1.get_yticks()): #prevent truncation of non-0s
        ylabels = [f'{int(num)}K' for num in ax1.get_yticks()/1000] #large y ticks
        ax1.set_yticklabels(ylabels) #large y ticks
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)
    for xtick, color in zip(ax1.get_xticklabels(),xtickcolors):
        xtick.set_color(color)
    for xtick, color in zip(ax2.get_xticklabels(),xtickcolors):
        xtick.set_color(color)
    ytickcolors1 = ['tomato' for tick in ax1.get_yticklabels()]
    for ytick, color in zip(ax1.get_yticklabels(),ytickcolors1):
        ytick.set_color(color)
    ytickcolors2 = ['maroon' for tick in ax2.get_yticklabels()]
    for ytick, color in zip(ax2.get_yticklabels(),ytickcolors2):
        ytick.set_color(color)
    plt.title(f'{plot_title}')
    plt.tight_layout()
    out_path = input_file.replace('.tsv',f'_usageSHM_{x.upper()}.png')
    plt.savefig(out_path, dpi=800)
    plt.close() 

    #also export SHM data as a TSV file
    # Remove duplicates based on the gene column
    df = df.drop_duplicates(subset=[gene_col])
    # Specify the filename for the exported data
    out_tsv = input_file.replace('.tsv', f'_usageSHM_{x.upper()}.tsv')
    # Select the columns to be exported
    export_data = df[[gene_col, 'gene_count', f'{x}_shm']]
    # Export the data to a .tsv file
    export_data.to_csv(out_tsv, sep='\t', index=False)

def thousands(lst):
    """
    This function takes a list of numerical values and returns True if all of them are divisible by 1000,
    otherwise it returns False.
    """
    return all(num % 1000 == 0 for num in lst if isinstance(num, (int, float)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description='Process input arguments')
    parser.add_argument('input', type=str, help='Input TSV file name')
    parser.add_argument('--segment', '-s', type=str, required=True, 
        choices=['V', 'D', 'J'], help='Specify segment type (V, D, or J)')
    parser.add_argument('--tool', '-t', type=str, required=True, 
            choices=['mixcr', 'igblast'], help='Specify the tool used to annotate/create the input (MiXCR or IgBLAST)')
    parser.add_argument('--read-max', type=int, default=250000,
        help='Maximum Y axis for read count (default: 250000)')
    parser.add_argument('--SHM-max', type=int, default=20,
        help='Maximum Y axis for somatic hypermutation (default: 20)')
    parser.add_argument('plot_title', type=str, help='Title for the plot')
    parser.add_argument('--max-genes', type=int, default=30,
        help=('Maximum number of genes (e.g. IGHV gene segments) to represent on the bar graph '
        '(default: 30)'))
    args = parser.parse_args()
    
    main(args)