import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import sys
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(description='Generate gene usage heatmaps from TSV data.')
    parser.add_argument('tsv_files', nargs='+', help='TSV files generated by the first script')
    parser.add_argument('-gene_names', action='store_true', help='Display gene names as row labels')
    parser.add_argument('-color_function', action='store_true', help='Color rows based on gene functionality')
    parser.add_argument('-skinny', action='store_true', help='Adjust heatmaps so that data cells are squares')
    parser.add_argument('-chain', required=True, choices=['heavy', 'kappa', 'lambda'], help='Chain type: heavy, kappa, or lambda')
    parser.add_argument('--vmin', type=float, help='Minimum value for the colorbar range')
    parser.add_argument('--vmax', type=float, help='Maximum value for the colorbar range')
    parser.add_argument('--annotate', action='store_true', help='Annotate heatmap cells with data values')
    args = parser.parse_args()

    # Define chain-specific colormap colors
    chain_cmap_high = {
        'heavy': '#9C8307',
        'kappa': '#0A87C5',
        'lambda': '#B71C36'
    }

    # Middle color for each chain 
    chain_cmap_middle = {
        'heavy': '#CEBA5A',
        'kappa': '#4FB8EC',
        'lambda': '#CC6677'
    }

    # Define functionality color mapping
    func_color_map = {
        'Functional': '#882255',
        'F': '#882255',
        'ORF': '#117733',
        'ORF_P': '#117733',
        'Pseudogene': '#332288',
        'P': '#332288'
    }

    chain = args.chain.lower()
    cmap_low_color = 'white'
    cmap_middle_color = chain_cmap_middle[chain]
    cmap_high_color = chain_cmap_high[chain]

    # Collect all data to determine global vmin and vmax
    all_data = []
    for tsv_file in args.tsv_files:
        df_combined = pd.read_csv(tsv_file, sep='\t', index_col=0, header=[0,1])
        data_cols = [col for col in df_combined.columns if col[1]=='normalized']
        df_plot = df_combined[data_cols].copy()
        df_plot.columns = df_plot.columns.droplevel(1)
        all_data.append(df_plot.values)

    # Compute global vmin and vmax
    all_data_array = np.concatenate(all_data)
    global_vmin = np.min(all_data_array)
    global_vmax = np.max(all_data_array)

    # Override vmin and vmax if specified by the user
    if args.vmin is not None:
        global_vmin = args.vmin

    if args.vmax is not None:
        global_vmax = args.vmax

    for tsv_file in args.tsv_files:
        # Get the gene_type from the TSV filename
        basename = os.path.basename(tsv_file)
        basename_no_ext = os.path.splitext(basename)[0]
        # Replace '.tsv' with '_usage.png' for output filename
        output_png = tsv_file.replace('.tsv','.png')

        # Assume the file is named like 'v_gene_usage.tsv', 'd_gene_usage.tsv', 'j_gene_usage.tsv'
        parts = basename_no_ext.split('_')
        if len(parts) >= 1:
            gene_type = parts[0]
        else:
            gene_type = 'unknown'

        # Read the TSV file
        df_combined = pd.read_csv(tsv_file, sep='\t', index_col=0, header=[0,1])

        # Extract the normalized counts
        data_cols = [col for col in df_combined.columns if col[1]=='normalized']
        df_plot = df_combined[data_cols].copy()
        df_plot.columns = df_plot.columns.droplevel(1)  # Drop 'normalized' level

        # Get gene functionality
        if ('functionality', '') in df_combined.columns:
            func_col = df_combined[('functionality', '')]
            gene_func_status = func_col.to_dict()
        elif 'functionality' in df_combined.columns.get_level_values(0):
            func_col = df_combined['functionality']
            if isinstance(func_col, pd.DataFrame):
                gene_func_status = func_col.iloc[:,0].to_dict()
            else:
                gene_func_status = func_col.to_dict()
        else:
            gene_func_status = {gene: 'Unknown' for gene in df_plot.index}

        # Decide whether to display gene names
        yticklabels = args.gene_names

        # Function to sort DataFrame by average usage
        def sort_genes(df):
            # Calculate average usage across datasets
            df['average_usage'] = df.mean(axis=1)
            # Add functionality status
            df['functionality'] = df.index.map(lambda gene: gene_func_status.get(gene, 'Unknown'))
            # Sort by average usage
            df_sorted = df.sort_values(by=['average_usage'], ascending=[False])
            # Keep the sorted index
            sorted_index = df_sorted.index
            # Drop helper columns
            df_sorted = df_sorted.drop(columns=['average_usage', 'functionality'])
            return df_sorted, sorted_index

        # Sort the DataFrame
        df_sorted, sorted_index = sort_genes(df_plot.copy())

        # Create a list of colors for the genes in df_sorted.index
        if args.color_function:
            gene_colors = []
            for gene in df_sorted.index:
                func = gene_func_status.get(gene, 'Unknown')
                color = func_color_map.get(func, 'dimgray')  # Default to dimgray if unknown
                gene_colors.append(color)
        else:
            gene_colors = None  # No coloring

        # Adjust figure size if -skinny is specified
        if args.skinny:
            num_rows = df_sorted.shape[0]
            num_cols = df_sorted.shape[1]
            cell_size = 0.25
            width = cell_size * num_cols + 2  # Adjusted for labels and legend
            height = cell_size * num_rows + 2
            figsize = (width, height)
        else:
            figsize = None

        # Define gradient colormap with three colors
        cmap = LinearSegmentedColormap.from_list(
            'custom_cmap',
            [cmap_low_color, cmap_middle_color, cmap_high_color],
            N=256
        )

        # Create the clustermap with black borders around the cells
        clustermap = sns.clustermap(
            df_sorted,
            row_cluster=False,
            col_cluster=False,
            row_colors=gene_colors,
            yticklabels=yticklabels,
            cmap=cmap,
            figsize=figsize,
            vmin=global_vmin,
            vmax=global_vmax,
            linewidths=0.5,        # Added to create borders
            linecolor='black'      # Set border color to black
        )

        # Force row labels to be horizontal
        clustermap.ax_heatmap.set_yticklabels(
            clustermap.ax_heatmap.get_ymajorticklabels(),
            rotation=0
        )

        # Ensure that the rightmost and bottommost borders have black lines
        for _, spine in clustermap.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_edgecolor('black')

        # Draw gridlines to cover the rightmost and bottommost borders
        clustermap.ax_heatmap.hlines(
            np.arange(0, df_sorted.shape[0]+1),
            *clustermap.ax_heatmap.get_xlim(),
            color='black',
            linewidth=0.5
        )
        clustermap.ax_heatmap.vlines(
            np.arange(0, df_sorted.shape[1]+1),
            *clustermap.ax_heatmap.get_ylim(),
            color='black',
            linewidth=0.5
        )

        # Annotate heatmap cells with data values if requested
        if args.annotate:
            # Invert y-axis to match the heatmap's orientation
            ax = clustermap.ax_heatmap
            data = df_sorted.values
            for y in range(data.shape[0]):
                for x in range(data.shape[1]):
                    value = data[y, x]
                    ax.text(
                        x + 0.5,
                        y + 0.5,
                        f"{value:.1f}",
                        ha='center',
                        va='center',
                        color='black'
                    )
            # Adjust the limits to ensure annotations are inside the heatmap
            ax.set_xlim(0, data.shape[1])
            ax.set_ylim(data.shape[0], 0)

        clustermap.fig.suptitle(f'{gene_type.upper()} Gene Usage', y=1.02)

        clustermap.savefig(output_png, bbox_inches='tight', dpi=800)
        plt.close(clustermap.fig)
        print(f'Saved as {output_png}')

if __name__ == '__main__':
    main()
