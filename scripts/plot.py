import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
import numpy as np

# Load the data
data = pd.read_csv('trait_correlation_data.csv')

# Define a function to classify the clusters
trait_to_cluster = dict(zip(data['Trait'], data['Cluster']))

def assign_cluster(trait):
    # Lookup the trait in the dictionary
    return trait_to_cluster.get(trait, 'Unknown')

# Apply the function to create the new Cluster column
data['Cluster'] = data['Trait'].apply(assign_cluster)

data.to_csv('output.txt', index=False, sep='\t')


# Predefined colors for the initial clusters
initial_colors = [
    '#ffcccb',  # Light red
    '#add8e6',  # Light blue
    '#90ee90',  # Light green
    '#ffffe0',  # Light yellow
    '#dda0dd'   # Light purple
]

# Identify unique clusters and separate into known and remaining
clusters_unique = data['Cluster'].unique()
num_initial = len(initial_colors)
remaining_clusters = clusters_unique[num_initial:]

# Assign initial colors to the first clusters
cluster_colors = {cluster: color for cluster, color in zip(clusters_unique[:num_initial], initial_colors)}

# Generate dynamic colors for remaining clusters
num_remaining = len(remaining_clusters)
colormap = plt.colormaps['tab10'](num_remaining)

for i, cluster in enumerate(remaining_clusters):
    cluster_colors[cluster] = colormap(i)

# Ensure 'Unknown' has a distinct color if not already assigned
if 'Unknown' not in cluster_colors:
    cluster_colors['Unknown'] = '#d3d3d3'  # Light gray

# Function to create a horizontally oriented plot with extended x-axis limits
def plot_correlations_extended(data, horizontal=False):

    if horizontal:
        # Sort the DataFrame by Trait column in reverse order to flip horizontally
        data = data[::-1].reset_index(drop=True)

    if horizontal:
        fig, ax = plt.subplots(figsize=(18, 12))  # Wider plot
        # Plot each correlation with vertical error bars
        ax.errorbar(data.index, data['Correlation_Mania'], yerr=data['SE_Mania'], fmt='o', color='magenta', label='Correlation with Mania', capsize=5)
        ax.errorbar(data.index, data['Correlation_BD'], yerr=data['SE_BD'], fmt='^', color='blue', label='Correlation with Bipolar Disorder', capsize=5)

        # Labels and title
        #ax.set_xlabel('Trait')
        ax.set_ylabel('Correlation')
        ax.axhline(0, color='grey', linewidth=0.8)  # Add a horizontal line at correlation=0 for reference

        ax.set_xticks(data.index)
        ax.set_xticklabels(data['Trait'], ha='right', rotation=45)  # Rotate labels for better readability

        # Extend the y-axis limits to -1 to 1.1
        ax.set_ylim([-1, 1.1])

        # Color code the background for each cluster
        for cluster, color in cluster_colors.items():
            indices = data[data['Cluster'] == cluster].index
            if len(indices) > 0:
                ax.axvspan(indices[0] - 0.5, indices[-1] + 0.5, color=color, alpha=0.3)

        # Adjust layout
        plt.subplots_adjust(bottom=0.25)

    else:
        fig, ax = plt.subplots(figsize=(12, 18))  # Taller plot
        # Plot each correlation with horizontal error bars
        ax.errorbar(data['Correlation_Mania'], data.index, xerr=data['SE_Mania'], fmt='o', color='magenta', label='Correlation with Mania', capsize=5)
        ax.errorbar(data['Correlation_BD'], data.index, xerr=data['SE_BD'], fmt='^', color='blue', label='Correlation with Bipolar Disorder', capsize=5)

        # Labels and title
        ax.set_ylabel('Trait')
        ax.set_xlabel('Correlation')
        ax.set_title('Correlations of Traits with Mania and Bipolar Disorder')
        ax.axvline(0, color='grey', linewidth=0.8)  # Add a vertical line at correlation=0 for reference

        ax.set_yticks(data.index)
        ax.set_yticklabels(ha='right')  # Align y-tick labels to the right

        # Extend the x-axis limits to -1 to 1.1
        ax.set_xlim([-1, 1.1])

        # Color code the background for each cluster
        for cluster, color in cluster_colors.items():
            indices = data[data['Cluster'] == cluster].index
            if len(indices) > 0:
                ax.axhspan(indices[0] - 0.5, indices[-1] + 0.5, color=color, alpha=0.3)

        # Adjust layout
        plt.subplots_adjust(left=0.25)

    # Add a legend
    if horizontal:
        ax.legend(loc='upper left', bbox_to_anchor=(0.07, 0.1), borderaxespad=0, frameon=True)

    if horizontal == False:
        ax.legend(loc='upper left', bbox_to_anchor=(0.07, 0.9), borderaxespad=0, frameon=True)

    # Improve layout to prevent cutting off labels
    plt.tight_layout()
    plt.savefig('trait_correlation_plot_with_clusters_background.png', format='png', bbox_inches='tight', dpi=300)

    #plt.show()

# Call the function with horizontal=True or False depending on desired orientation
plot_correlations_extended(data, horizontal=True)
