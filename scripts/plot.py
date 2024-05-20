# Tested with Python 3.11.6

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import matplotlib.colors as mcolors

print("Starting Plotting")


with open('dir.txt', 'r') as f:
    working_directory = f.readline().strip()
    manifest = f.readline().strip()
    latent = f.readline().strip()

os.chdir(working_directory)  # Change the working directory to the one read from the first line

# Load the data
data = pd.read_csv('trait_correlation_data.csv')
data = data.iloc[::-1].reset_index(drop=True)
data = data[data['Trait'] != latent].reset_index(drop=True)

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
    'lightcoral',  # Light red equivalent
    'blue',   # Light blue
    'lightgreen',  # Light green
    'lightyellow', # Light yellow
    'plum'         # Light purple equivalent
]

# Convert initial colors to RGBA
initial_colors_rgba = [mcolors.to_rgba(color) for color in initial_colors]

# Identify unique clusters and separate into known and remaining
clusters_unique = data['Cluster'].unique()
num_initial = len(initial_colors)
remaining_clusters = clusters_unique[num_initial:]

# Assign initial colors to the first clusters
cluster_colors = {cluster: color for cluster, color in zip(clusters_unique[:num_initial], initial_colors_rgba)}

# Generate dynamic colors for remaining clusters
num_remaining = len(remaining_clusters)
dynamic_colors = plt.cm.tab10(np.linspace(0, 1, num_remaining))

for i, cluster in enumerate(remaining_clusters):
    cluster_colors[cluster] = dynamic_colors[i]

# Ensure 'Unknown' has a distinct color if not already assigned
if 'Unknown' not in cluster_colors:
    cluster_colors['Unknown'] = mcolors.to_rgba('grey')  # Light gray

# Function to create a horizontally oriented plot with extended x-axis limits
def plot_correlations_extended(data, horizontal=False):

    data['Trait'] = data['label_sign']

    if horizontal:
        # Sort the DataFrame by Trait column in reverse order to flip horizontally
        data = data[::-1].reset_index(drop=True)

    if horizontal:
        fig, ax = plt.subplots(figsize=(18, 12))  # Wider plot
        # Plot each correlation with vertical error bars
        ax.errorbar(data.index, data.iloc[:, 1], yerr=data.iloc[:, 3], fmt='o', color='magenta', label=f'Correlation with {latent}', capsize=5)
        ax.errorbar(data.index, data.iloc[:, 2], yerr=data.iloc[:, 4], fmt='^', color='blue', label=f'Correlation with {manifest}', capsize=5)

        # Labels and title
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
                ax.axvspan(indices[0] - 0.5, indices[-1] + 0.5, facecolor=color, alpha=0.3)

        # Adjust layout
        plt.subplots_adjust(bottom=0.25)

    else:
        fig, ax = plt.subplots(figsize=(8, 12))  # Taller plot
        # Plot each correlation with horizontal error bars
        ax.errorbar(data.iloc[:, 1], data.index, xerr=data.iloc[:, 3], fmt='o', color='magenta', label=f'Correlation with {latent}', capsize=5)
        ax.errorbar(data.iloc[:, 2], data.index, xerr=data.iloc[:, 4], fmt='^', color='blue', label=f'Correlation with {manifest}', capsize=5)

        # Labels and title
        ax.set_xlabel('Correlation')
        ax.axvline(0, color='grey', linewidth=0.8)  # Add a vertical line at correlation=0 for reference

        ax.set_yticks(data.index)
        ax.set_yticklabels(data['Trait'], ha='right', size = 12)  # Align y-tick labels to the right

        # Extend the x-axis limits to -1 to 1.1
        ax.set_xlim([-1, 1.1])

        # Color code the background for each cluster
        for cluster, color in cluster_colors.items():
            indices = data[data['Cluster'] == cluster].index
            if len(indices) > 0:
                ax.axhspan(indices[0] - 0.5, indices[-1] + 0.5, facecolor=color, alpha=0.3)

        # Adjust layout
        plt.subplots_adjust(left=0.25)

    # Add a legend
    if horizontal:
        ax.legend(loc='upper left', bbox_to_anchor=(0.07, 0.1), borderaxespad=0, frameon=True)

    if not horizontal:
        ax.legend(loc='upper left', bbox_to_anchor=(0.01, 0.975), borderaxespad=0, frameon=True)

    # Improve layout to prevent cutting off labels
    plt.tight_layout()
    plt.savefig('trait_correlation_plot_with_clusters_background.png', format='png', bbox_inches='tight', dpi=300)

    plt.show()

# Call the function with horizontal=True or False depending on desired orientation
plot_correlations_extended(data, horizontal=False)

print("Finished Plotting")
