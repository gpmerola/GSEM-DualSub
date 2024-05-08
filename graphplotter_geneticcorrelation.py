import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file into a pandas DataFrame
merged_df = pd.read_csv("trait_correlation_data.csv")

# Define cluster assignments (you need to define these based on your specific requirements)
clusters = ['Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Psychiatric', 'Risky Behaviour', 'Risky Behaviour', 'Risky Behaviour', 'Risky Behaviour', 'Risky Behaviour', 'Risky Behaviour', 'Social', 'Social', 'Social', 'Somatic', 'Somatic', 'Somatic', 'Somatic', 'Somatic', 'Somatic', 'Somatic', 'Somatic']
inverted_clusters = clusters[::-1]

# Assign clusters to DataFrame
merged_df['Cluster'] = inverted_clusters

# Define bar width and colors
maniacol = '#B22222'  # Firebrick
bdcolor = '#4169E1'  # RoyalBlue
bar_width = 0.6
bar_spacing = 2  # Space between individual bars within a cluster

# Dynamic adjustment for spacing between clusters
cluster_spacing = 2.0  # Increased space between clusters

# Create a list of unique clusters
unique_clusters = merged_df['Cluster'].unique()

# Adjust indices dynamically for bar plotting to create more space between pairs and clusters
cluster_indices = []
current_position = 0
for cluster in unique_clusters:
    cluster_mask = merged_df['Cluster'] == cluster
    num_traits_in_cluster = cluster_mask.sum()
    for _ in range(num_traits_in_cluster):
        cluster_indices.append(current_position)
        current_position += bar_spacing  # Adjust for the next bar
    current_position += cluster_spacing  # Extra spacing between different clusters

# Adjust figure size dynamically based on the total height needed
fig_height = len(cluster_indices) * (bar_width + bar_spacing) / 2
plt.figure(figsize=(12, 8))

# Plot bars with dynamic positioning
for i, (index, row) in enumerate(zip(cluster_indices, merged_df.itertuples())):
    plt.barh(index, row.Correlation_Mania, height=bar_width, color=maniacol, alpha=0.7, label='Mania' if i == 0 else "", edgecolor='black', linewidth=1)
    plt.barh(index + bar_width, row.Correlation_BD, height=bar_width, color=bdcolor, alpha=0.7, label='Bipolar Disorder' if i == 0 else "", edgecolor='black', linewidth=1)

    # Add error bars with adjusted positions
    plt.errorbar(merged_df['Correlation_Mania'], cluster_indices, xerr=merged_df['SE_Mania'], fmt='none', ecolor='black', capsize=2, elinewidth=0.5)
    plt.errorbar(merged_df['Correlation_BD'], [x + bar_width for x in cluster_indices], xerr=merged_df['SE_BD'], fmt='none', ecolor='black', capsize=2, elinewidth=0.5)

# Customize plot
plt.yticks([x + bar_width / 2 for x in cluster_indices], merged_df['Trait'])
plt.ylabel('Trait')
plt.xlabel('Correlation')
plt.legend()

# Adjust layout to minimize whitespace and ensure visibility
plt.tight_layout()
plt.grid(False)

# Save the plot as a PNG file
plt.savefig('trait_correlation_plot.png', format='png', bbox_inches='tight', dpi=300)

# Show the plot
plt.show()