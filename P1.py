import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import LabelEncoder
from skbio import diversity
from skbio.diversity import beta_diversity
from skbio.stats.ordination import PCoA
from scipy.stats import mannwhitneyu  # Import Mann-Whitney U test

def main():
    """Main function to execute the entire pipeline."""

    # Update and install prerequisites
    os.system("sudo apt upgrade")
    os.system("sudo apt update")
    os.system("sudo apt install build-essential")
    os.system("sudo apt install R-base-core")

    # Prompt for input files
    fastq1 = input("Enter the first FASTQ file: ")
    fastq2 = input("Enter the second FASTQ file: ")

    # Perform profiling with Greengenes and SILVA databases
    for database in ("G", "S"):
        output_dir = f"Results_{'Greengenes' if database == 'G' else 'SILVA'}"
        os.system(f"PM-parallel-meta -R {fastq1} {fastq2} -o {output_dir} -l 150 -D {database}")

    # Create list.txt with paths to classification.txt files
    with open("list.txt", "w") as f:
        for database in ("Greengenes", "SILVA"):
            output_path = f"Results_{database}/classification.txt"
            f.write(f"{database}\t{output_path}\n")

    # Perform taxonomic classification and functional profiling
    os.system("PM-select-taxa -l list.txt -o taxa.txt -L 5")  # Genus level
    os.system("PM-select-func -l list.txt -o func.txt -L 2")  # KEGG pathway level 2

    # Load microbiome data for diversity analysis (directly from FASTQ files)
    # Use the output files from Parallel-META profiling
    greengenes_file = f"Results_Greengenes/classification.txt"
    silva_file = f"Results_SILVA/classification.txt"

    # Load the data using appropriate methods for Parallel-META output
    # (Replace this section with the specific loading method for Parallel-META's classification.txt files)
    df_greengenes = load_greengenes_data(greengenes_file)  # Replace with actual loading
    df_silva = load_silva_data(silva_file)  # Replace with actual loading

    # Combine the data from both databases (if desired)
    df = combine_data(df_greengenes, df_silva)  # Replace with actual combination

    # Data preprocessing
    min_abundance_threshold = 10
    df_filtered = df[df.sum(axis=1) >= min_abundance_threshold]
    target_depth = 1000
    np.random.seed(42)
    sample_sizes = {sample: target_depth for sample in df_filtered.columns}
    df_rarefied, _ = df_filtered.subsample(sample_sizes, axis='sample')

    # Alpha and beta diversity analysis
    alpha_diversity = diversity.alpha_diversity('shannon', df_rarefied.values.T)
    bc_dm = beta_diversity('braycurtis', df_rarefied.values.T)

    #Principal Coordinate Analysis (PCoA) for Visualization
pcoa = PCoA(bc_dm)
pcoa.plot(df.sample_metadata, 'SampleType', cmap='viridis')

#Visualize Alpha Diversity
plt.figure(figsize=(8, 6))
plt.hist(alpha_diversity, bins=20, color='skyblue', edgecolor='black')
plt.xlabel('Alpha Diversity (Shannon Index)')
plt.ylabel('Frequency')
plt.title('Alpha Diversity Distribution')
plt.show()

#Mann-Whitney U Test 
# Load sample metadata (modify the path accordingly)
metadata_file = input("Enter the path to the sample metadata file (e.g., 'sample_metadata.csv'): ")
sample_metadata = pd.read_csv(metadata_file)

#Mann-Whitney U test comparing alpha diversity between two groups
group1_alpha = alpha_diversity[sample_metadata['Group'] == 'Control']
group2_alpha = alpha_diversity[sample_metadata['Group'] == 'Treatment']
_, alpha_p_value = mannwhitneyu(group1_alpha, group2_alpha)

#Interpretation for Alpha Diversity
if alpha_p_value < 0.05:
    print(f"Alpha diversity is significantly different between groups (p-value = {alpha_p_value:.4f}).")
else:
    print(f"No significant difference in alpha diversity between groups (p-value = {alpha_p_value:.4f}).")

#Mann-Whitney U Beta
# Assuming you have another metadata column for beta diversity comparison (e.g., 'TreatmentGroup')
beta_group1 = bc_dm[sample_metadata['Group'] == 'Control']
beta_group2 = bc_dm[sample_metadata['Group'] == 'Treatment']
_, beta_p_value = mannwhitneyu(beta_group1, beta_group2)

#Interpretation for Beta Diversity
if beta_p_value < 0.05:
    print(f"Beta diversity is significantly different between groups (p-value = {beta_p_value:.4f}).")
else:
    print(f"No significant difference in beta diversity between groups (p-value = {beta_p_value:.4f}).")

#Principal Coordinate Analysis (PCoA) for Visualization
pcoa = PCoA(bc_dm)

#Plot the PCoA for Beta Diversity
plt.figure(figsize=(8, 6))
plt.scatter(pcoa.samples['PC1'], pcoa.samples['PC2'], c=df.sample_metadata['SampleType'], cmap='viridis')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCoA Plot of Beta Diversity')
plt.colorbar(label='Sample Type')
plt.show()
if __name__
