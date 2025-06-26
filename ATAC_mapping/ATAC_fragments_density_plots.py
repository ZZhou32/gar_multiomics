import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

# ========================================================================
# CONFIGURATION (Updated with your paths)
# ========================================================================
input_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_bam/fragments"
output_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_bam/fragments/plots"
combined_plot_name = "combined_density_plots.png"
max_columns = 3                      # Number of columns in combined plot
dpi = 300                            # Image resolution
min_frag_length = 0                  # Minimum fragment length to plot
max_frag_length = 1500               # Maximum fragment length to plot

# Modern style configuration
sns.set_style("white")
plt.rcParams.update({
    'grid.color': '0.8',
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'axes.edgecolor': '0.3',
    'axes.linewidth': 0.5,
    'xtick.color': '0.3',
    'ytick.color': '0.3',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans']
})

# ========================================================================
# PROCESS FILES AND CREATE PLOTS
# ========================================================================
try:
    os.makedirs(output_folder, exist_ok=True)
except PermissionError:
    print(f"ERROR: Cannot create output directory {output_folder}")
    exit(1)

all_samples = []
processed_files = 0

print(f"Processing files in: {input_folder}")
print(f"Saving outputs to: {output_folder}")

for filename in sorted(os.listdir(input_folder)):
    if filename.endswith(".txt"):
        filepath = os.path.join(input_folder, filename)
        try:
            # Load data with error checking
            if os.path.getsize(filepath) == 0:
                print(f"Skipping empty file: {filename}")
                continue
                
            data = np.loadtxt(filepath)
            if data.size == 0:
                print(f"Skipping empty data file: {filename}")
                continue
                
            # Handle both single-column and multi-column files
            if data.ndim == 1:
                lengths = data
                counts = np.ones_like(lengths)
            else:
                counts, lengths = data[:, 0], data[:, 1]
            
            sample_name = os.path.splitext(filename)[0]
            
            # Calculate density with protection against zero division
            total_counts = np.sum(counts)
            if total_counts <= 0:
                print(f"Skipping {filename}: no valid counts")
                continue
                
            if len(lengths) > 1:
                bin_width = lengths[1] - lengths[0]
            else:
                bin_width = 1
                
            density = counts / (total_counts * max(bin_width, 1))
            
            # Store for combined plot
            all_samples.append((sample_name, lengths, density))
            
            # Create individual plot
            fig, ax = plt.subplots(figsize=(8, 5))
            
            # Main plot
            ax.plot(lengths, density, color='#1f77b4', linewidth=1.5)
            ax.fill_between(lengths, density, color='#1f77b4', alpha=0.3)
            
            # Axis configuration
            ax.set_xlim(min_frag_length, max_frag_length)
            y_max = max(0.016, np.max(density)*1.1)
            ax.set_ylim(0, y_max)
            
            # Custom ticks
            ax.set_xticks(np.arange(min_frag_length, max_frag_length+1, 200))
            ax.set_yticks(np.linspace(0, y_max, 9))
            
            # Labels and title
            ax.set_xlabel("Fragment length (bp)", fontsize=12, labelpad=10)
            ax.set_ylabel("Density", fontsize=12, labelpad=10)
            ax.set_title(sample_name.replace("_", " "), fontsize=14, pad=20)
            
            # Save individual plot
            output_path = os.path.join(output_folder, f"{sample_name}_density.png")
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
            plt.close()
            processed_files += 1
            print(f"Created: {output_path}")
            
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")

# ========================================================================
# CREATE COMBINED PLOT
# ========================================================================
if len(all_samples) > 0:
    print(f"\nCreating combined plot with {len(all_samples)} samples...")
    
    n_cols = min(max_columns, len(all_samples))
    n_rows = int(np.ceil(len(all_samples) / n_cols))
    
    fig = plt.figure(figsize=(8*n_cols, 5*n_rows), dpi=dpi)
    gs = GridSpec(n_rows, n_cols, figure=fig)
    
    fig.suptitle("ATAC-seq Fragment Length Distributions", fontsize=16, y=1.02)
    
    for i, (name, lengths, density) in enumerate(all_samples):
        ax = fig.add_subplot(gs[i//n_cols, i%n_cols])
        
        # Plot with consistent styling
        ax.plot(lengths, density, color='#1f77b4', linewidth=1.5)
        ax.fill_between(lengths, density, color='#1f77b4', alpha=0.3)
        
        # Consistent axis limits
        ax.set_xlim(min_frag_length, max_frag_length)
        ax.set_ylim(0, 0.016)  # Force same Y-axis for all
        
        # Formatting
        ax.set_title(name.replace("_", " "), fontsize=12)
        ax.set_xlabel("Fragment length (bp)", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(np.arange(min_frag_length, max_frag_length+1, 200))
    
    # Hide empty subplots if any
    for j in range(i+1, n_rows*n_cols):
        fig.add_subplot(gs[j//n_cols, j%n_cols]).axis('off')
    
    plt.tight_layout()
    combined_path = os.path.join(output_folder, combined_plot_name)
    plt.savefig(combined_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Created combined plot: {combined_path}")
else:
    print("\nNo valid samples found to plot!")

print(f"\nProcessing complete. Successfully processed {processed_files} files.")