# Gene Expression Analysis and Visualization
This repository contains a comprehensive pipeline for analyzing and visualizing differential gene expression data. The script performs quality control, differential expression analysis, and various visualizations to help interpret the results.

# Features
- Quality Control: Identifies missing values and duplicate entries in the dataset.
- Differential Expression Analysis: Filters genes based on log fold change and p-value thresholds.
- Visualization of DEGs: Generates various plots to visualize differentially expressed genes (DEGs), including: Volcano Plot, Scatter Plot, Boxplot, PCA Plot, MA Plot, Sample Histograms, Heatmap
- Pathway Enrichment Analysis: Identifies enriched pathways for the DEGs.
- Visualization of Pathway Enrichment: Creates visualizations for pathway enrichment scores and p-values.

# Installation

1. Clone the repository:
     ```
     git clone <repository_url>
     cd <repository_directory>

     ```
2. Install the required dependencies:
     ```
     pip install -r requirements.txt
     ```

# Usage

1. Ensure your gene expression data is in an Excel file format (.xlsx).
2. Run the script:
     ```
     python Gene\ Expression\ Analysis.py
     ```
3. Follow the prompts to input the path to your gene expression data file.


## Visualization Outputs

### Volcano Plot

Visualizes the relationship between log2 fold change and -log10(p-value) for each gene, highlighting significant DEGs.

### Scatter Plot

Displays expression values for each sample, with DEGs highlighted.

### Boxplot

Shows the distribution of expression values across different samples, highlighting significant DEGs.

### PCA Plot

Plots the principal components of the gene expression data to visualize sample clustering.

### MA Plot

Visualizes the relationship between log2 fold change and mean expression values, highlighting significant DEGs.

### Sample Histograms

Displays histograms of expression values for each sample.

### Heatmap

Creates a heatmap of DEGs for randomly selected samples.

### Pathway Enrichment Analysis

Performs pathway enrichment analysis using g and visualizes enrichment scores and p-values.

## Example Data

If you would like to test the script with example data, ensure your data file contains columns for gene identifiers, log fold change, p-values, and expression values for different samples.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue if you encounter any problems or have suggestions for improvements.

## Contact

For any issues or inquiries, please contact [Saurabh Kumar Singh] at [heyyysaurabh@gmail.com].
