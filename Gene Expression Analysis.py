import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from gprofiler import GProfiler

#Loading Data
def load_data(input_file):
    # Read the Excel file into a DataFrame
    df = pd.read_excel(input_file)
    print("Data Loaded Successfully!\n")
    return(df)

#Quality Check for the Data
def quality_control(dataf):

    #Information about the Dataset
    print("\nData Information\n")
    print(dataf.info())

    #Missing Values
    missing_values = dataf.isnull().sum()
    qc_summary = pd.DataFrame({'Missing Values': missing_values})
    print(qc_summary)

    #Duplicates
    duplicates = dataf[dataf.duplicated()]
    if duplicates.empty:
        print("\nNo Duplicate values present.\n")

#Differential Expression Analysis
def deg_analysis(df):

    # Select the required columns
    selected_columns = ['gene_id', 'SYMBOL', 'log10mean', 'log2FC', 'pvalue'] + list(df.columns[-6:])

    # Filter rows based on conditions (abs(log2FC) > 1 and pvalue <= 0.05)
    filtered_df = df[(abs(df['log2FC']) > 1) & (df['pvalue'] <= 0.05)]

    # Extract selected columns from the filtered DataFrame
    result_df = filtered_df[selected_columns]
    result_df.to_excel("differentially_expressed_genes_data.xlsx", index=False)
    print("\nAnalysis Completed Successfully!\n")
    print(f"\nNumber of Significant Genes found : {len(result_df)}\n")
    print(f"Results saved to: differentially_expressed_genes_data.xlsx\n")
    return result_df

#creating volcano Plot
def create_volcano_plot(dataf, fold_change_col, pvalue_col, threshold_fold_change=1, threshold_pvalue=0.05, title="Volcano Plot"):

    plt.figure(figsize=(10, 6))

    # Highlight points that meet significance thresholds
    significant_points = dataf[(abs(dataf[fold_change_col]) >= threshold_fold_change) & (dataf[pvalue_col] <= threshold_pvalue)]
    plt.scatter(significant_points[fold_change_col], -1 * (significant_points[pvalue_col].apply(lambda x: -1 * (10**-x))),
                color='red', label='Significant', alpha=0.7)

    
    plt.scatter(dataf[fold_change_col], -1 * (dataf[pvalue_col].apply(lambda x: -1 * (10**-x))),
                color='gray', label='Not Significant', alpha=0.5)

    plt.title(title)
    plt.xlabel('Log2(Fold Change)')
    plt.ylabel('-log10(P-value)')
    plt.axhline(-1 * (10**(-threshold_pvalue)), color='black', linestyle='--', linewidth=1, label=f'P-value = {threshold_pvalue}')
    plt.axvline(threshold_fold_change, color='black', linestyle='--', linewidth=1, label=f'log2Fold Change = {threshold_fold_change}')
    plt.axvline(-threshold_fold_change, color='black', linestyle='--', linewidth=1)
    plt.legend()
    plt.show()

def extract_full_dataset(dataf):
    desired_columns = ['gene_id', 'SYMBOL', 'log10mean', 'log2FC', 'pvalue']

    # Identify the last 6 columns dynamically
    additional_columns = dataf.columns[-6:].tolist()

    # Concatenate the desired and additional columns
    selected_columns = desired_columns + additional_columns

    # Extract the selected columns from the DataFrame
    extracted_df = dataf[selected_columns]
    return extracted_df

#Create ScatterPlot
def create_scatterplot(full_dataset, deg_dataset, deg_column='log2FC'):
    condition_columns = full_dataset.columns[-6:]

    # Plot scatterplots for the complete dataset
    plt.figure(figsize=(10, 6))

    for condition_column in condition_columns:
        plt.scatter(full_dataset[deg_column], full_dataset[condition_column], alpha=0.5, label=condition_column)

    # Highlight DEGs by using a different color or marker
    plt.scatter(deg_dataset[deg_column], deg_dataset[condition_columns[0]], color='red', marker='*', label='DEGs')

    plt.title('Scatter Plot for Complete Dataset with DEGs Highlighted')
    plt.xlabel(f'DEG Values ({deg_column})')
    plt.ylabel('Expression Values')
    plt.legend()
    plt.show()

#Create BoxPlot
def create_boxplot(deg_dataset):
     #Extract the last 6 columns as sample columns
    sample_columns = deg_dataset.columns[-6:]

    # Set custom colors for each sample
    sample_colors = ['skyblue', 'lightgreen', 'lightcoral', 'gold', 'lightcyan', 'lightpink']

    # Set the figure size
    plt.figure(figsize=(12, 8))

    # Create a boxplot for the sample values with custom colors and labels
    sns.boxplot(data=deg_dataset[sample_columns], palette=sample_colors, linewidth=1.5)
    
    # Add labels for each sample
    for i, col in enumerate(sample_columns):
        plt.text(i, -0.35, col, ha='center', va='center', color='black', fontweight='bold')

    plt.yscale('log')
    plt.title('Boxplot for Sample Values')
    plt.xlabel('Samples')
    plt.ylabel('Expression Values')
    plt.show()

#Create PCA Plot
def create_pca_plot(dataframe):
    # Transpose the DataFrame to have samples as rows and genes as columns
    transposed_df = dataframe.iloc[:, 5:].transpose()

    # Extract gene expression values (DEGs) and sample names
    genes = transposed_df.values  # Transposed DataFrame
    sample_names = transposed_df.index  # Sample names as index

    # Standardize the data (mean=0 and variance=1)
    genes_standardized = StandardScaler().fit_transform(genes)

    # Apply PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(genes_standardized)
    principal_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    principal_df['Sample'] = sample_names

    # Plot PCA
    plt.figure(figsize=(12, 8))
    plt.scatter(principal_df['PC1'], principal_df['PC2'], edgecolors='w', alpha=0.7, label='Samples')

    # Add labels for each sample
    for i, sample in enumerate(sample_names):
        plt.text(principal_df.loc[i, 'PC1'], principal_df.loc[i, 'PC2'], sample, fontsize=8, ha='right')

    plt.title('PCA Plot for Differentially Expressed Genes')
    plt.xlabel('Principal Component 1 (PC1)')
    plt.ylabel('Principal Component 2 (PC2)')
    plt.legend()
    plt.grid(True)
    plt.show()

#Create MA Plot
def create_ma_plot(dataframe, log2foldchange_col='log2FC', pvalue_col='pvalue', significance_threshold=0.05):
    # Extract log2 fold change and p-value columns
    log2foldchange = dataframe[log2foldchange_col]
    pvalue = dataframe[pvalue_col]

    # Create an MA plot
    plt.figure(figsize=(10, 6))
    plt.scatter(log2foldchange, -1 * (pvalue.apply(lambda x: -1 * (10 ** -x))), color='blue', alpha=0.7)

    # Highlight significant DEGs
    significant_degs = dataframe[dataframe[pvalue_col] < significance_threshold]
    plt.scatter(significant_degs[log2foldchange_col],
                -1 * (significant_degs[pvalue_col].apply(lambda x: -1 * (10 ** -x))),
                color='red', marker='*', label='Significant DEGs')

    # Add labels and title
    plt.title('MA Plot for Differentially Expressed Genes')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.legend()
    plt.show()

#Creating Sample Histogram
def plot_sample_histograms(dataframe, bins=30):
    # Extract sample columns (assuming the last six columns are samples)
    sample_columns = dataframe.columns[-6:]

    # Plot a histogram for each sample
    plt.figure(figsize=(12, 8))
    for sample_column in sample_columns:
        plt.hist(dataframe[sample_column], bins=bins, alpha=0.7, label=sample_column)

    plt.title('Histograms for Each Sample')
    plt.xlabel('Expression Values')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

#Creating Heatmap
def plot_heatmap(deg_dataset):
    sample_columns = deg_dataset.columns[5:]

    # Determine the number of samples to select (minimum of 100 or the available samples)
    num_samples_to_select = min(100, len(sample_columns))

    # Randomly select samples if available, otherwise use all samples
    selected_samples = np.random.choice(sample_columns, size=num_samples_to_select, replace=False)

    # Select the first 100 entries for each randomly selected sample
    df_subset = deg_dataset[['gene_id', 'SYMBOL', 'log10mean', 'log2FC', 'pvalue'] + list(selected_samples)].head(100)

    # Set the gene_id column as the index (assuming it's unique for each gene)
    df_subset.set_index('gene_id', inplace=True)

    # Select only the columns related to samples
    samples_data = df_subset.iloc[:, 4:]

    # Create a heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(samples_data, cmap='viridis', annot=False, linewidths=.5)
    plt.title('Heatmap of DEGs for Randomly Selected Samples')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.show()

#Visualising Data
def visualize_data(full_dataset, deg_dataset):
    create_scatterplot(full_dataset, deg_dataset)
    create_boxplot(deg_dataset)
    create_pca_plot(deg_dataset)
    create_ma_plot(full_dataset)
    plot_sample_histograms(deg_dataset)
    plot_heatmap(deg_dataset)

#Pathway Enrichment Analysis
def perform_pathway_enrichment_analysis_deg(dataset, organism='hsapiens', excel_output='enrichment_results.xlsx'):
    
    gene_list = dataset['SYMBOL'].dropna().tolist()

    # Check if the gene_list is not empty before proceeding
    if gene_list:
        # Initialize the g:Profiler client
        gp = GProfiler(return_dataframe=True)

        # Perform enrichment analysis
        enrichment_results = gp.profile(organism=organism, query=gene_list, no_evidences=False)

        # Write results to an Excel file
        enrichment_results.to_excel(excel_output, index=False)

        print(f"Enrichment results written to {excel_output} (Excel)\n")
    else:
        print("Gene list is empty or contains only 'nan' values. Check your data.")

#Determining Pathways for Genes
def pathways_for_genes(deg_dataset, organism='hsapiens', excel_output='pathways_for_genes.xlsx'):
    # Initialize the g:Profiler client
    gene_list= deg_dataset['SYMBOL'].dropna().tolist()

    gp = GProfiler(return_dataframe=True)

    # Perform enrichment analysis for the given gene list
    enrichment_results = gp.profile(organism=organism, query=gene_list, no_evidences=False)

    # Extract pathways from the results
    pathways = enrichment_results[enrichment_results['source'] == 'KEGG']

    # Check if the DataFrame is not empty before exporting to Excel
    if not pathways.empty:
        # Write results to an Excel file
        pathways.to_excel(excel_output, index=False)

        print(f"Pathways for genes written to {excel_output} (Excel)")
    else:
        print("No pathways found for the given gene list.")

#Visualizing Enrichment Scores
def visualize_enrichment_scores(file_path):
    # Read the Excel file
    df = pd.read_excel(file_path)

    # Sort the DataFrame by precision
    df_sorted = df.sort_values(by='precision', ascending=False)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.bar(df_sorted['name'], df_sorted['precision'], color='skyblue')
    plt.xticks(rotation=90)
    plt.xlabel('Pathway Name')
    plt.ylabel('Enrichment Precision')
    plt.title('Enrichment Precision for Each Pathway')
    plt.tight_layout()
    plt.show()

#P-Value ScatterPlot Based on Genes
def visualize_p_value_scatter(file_path):
    # Read the Excel file
    df = pd.read_excel(file_path)

    # Scatter plot
    plt.figure(figsize=(12, 8))
    sns.scatterplot(x='p_value', y='precision', data=df, hue='significant', palette='viridis', s=100)
    # Annotate data points with pathway names
    for _, row in df.iterrows():
        plt.annotate(row['name'], (row['p_value'], row['precision']), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('P-Value (log scale)')
    plt.ylabel('Enrichment Precision (log scale)')
    plt.title('Scatter Plot of P-Values vs. Enrichment Precision')
    plt.legend(title='Significant')
    plt.tight_layout()
    plt.show()

#Visualising Pathway
def visualise_pathways(pathways_for_genes_file_path):
    visualize_enrichment_scores(pathways_for_genes_file_path)
    visualize_p_value_scatter(pathways_for_genes_file_path)

#main
def main():
    #Step1
    #Load Data from Excel file
    file_path = input("Enter the path to the gene expression data file: ")
    dataf = load_data(file_path)

    #Step2
    #Quality Control
    print("\nPerforming Quality Control\n")
    quality_control(dataf)

    #Step3
    #Differential Expression Analysis
    print("\nPerforming Differential Expression Analysis\n")
    deg_dataset=deg_analysis(dataf)
    
    #Visulising DEGs 
    create_volcano_plot(dataf,"log2FC","pvalue")
    full_dataset=extract_full_dataset(dataf)
    
    #Step4
    #Visualize Data
    visualize_data(full_dataset,deg_dataset)

    #Step5
    #Pathway Enrichment
    print("\nPerforming Pathway Enrichment Analysis!\n")
    perform_pathway_enrichment_analysis_deg(deg_dataset)
    print("\nPerforming Gene Pathway Analysis\n")
    pathways_for_genes(deg_dataset)
    
    pathways_for_genes_file_path='./pathways_for_genes.xlsx'
    visualise_pathways(pathways_for_genes_file_path)

try:
    main()
except ValueError as e:
    print(f"Error: {e}")