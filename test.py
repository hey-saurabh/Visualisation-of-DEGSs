import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns

def quality_control(data, column_indices_of_interest):
    """
     Perform quality control on gene expression data.

     Parameters:
     - data (pd.DataFrame): Gene expression data.

     Returns:
     - pd.DataFrame: Summary statistics for quality control.
    """
    # Check for missing values
    missing_values = data.isnull().sum()
    print(type(missing_values))
    print(missing_values)

    subset_data = data.iloc[:, column_indices_of_interest]
    # Check data distribution for a subset of genes (you can modify this based on your dataset)
    subset_data.boxplot()
    plt.title('Distribution of Gene Expression (Subset)')
    plt.ylabel('Expression')
    plt.show(block=False)
    plt.savefig("./")
    input()

    
    # Compile summary statistics for QC
    qc_summary = pd.DataFrame({
        'Missing Values': missing_values
    })

    return qc_summary

def process(df):

    final_dict = {}

    not_included_list = ['Row.names', 'tx_id', 'gene_id', 'tx_name', 'SYMBOL', 'log10mean',
        'keep', 'stat', 'log2FC', 'pvalue', 'locfdr', 'qvalue']

    data = {}
    data['log2FC'] = list(df['log2FC'])
    data['pValue'] = list(df['pvalue'])

    first_3 = []
    second_3 = []

    for keys, values in df.items():
        if keys not in not_included_list:
            data[keys] = list(values)

    other_columns = [key for key in data.keys() if key not in ["log2FC", "pValue"]]

    significant_result=[]
    for i, (log2fc, pvalue) in enumerate(zip(data['log2FC'], data['pValue'])):
            
        # if abs(log2fc) > 1 and pvalue <= 0.5:
        first_3.extend(data[column][i] for column in other_columns[:3])
        second_3.extend(data[column][i] for column in other_columns[3:])

        final_dict['first_set'] = first_3
        final_dict['second_set'] = second_3

        result= ttest_ind(final_dict['first_set'], final_dict['second_set'])
        
        if (result.pvalue <=0.05 and abs(result.statistic)>1):
            significant_result.append({'Group 1':final_dict['first_set'], 
                                       'Group 2':final_dict['second_set'],
                                       'pvalue':result.pvalue,
                                       'log2FC': result.statistic})

        # else:
            # continue
        
    # df2 = pd.DataFrame(final_dict)
    # print(df2)
    if significant_result:
        return pd.DataFrame(significant_result)
    else:
        print("No significant results found.")
        return pd.DataFrame()


    # return df2

def create_scatterplot(data, x_column_idx, y_column_idx, title):
    """
    Create a scatter plot.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - x_column_idx (int): Column index for the x-axis.
    - y_column_idx (int): Column index for the y-axis.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    plt.scatter(data.iloc[:, x_column_idx], data.iloc[:, y_column_idx])
    plt.title(title)
    plt.xlabel(f'Genes of Subset 1')
    plt.ylabel(f'Genes of Subset 2')
    plt.show()

def create_boxplot(data,column_indices_of_interest , title):
    """
    Create a boxplot.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - x_column_idx (int): Column index for the x-axis.
    - y_column_idx (int): Column index for the y-axis.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    subset_data = data.iloc[:, column_indices_of_interest]
    plt.figure(figsize=(10, 6))
    subset_data.boxplot()
    plt.title(title)
    plt.show()

def create_pca_plot(data, column_indices_of_interest, title):

    """
    Create a PCA plot.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - column_indices_of_interest (list of int): Column indices to include in the PCA.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    features = data.iloc[:, column_indices_of_interest]
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(features)

    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    plt.figure(figsize=(10, 6))
    plt.scatter(pca_df['PC1'], pca_df['PC2'])
    plt.title(title)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.show()

def create_ma_plot(data, log2fc_column, average_expression_column, title):
    """
    Create an MA plot.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - log2fc_column (str): Column name for the log2 fold change values.
    - average_expression_column (str): Column name for the average expression values.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    plt.figure(figsize=(10, 6))
    plt.scatter(data[average_expression_column], data[log2fc_column], alpha=0.5)
    plt.title(title)
    plt.xlabel('Average Expression')
    plt.ylabel('Log2 Fold Change')
    plt.show()

def create_pvalue_histogram(data, pvalue_column_idx, group1_column_indices, group2_column_indices, group1_label, group2_label, title):
    """
    Create a histogram of p-values between two sets of groups defined by column indices.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - pvalue_column_idx (int): Index of the column for the p-values.
    - group1_column_indices (list of int): Indices indicating the columns for the first group.
    - group2_column_indices (list of int): Indices indicating the columns for the second group.
    - group1_label (str): Label for the first group.
    - group2_label (str): Label for the second group.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    group1_pvalues = data.iloc[:, pvalue_column_idx]
    for index in group1_column_indices:
        group1_pvalues *= data.iloc[:, index]

    group2_pvalues = data.iloc[:, pvalue_column_idx]
    for index in group2_column_indices:
        group2_pvalues *= data.iloc[:, index]

    plt.figure(figsize=(10, 6))
    plt.hist([group1_pvalues, group2_pvalues], bins=30, edgecolor='black', alpha=0.7, label=[group1_label, group2_label])
    plt.title(title)
    plt.xlabel('P Values')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

def create_histogram(data, subset1_column_indices, subset2_column_indices, subset1_label, subset2_label, title):
    """
    Create histograms for two subsets defined by column indices.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - subset1_column_indices (list of int): Indices indicating the columns for the first subset.
    - subset2_column_indices (list of int): Indices indicating the columns for the second subset.
    - subset1_label (str): Label for the first subset.
    - subset2_label (str): Label for the second subset.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    subset1_data = data.iloc[:, subset1_column_indices]
    subset2_data = data.iloc[:, subset2_column_indices]

    plt.figure(figsize=(10, 6))
    plt.hist([subset1_data.values.flatten(), subset2_data.values.flatten()], bins=30, edgecolor='black', alpha=0.7, label=[subset1_label, subset2_label])
    plt.title(title)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

def create_heatmap(data, sample_indices, title):
    """
    Create a heatmap for six sample values using their indices.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data.
    - sample_indices (list of int): Indices of the samples to include in the heatmap.
    - title (str): Title of the plot.

    Returns:
    - None (displays the plot).
    """
    samples_data = data.iloc[:,sample_indices]
    plt.figure(figsize=(10, 6))
    sns.heatmap(samples_data, cmap='coolwarm', annot=True, fmt=".2f", linewidths=.5)
    plt.title(title)
    plt.xlabel('Genes')
    plt.ylabel('Samples')
    plt.show()

def main():

    df=pd.read_excel(r"./Book3modified.xlsx")

    column_indices_of_interest = [12, 13, 14, 15, 16, 17]
    # Quality control
    qc_summary = quality_control(df, column_indices_of_interest)
    print("Quality Control Summary:")
    print(qc_summary)

    # sub=process(df)
    # print(type(sub))
    # print(sub)

    # x_col_indx=[12,13,14]
    # y_col_indx=[15,16,17]
    # p_value_index=9
    # create_scatterplot(df,x_col_indx,y_col_indx,"Scatter-Plot Of Two Subsets of Genes")
    # create_boxplot(df,column_indices_of_interest,"Box Plot of Two subsets of Genes")
    # create_pca_plot(df,column_indices_of_interest,"PCA Plot")
    # create_ma_plot(df,"log2FC","log10mean", "MA Plot")
    # create_pvalue_histogram(df,p_value_index,x_col_indx,y_col_indx,"First subset of Genes","Second subset of Genes","P-Value Histogram for the Two Subsets")
    # create_histogram(df,x_col_indx,y_col_indx,"First subset of Genes","Second subset of Genes","Histogram for the Two Subsets")
    # create_heatmap(df,column_indices_of_interest,"Heatmap For the Samples")
try:
    main()
except ValueError as e:
    print(f"Error: {e}")







# def differential_expression_analysis(data, group1_columns_idx, group2_columns_idx, pvalue_column_idx, log2FC_column_idx):
#     significant_results = []
#     for idx1 in group1_columns_idx:
#         for idx2 in group2_columns_idx:
#             subset1 = data.iloc[:, idx1]
#             subset2 = data.iloc[:, idx2]
            
#             result = ttest_ind(subset1, subset2)
            

#             if result.pvalue <= 0.05 and abs(result.statistic) > 1:
#                 significant_results.append({
#                     'Group1_Column_Index': idx1,
#                     'Group2_Column_Index': idx2,
#                     'P-Value': result.pvalue,
#                     'Log2FC': result.statistic
#                 })
#     if significant_results:
#         return pd.DataFrame(significant_results)
#     else:
#         print("No significant results found.")
#         return pd.DataFrame()


"""
# # Specify the column indices for the analysis
# condition_columns_indices2 = [(12, 15), (13, 16), (14, 17)]  # Adjust based on your actual column indices
# log2FC_column_index = 8  # Adjust based on your actual column index
# pvalue_column_index = 9  # Adjust based on your actual column index

# # Perform differential expression analysis
# differentially_expressed_genes = differential_expression_analysis(
#     t, condition_columns_indices2, log2FC_column_index, pvalue_column_index)
# """
# significant_results_df = differential_expression_analysis(df, 
#                                                         group1_columns_idx=[12,13,14], 
#                                                         group2_columns_idx=[15,16,17],
#                                                         pvalue_column_idx=9,
#                                                         log2FC_column_idx=8)
