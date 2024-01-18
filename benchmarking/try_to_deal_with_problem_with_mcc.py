import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef


my_prediction = pd.read_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\predictions_vs_real_with_variant_all_genes_updated_201123.csv")

for gene in my_prediction["gene"].unique():
    print(f"Gene: {gene}")
    print(my_prediction[my_prediction["gene"] == gene]["pathogenicity"].value_counts())

tools = ["VEST"]

merged = my_prediction
for tool in tools:
    tool_predictions = pd.read_csv(f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\{tool}_on_dvd_data_predictions.csv")
    # Merge the two files
    merged = pd.merge(merged, tool_predictions, on=["gene", "variant"])
    # change column name of pathogenicity_x to pathogenicity
    merged = merged.rename(columns={"pathogenicity_x": "pathogenicity", "predictions_x": "predictions"})

# Print how many variants there are from each gene
merged["gene"].value_counts()
my_prediction["gene"].value_counts()


# separate the merged dataframe according to gene
genes = merged["gene"].unique()
print(f"Number of genes: {len(genes)}")

# Calculate MCC for each gene specific predictor for mutpred

# Build dictionary with gene names as keys.
mccs = {gene: 0 for gene in genes}

for tool in tools:

    gene_df = merged[merged["gene"] == "SLC26A4"]
    # Assuming you have a DataFrame called 'data' with 'prediction' and '{tool}_pathogenicity' columns

    # # Calculate MCC
    # mcc = matthews_corrcoef(gene_df['pathogenicity'], gene_df[f"{tool}_pathogenicity"])

    y_true = gene_df['pathogenicity']
    y_pred = gene_df[f"{tool}_pathogenicity"]

    # Calculate confusion matrix
    TP = np.sum((y_true == 1) & (y_pred == 1))
    FP = np.sum((y_true == 0) & (y_pred == 1))
    FN = np.sum((y_true == 1) & (y_pred == 0))
    TN = np.sum((y_true == 0) & (y_pred == 0))

    # Calculate MCC
    mcc = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    print(f"MCC: {mcc}")

    mcc = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    # Get gene name for the use for the dictionary
    gene = gene_df['gene'].unique()[0]

    # Append mcc to dictionary
    mccs[gene] = mcc

    # print(f"MCCs of {tool} predictions for each gene:")
    # for gene in mccs:
    #     print(f"{gene}: {mccs[gene]}")


print(f"MCCs of predictions for each gene:")
for gene in mccs:
    print(f"{gene}: {mccs[gene]}")

# get MCC for my predictions
mccs = {gene: 0 for gene in genes}

## Calculate for my predictions
gene_df = merged[merged["gene"] == "SLC26A4"]

# # Calculate MCC
# mcc = matthews_corrcoef(gene_df['pathogenicity'], gene_df['predictions'])

y_true = gene_df['pathogenicity']
y_pred = gene_df['predictions']

# Calculate confusion matrix
TP = np.sum((y_true == 1) & (y_pred == 1))
FP = np.sum((y_true == 0) & (y_pred == 1))
FN = np.sum((y_true == 1) & (y_pred == 0))
TN = np.sum((y_true == 0) & (y_pred == 0))

# Calculate MCC
mcc = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

# Get gene name for the use for the dictionary
gene = gene_df['gene'].unique()[0]

# Append mcc to dictionary
mccs[gene] = mcc

print(f"MCCs of predictions for each gene:")
for gene in mccs:
    print(f"{gene}: {mccs[gene]}")
