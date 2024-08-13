import pandas as pd

predictions_file = pd.read_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\"
                               "predictions_of_other_tools_for_all_dvd_variants.csv")

# Prediction file of the tool we want to add to the predictions file
tool_predictions = pd.read_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\AlphaMissense\\"
                               "AlphaMissense_aa_substitutions_dvd_genes.csv")

path_to_combined = "C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\predictions_of_other_tools_for_all_dvd_variants_with_alphamissense.csv"


# Add the predictions to the predictions file, go over each row in the predictions file and add the prediction
# of the tool to that variant

# # REVEL
# for index, row in predictions_file.iterrows():
#     variant = row["variant"]
#     gene = row["gene"]
#     # If there is a prediction for the variant in the REVEL file, add it to the predictions file
#     if not tool_predictions[(tool_predictions["variant"] == variant) & (tool_predictions["gene"] == gene)].empty:
#         revel_prediction = tool_predictions[(tool_predictions["variant"] == variant) & (tool_predictions["gene"] == gene)]
#         predictions_file.at[index, "revel_score"] = revel_prediction["revel_score"].values[0]
#         print("Added REVEL prediction for variant: " + variant + " in gene: " + gene)
#     else:
#         predictions_file.at[index, "revel_score"] = None
#         print("No REVEL prediction for variant: " + variant + " in gene: " + gene)
#
# predictions_file.to_csv("C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\benchmarking\\"
#                         "predictions_of_other_tools_for_all_dvd_variants_with_revel.csv", index=False)


# AlphaMissense
# Make sure the genes are only the relevant ones. Otherwise, it's very costly.
# Change the name of "protein_variant" to "variant" in the tool predictions file
tool_predictions = tool_predictions.rename(columns={'protein_variant': 'variant'})

combined_data = predictions_file.merge(tool_predictions, on=['uniprot_id', 'variant'], how='left')
combined_data.to_csv(path_to_combined, index=False)

