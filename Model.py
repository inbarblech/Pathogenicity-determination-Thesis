import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from imblearn.over_sampling import SMOTE
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

filename = "WFS1_with_position.csv"
data = pd.read_csv(f"gene_specific_df/{filename}")
data.head(10)
dot_index = filename.index('_')
gene = filename[:dot_index]
if gene == "combined" or gene == "combined_with_source":
    gene = "7 genes"
if gene == "features":
    gene = "190 genes"

data = pd.get_dummies(data, columns=["secondary_structure"])

mapping = {"benign": 0, "pathogenic": 1}
data["pathogenicity"] = data["pathogenicity"].map(mapping)

data = data.drop(
    labels=["uniprot_id", "stability_WT", "stability_MUT", "hydrophobicity_WT", "hydrophobicity_MUT", "volume_WT",
            "volume_MUT", "sequence_length", "oda_MUT", "oda_WT", "sasa_WT", "sasa_MUT", "RSA_MUT", "gene",
            "protein_contain_transmembrane", "is_residue_transmembranal", "aa_WT", "aa_MUT"], axis=1, inplace=False)
X = data.drop(labels="pathogenicity", axis=1, inplace=False)
y = data["pathogenicity"]

kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

tps = []
fps = []
fns = []
tns = []
errors = []
mistakes = 0

for i, (train_index, test_index) in enumerate(kf.split(X, y)):
    # Split the data into train and test sets
    print(f"Fold {i}")
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]  # Index DataFrame using iloc
    y_train, y_test = y[train_index], y[test_index]

    # drop the variant column and save it for later (for writing to csv)
    X_test_variant = X_test["variant"]
    X_train = X_train.drop(labels="variant", axis=1, inplace=False)
    X_test = X_test.drop(labels="variant", axis=1, inplace=False)

    # Apply SMOTE oversampling to the training set only
    smote = SMOTE(sampling_strategy='auto', random_state=42)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

    # class_distribution = y_train_resampled.value_counts()
    # print(f"Training set: (SMOTE)\n{class_distribution}")

    # class_distribution = y_test.value_counts()
    # print(f"Test set: \n{class_distribution}")

    # Train the model
    xgb_classifier = xgb.XGBClassifier(learning_rate=0.1, max_depth=3,
                                       n_estimators=100, random_state=42)

    xgb_classifier.fit(X_train_resampled, y_train_resampled)

    # predictions = [value for value in y_pred]

    # for pred, actual in zip(y_pred, y_test):
    #     print(f"Predicted: {pred} and Real: {actual}")
    # Make predictions on the test set
    # y_pred = model.predict(X_test)

    # Make predictions on the test set
    y_pred = xgb_classifier.predict(X_test)
    y_pred = np.round(y_pred)  # Convert predicted probabilities to binary predictions

    tp = sum((y_test == 1) & (y_pred == 1))
    fp = sum((y_test == 0) & (y_pred == 1))
    fn = sum((y_test == 1) & (y_pred == 0))
    tn = sum((y_test == 0) & (y_pred == 0))

    print(f"TP: {tp}, FP: {fp}, TN: {tn}, FN: {fn}")

    tps.append(tp)
    fps.append(fp)
    fns.append(fn)
    tns.append(tn)

    print(f"tps: {sum(tps)}, fps: {sum(fps)}, fns: {sum(fns)}, tns: {sum(tns)}")

    print(f"Prediction: {y_pred}. Reality: {y_test.values}")

    # Add the real values to the X_test for comparison
    y_test = pd.DataFrame({'pathogenicity': y_test})
    y_pred = pd.DataFrame({'predictions': y_pred})
    X_test_with_labels = pd.concat([X_test, y_test, y_pred], axis=1)
    X_test_with_labels = pd.concat([X_test_with_labels, X_test_variant], axis=1)
    print(X_test_with_labels)

    # Write X_test with labels to csv
    X_test_with_labels.to_csv(f"predictions_vs_real/{gene}_fold_{i}_prediction_vs_real.csv", index=False)

# Calculate MCC
TP = sum(tps)
FP = sum(fps)
FN = sum(fns)
TN = sum(tns)
mcc = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

print(f"TP: {sum(tps)}, FP: {sum(fps)}, TN: {sum(fns)}, FN: {sum(tns)}")

sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)
precision = TP / (TP + FP)
accuracy = (TP + TN) / (TP + TN + FP + FN)
print(f"Results for {gene}:")
print(f"Sensitivity (Recall): {sensitivity}")
print(f"Specificity: {specificity}")
print(f"Precision: {precision}")
print(f"MCC for {gene}: {mcc}")
print(f"Accuracy: {accuracy}")