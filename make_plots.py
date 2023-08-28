import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_percentage_by_pathogenicity(data, feature):
    """Get the percentage of each instance of feature, in each pathogenicity group
    Example: feature = 'secondary_structure' ->
    pathogenicity  secondary_structure percentage
    benign         alpha                 0.000000
                    beta                  0.000000
                     loop                  0.000000
                    turn                  0.000000
    pathogenic     alpha                 0.000000
                    beta                  0.000000
                    loop                  0.000000
                    turn                  0.000000

    Usage: get_percentage_by_pathogenicity(data_file, 'secondary_structure')
    """
    df = pd.read_csv(data)
    # split the dataset according to pathogenicity
    df_pathogenic = df[df['pathogenicity'] == "pathogenic"]
    df_non_pathogenic = df[df['pathogenicity'] == "benign"]

    # Get the percentage of each instance of feature, in each pathogenicity group
    a = df_pathogenic.groupby(feature).count() / len(df_pathogenic)
    b = df_non_pathogenic.groupby(feature).count() / len(df_non_pathogenic)
    print(a)
    print(b)


if __name__ == "__main__":
    data_file = "C:\\Users\\InbarBlech\\OneDrive - mail.tau.ac.il\\Documents\\Thesis\\Findings\\features.csv"
    # Transform to data frame
