import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score


def get_data(data_csv_path: str) -> pd.DataFrame:
    """Get the data from the csv file."""
    data = pd.read_csv(data_csv_path)
    return data


