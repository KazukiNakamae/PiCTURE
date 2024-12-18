import numpy as np
import pandas as pd
import sklearn.metrics
import argparse

def calculate_metric_with_sklearn(logits: np.ndarray, labels: np.ndarray):
    if logits.ndim == 3:
        # Reshape logits to 2D if needed
        logits = logits.reshape(-1, logits.shape[-1])
    predictions = np.argmax(logits, axis=-1)
    valid_mask = labels != -100  # Exclude padding tokens (assuming -100 is the padding token ID)
    valid_predictions = predictions[valid_mask]
    valid_labels = labels[valid_mask]
    return {
        "accuracy": sklearn.metrics.accuracy_score(valid_labels, valid_predictions),
        "f1": sklearn.metrics.f1_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
        "matthews_correlation": sklearn.metrics.matthews_corrcoef(
            valid_labels, valid_predictions
        ),
        "precision": sklearn.metrics.precision_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
        "recall": sklearn.metrics.recall_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
    }

def main(file_path):
    # ファイルの読み込み
    data = pd.read_csv(file_path)

    # logitsとlabelsの抽出
    logits = data[['logits_0', 'logits_1']].values
    labels = data['labels'].values

    # 関数を実行し、結果を表示
    metrics = calculate_metric_with_sklearn(logits, labels)
    print("Metrics:", metrics)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate metrics from logits and labels in a CSV file.')
    parser.add_argument('file_path', type=str, help='Path to the CSV file containing logits and labels.')
    args = parser.parse_args()
    main(args.file_path)
