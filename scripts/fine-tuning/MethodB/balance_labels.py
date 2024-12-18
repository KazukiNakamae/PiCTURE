import pandas as pd
import numpy as np
import argparse

def balance_labels(df, seed=42):
    label_counts = df['label'].value_counts()
    min_count = label_counts.min()
    balanced_df = df.groupby('label').apply(lambda x: x.sample(min_count, random_state=seed)).reset_index(drop=True)
    return balanced_df

def main(dev_path, test_path, train_path, output_dir, seed):
    # Load the CSV files
    dev_df = pd.read_csv(dev_path)
    test_df = pd.read_csv(test_path)
    train_df = pd.read_csv(train_path)

    # Balance the labels for each DataFrame
    balanced_dev_df = balance_labels(dev_df, seed)
    balanced_test_df = balance_labels(test_df, seed)
    balanced_train_df = balance_labels(train_df, seed)

    # Save the balanced DataFrames to new CSV files
    balanced_dev_df.to_csv(f'{output_dir}/dev.csv', index=False)
    balanced_test_df.to_csv(f'{output_dir}/test.csv', index=False)
    balanced_train_df.to_csv(f'{output_dir}/train.csv', index=False)

    print(f"Balanced datasets saved to {output_dir} as 'balanced_dev.csv', 'balanced_test.csv', and 'balanced_train.csv'.")
    print(f"Balanced Dev Shape: {balanced_dev_df.shape}")
    print(f"Balanced Test Shape: {balanced_test_df.shape}")
    print(f"Balanced Train Shape: {balanced_train_df.shape}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Balance labels in CSV files.')
    parser.add_argument('--dev', type=str, required=True, help='Path to the development CSV file')
    parser.add_argument('--test', type=str, required=True, help='Path to the test CSV file')
    parser.add_argument('--train', type=str, required=True, help='Path to the train CSV file')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output CSV files')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    
    args = parser.parse_args()
    main(args.dev, args.test, args.train, args.output_dir, args.seed)
