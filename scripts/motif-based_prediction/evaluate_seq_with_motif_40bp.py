import pandas as pd
import argparse

# IUPAC nucleotide code dictionary for pattern matching
iupac_codes = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def matches_pattern(seq, pattern):
    return all(seq[i] in iupac_codes[pattern[i]] for i in range(len(pattern)))

def process_file(file_path, pattern, output_file_path):
    # Load the uploaded CSV file
    data = pd.read_csv(file_path)

    # Initialize columns for logits
    data['logits_0'] = 1
    data['logits_1'] = 0

    # Check each sequence for the pattern at positions 20, 21, 22
    for index, row in data.iterrows():
        sequence = row['sequence']
        if len(sequence) >= 22 and matches_pattern(sequence[19:22], pattern):
            data.at[index, 'logits_0'] = 0
            data.at[index, 'logits_1'] = 1

    # Add labels column
    data['labels'] = data['label']

    # Select necessary columns and save to a new CSV file
    output_data = data[['logits_0', 'logits_1', 'labels']]
    output_data.to_csv(output_file_path, index=False)

    # Count non-motif and motif regions
    non_motif_count = len(data[(data['logits_0'] == 1) & (data['labels'] == 1)])
    motif_count = len(data[(data['logits_1'] == 1) & (data['labels'] == 1)])

    # Print the counts
    print(f'Non-motif regions count: {non_motif_count}')
    print(f'Motif regions count: {motif_count}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a CSV file to evaluate DNA sequences.')
    parser.add_argument('file_path', type=str, help='Path to the input CSV file containing sequences.')
    parser.add_argument('pattern', type=str, help='DNA motif pattern to match.')
    parser.add_argument('output_file_path', type=str, help='Path to save the output CSV file.')

    args = parser.parse_args()

    process_file(args.file_path, args.pattern, args.output_file_path)
