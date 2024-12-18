import os
import csv
import argparse
from sklearn.model_selection import train_test_split

def read_fasta(file_path):
    sequences = []
    with open(file_path, "r") as file:
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return sequences

def remove_duplicates(seq_pos, seq_neg):
    print("Remove sequences exist in the both labels from positive dataset")
    return [item for item in seq_pos if item not in seq_neg]

def split_data(sequences, train_size=0.7, random_state=0):
    train, temp = train_test_split(sequences, train_size=train_size, random_state=random_state)
    develop, test = train_test_split(temp, train_size=0.5, random_state=random_state)
    return train, develop, test

def write_data(pos_list, neg_list, out_fn):
    print("Add labels; positive=1, negative=0")
    data = [{"sequence": seq, "label": 1} for seq in pos_list]
    data += [{"sequence": seq, "label": 0} for seq in neg_list]
    
    field_names = ["sequence", "label"]
    with open(out_fn, "w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter=',')
        writer.writeheader()
        writer.writerows(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some FASTA files.")
    parser.add_argument("--positive", type=str, required=True, help="Path to the positive FASTA file")
    parser.add_argument("--negative", type=str, required=True, help="Path to the negative FASTA file")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for data splitting")
    parser.add_argument("--outdir", type=str, default=".", help="output directory")
    args = parser.parse_args()

    # メイン処理
    seq_pos = read_fasta(args.positive)
    seq_neg = read_fasta(args.negative)

    print(f"Sequences with positive label: {len(seq_pos)}")
    print(f"Sequences with negative label: {len(seq_neg)}")

    # Negativeデータセットにも存在している配列をPositiveデータセットから除外
    seq_pos = remove_duplicates(seq_pos, seq_neg)
    print(f"Filtered sequences with positive label: {len(seq_pos)}")

    # シャッフルしつつデータを分割（train:dev:test=7:1.5:1.5）
    seq_pos_train, seq_pos_develop, seq_pos_test = split_data(seq_pos, random_state=args.seed)
    seq_neg_train, seq_neg_develop, seq_neg_test = split_data(seq_neg, random_state=args.seed)

    # ラベル（positive=1/negative=0）を付加した形式でファイルを出力
    write_data(seq_pos_train, seq_neg_train, os.path.join(args.outdir, "train.csv"))
    write_data(seq_pos_develop, seq_neg_develop, os.path.join(args.outdir, "dev.csv"))
    write_data(seq_pos_test, seq_neg_test, os.path.join(args.outdir, "test.csv"))
