# read a fa file and process each line

import csv
import os

file_pos = open("./dataset/positive/SRR11561298.hg38.identified.snp.fltr.CtoT.fa", "r")
file_neg = open("./dataset/negative/SRR11561314.hg38.identified.snp.fltr.CtoT.fa", "r")

seq_pos = []
for line in file_pos:
    if line[0] == ">":
        continue
    else:
        line2 = line.split("\n")[0]
        seq_pos.append(line2)


seq_neg = []
for line in file_neg:
    if line[0] == ">":
        continue
    else:
        line2 = line.split("\n")[0]
        # line3 = line2[30:-30]
        seq_neg.append(line2)

seq_pos = [item for item in seq_pos if item not in seq_neg]


data = []
for seq in seq_pos:
    data.append({"sequence":seq,"label":1})

for seq in seq_neg:
    data.append({"sequence":seq,"label":0})

import random
random.shuffle(data)

field_names = ["sequence","label"]
with open("/Users/suzuki/Desktop/train.csv","w") as output_file:
    writer = csv.DictWriter(output_file,fieldnames=field_names, delimiter=',')
    writer.writeheader()
    writer.writerows(data)