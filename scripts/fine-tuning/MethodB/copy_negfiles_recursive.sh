#!/bin/bash

# ファイルリスト
file_list=("SRR11561325" "SRR11561326" "SRR11561303" "SRR11561314" "SRR11561278" "SRR11561289")

# 指定したディレクトリ
source_directory="/Volumes/TBT3_SSD/report/*/all/sequence"

# negativeデータディレクトリ
destination_directory=$(pwd)"/negative"

# ファイルのコピー
for file_id in "${file_list[@]}"; do
    # findコマンドでファイルを検索し、見つかったファイルをカレントディレクトリにコピー
    find ${source_directory} -type f -name "${file_id}.hg38.identified.snp.fltr.CtoT.fa" -exec cp {} "$destination_directory" \;
    # コピーが成功したかをチェック
    if [ -e "${destination_directory}/${file_id}.hg38.identified.snp.fltr.CtoT.fa" ]; then
        echo "Copied: ${file_id}.hg38.identified.snp.fltr.CtoT.fa to $destination_directory"
    else
        echo "File not found: ${file_id}.hg38.identified.snp.fltr.CtoT.fa"
    fi
done

