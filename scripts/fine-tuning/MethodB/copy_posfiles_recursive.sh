#!/bin/bash

# ファイルリスト
file_list=("SRR11561273" "SRR11561299" "SRR11561297" "SRR11561298" "SRR11561323" "SRR11561324" "SRR8096262")

# 指定したディレクトリ
source_directory="/Volumes/TBT3_SSD/report/*/all/sequence"

# negativeデータディレクトリ
destination_directory=$(pwd)"/positive"

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

