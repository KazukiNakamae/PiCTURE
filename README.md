# PiCTURE
**Pi**peline for **C**RISPR-induced **T**ranscriptome-wide **U**nintended **R**NA **E**diting

## Directory Structure

```
- PiCTURE
  -- Docker              # A collection of scripts to execute each step of PiCTURE on Docker
  -- Singularity         # A collection of scripts to execute each step of PiCTURE on Singularity (Apptainer)
  -- scripts             # Auxiliary scripts for processing PiCTURE outputs
    --- visualization    # Auxiliary scripts for additional visualization
      ---- chr_dist_from_vcf  # Used for visualizing Chromosome-scale distribution from VCF files
    --- fine-tuning      # Scripts for creating datasets for fine-tuning
      ---- MethodA       # Scripts used to create STL models
      ---- MethodB       # Scripts used to create STL models
    --- motif-based_prediction # Scripts for counting specific motifs in sequence data and performing motif-based predictions
```

## Run PiCTURE with Docker

### Download the PiCTURE Pipeline using `preparation.sh`

Download the latest version of PiCTURE and set the Docker directory as the current directory.

```bash
git clone https://github.com/KazukiNakamae/PiCTURE.git
cd PiCTURE/Docker
```

### Prepare the PiCTURE Pipeline using `preparation.sh`

The PiCTURE pipeline is designed to operate entirely within the output directory.

After setting the Docker directory as the current directory, execute the following command to create the output directory and Docker image:

```bash
chmod +x preparation.sh
sudo bash preparation.sh <output directory name>
```

That's it! This output file only needs to be created once and can be reused for multiple samples using the `run.sh` script below.

### Run the PiCTURE Pipeline

For each sample, execute the following command to run the pipeline up to BQSR:

```bash
chmod +x run.sh
sudo bash run.sh \
  <raw forward.fastq> \
  <raw reverse.fastq> \
  <sample name> \
  <output directory name> \
  <memory>
```

The arguments are described as follows:

- `<raw forward.fastq>` & `<raw reverse.fastq>`: Paths to the raw FASTQ data of RNA-seq. Gzip files are not supported.
- `<sample name>`: The sample name. Any unique string is acceptable, but it is recommended to use the SRA Run ID.
- `<output directory name>`: The name of the output directory. This must match the one used in `preparation.sh`.
- `<memory>`: Specifies the amount of memory to use.

#### Pipeline A: Processing a Single RNA-seq Dataset

Perform variant calling and generate a VCF file from a single RNA-seq dataset.

```bash
chmod +x variant_identification_from_singledb.sh
sudo bash variant_identification_from_singledb.sh \
  <sample name> \
  <output directory name>
```

The arguments are described as follows:

- `<sample name>`: The sample name. Use the same name as in `run.sh`.
- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.

Once complete, execute the following command for SNP detection and motif extraction:

```bash
chmod +x motif_estimation.sh
sudo bash motif_estimation.sh \
  <sample name> \
  <output directory name> \
  <VAF threshold>
```

The arguments are described as follows:

- `<sample name>`: The sample name. Use the same name as in `run.sh`.
- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<VAF threshold>`: The threshold for classifying VAF. Input a value between 0.0 and 1.0.

To test various thresholds, simply modify the threshold value above. It is not necessary to create a separate output directory for each threshold.

Afterward, package the result files into a `tar.gz` format using the following command:

```bash
chmod +x get_result.sh
sudo bash get_result.sh \
  <sample name> \
  <output directory name> \
  <result_name>
```

The arguments are described as follows:

- `<sample name>`: The sample name. Use the same name as in `run.sh`.
- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<result_name>`: The name of the packaged tar.gz data. A file named `<result_name>.tar.gz` will be created under `<output directory name>/report`.

The tar.gz file includes all data for each threshold tested for the sample.

#### Pipeline B: Processing Multiple RNA-seq Datasets

Perform variant calling and generate VCF files from multiple RNA-seq datasets.

```bash
chmod +x variant_identification_from_multidb.sh
sudo bash variant_identification_from_multidb.sh \
  <output directory name> \
  <group name> \
  <sample name 1> \
  <sample name 2> \
  <sample name 3> \
  ... \
  <sample name n>
```

The arguments are described as follows:

- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<group name>`: The analysis group name. Specify any unique name.
- `<sample name n>`: The sample names to be analyzed. Use the same names as in `run.sh`. There is no limit on the number of samples.

Next, execute the following command for SNP detection and motif extraction:

```bash
chmod +x motif_estimation.sh
sudo bash motif_estimation.sh \
  <group name> \
  <output directory name> \
  <VAF threshold>
```

The arguments are described as follows:

- `<group name>`: The analysis group name. Use the same name as in `variant_identification_from_multidb.sh`.
- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<VAF threshold>`: The threshold for classifying VAF. Input a value between 0.0 and 1.0.

Finally, package the results into a `tar.gz` format using:

```bash
chmod +x get_result.sh
sudo bash get_result.sh \
  <group name> \
  <output directory name> \
  <result_name>
```

The arguments are described as follows:

- `<group name>`: The analysis group name. Use the same name as in `variant_identification_from_multidb.sh`.
- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<result_name>`: The name of the packaged tar.gz data. A file named `<result_name>.tar.gz` will be created under `<output directory name>/report`.

The tar.gz file includes all data for each threshold tested for the group.

#### Pipeline C: Extracting Shared Variants Across Datasets

Requires outputs from Pipeline A or B for at least two datasets. Extract shared variants across multiple VCF files.

```bash
chmod +x get_intersection_variants.sh
sudo bash get_intersection_variants.sh \
  <output directory name> \
  <set name> \
  <sample name or group name 1> \
  <sample name or group name 2> \
  ... \
  <sample name or group name n>
```

The arguments are described as follows:

- `<output directory name>`: The output directory name. Use the same name as in `run.sh`.
- `<set name>`: The name of the shared variant set. Specify a unique name different from `<sample name>` or `<group name>`.
- `<sample name or group name n>`: The sample or group names to be analyzed. There is no limit on the number of items.

Then perform SNP detection and motif extraction:

```bash
chmod +x motif_estimation.sh
sudo bash motif_estimation.sh \
  <set name> \
  <output directory name> \
  <VAF threshold>
```

Finally, package the results into a `tar.gz` format:

```bash
chmod +x get_result.sh
sudo bash get_result.sh \
  <set name> \
  <output directory name> \
  <result_name>
```

### Auxiliary Scripts

Refer to the respective `README.md` files located under the `scripts` directory for details.

<details>
<summary>For Developers</summary>

## ディレクトリ構成

```
- PiCTURE
-- Docker # Docker上でPiCTUREの各ステップを実行するためのスクリプト集
-- Singularity # Singularity（Apptainer）上でPiCTUREの各ステップを実行するためのスクリプト集
-- scripts # PiCTUREの出力を加工するための補助スクリプト
--- visualization # 追加のvisualizationを行うための補助スクリプト
---- chr_dist_from_vcf # VCFファイルを元にChromosome-scale ditributionを描写するために利用
--- fine-tuning # ファインチューニングに利用するためのデータセット作成スクリプト
---- MethodA # STL modelの作成に利用したスクリプト
---- MethodB # STL modelの作成に利用したスクリプト
--- motif-based_prediction # 配列データ内の特定モチーフのカウントやモチーフベースの予測を行うためのスクリプト
```

## Run PiCTURE with Docker

### Download PiCTURE pipeline using preparation.sh

最新のPiCTUREをダウンロードして、Dockerディレクトリをカレンとディレクトリに設定してください。

```
git clone https://github.com/KazukiNakamae/PiCTURE.git;
cd PiCTURE/Docker;
```

### Parepare PiCTURE pipeline using preparation.sh

PiCTUREパイプラインは出力ディレクトリ内で処理が完結するように設計されています。

Dockerディクレクトリをカレントディレクトリに設定後、次のコマンドを実行して出力ディレクトリとDockerイメージの作成を行なってください。

```
chmod +x preparation.sh;
sudo bash preparation.sh　<output directory name>;
```

完了です。この出力ファイルはサンプルごとに何度も作る必要はなく、一度作ってしまえばあとは下記のrun.shで使いまわせます。

### Run PiCTURE pipeline

各サンプルごとに次のコマンドを実行してBQSRまでを実施してください。

```
chmod +x run.sh;
sudo bash run.sh　\
<raw forward.fastq> \
<raw reverse.fastq> \
<sample name> \
<output directory name> \
<memory>;
```

それぞれの引数の説明はこちらになります。
```
<raw forward.fastq> & <raw reverse.fastq>: RAN-seqのraw fastqデータのパスを入力します。gzipファイルは受け付けていません。
<sample name>: サンプル名です。ユニークな文字列であれば何度もよいですが、基本的にはSRAのRun IDを入力することを推奨します。
<output directory name>: 出力ディレクトリ名です。preparation.shで入力したものと同一である必要があります。
<memory>： 使用するメモリ指定です。
```

#### PipelineA: 単一のRNA-seqデータを取得する場合

バリアントコールして、単一のRNA-seqデータからvcfファイルを生成します。

```
chmod +x variant_identification_from_singledb.sh;
sudo bash variant_identification_from_singledb.sh　\
<sample name> \
<output directory name>;
```

それぞれの引数の説明はこちらになります。
```
<sample name>: サンプル名です。run.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
```

完了したら次のコマンドを実行して、SNP検出とモチーフ抽出を行います。

```
chmod +x motif_estimation.sh;
sudo bash motif_estimation.sh　\
<sample name> \
<output directory name> \
<VAF threshold>;
```

それぞれの引数の説明はこちらになります。
```
<sample name>: サンプル名です。run.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<VAF threshold>: VAFを分類するときの閾値です。0.0-1.0までの数値を入力してください。
```

様々な閾値を試したい場合は、上記の閾値を変更して処理してください。閾値ごとに別々に出力ディレクトリを作成する必要はありません。

完了したら次のコマンドを実行して、結果ファイルをtar.gz形式にパッケージングします。

```
chmod +x get_result.sh;
sudo bash get_result.sh　\
<sample name> \
<output directory name> \
<result_name>;
```

それぞれの引数の説明はこちらになります。
```
<sample name>: サンプル名です。run.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<result_name>: パッケージされたtar.gzデータの名前です。`<output directory name> /report`上に`<result_name>.tar.gz`というファイルが作成されます。
```

tar.gzファイルにはサンプル名を含む各閾値ごとのデータを全て含みます。

#### PipelineB: 複数のRNA-seqデータを取得する場合

複数のバリアントコールして、それらRNA-seqデータからvcfファイルを生成します。

```
chmod +x variant_identification_from_multidb.sh;
sudo bash variant_identification_from_multidb.sh　\
<output directory name> \
<group name>
<sample name 1> \
<sample name 2> \
<sample name 3> \
...
<sample name n>;
```

それぞれの引数の説明はこちらになります。
```
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<group name>: 解析グループ名です。重複のない任意の名前を指定してください。
<sample name n>: 解析対象とするサンプル名です。run.shと同じものを入力してください。対象とできるサンプル数は無制限となっています。
```

完了したら次のコマンドを実行して、SNP検出とモチーフ抽出を行います。

```
chmod +x motif_estimation.sh;
sudo bash motif_estimation.sh　\
<group name> \
<output directory name> \
<VAF threshold>;
```

それぞれの引数の説明はこちらになります。
```
<group name>: 解析グループ名です。variant_identification_from_multidb.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<VAF threshold>: VAFを分類するときの閾値です。0.0-1.0までの数値を入力してください。
```

様々な閾値を試したい場合は、上記の閾値を変更して処理してください。閾値ごとに別々に出力ディレクトリを作成する必要はありません。

完了したら次のコマンドを実行して、結果ファイルをtar.gz形式にパッケージングします。

```
chmod +x get_result.sh;
sudo bash get_result.sh　\
<group name> \
<output directory name> \
<result_name>;
```

それぞれの引数の説明はこちらになります。
```
<group name>: 解析グループ名です。variant_identification_from_multidb.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<result_name>: パッケージされたtar.gzデータの名前です。`<output directory name> /report`上に`<result_name>.tar.gz`というファイルが作成されます。
```

tar.gzファイルにはサンプル名を含む各閾値ごとのデータを全て含みます。

#### PipelineC: 異なるデータの共通のバリアントデータから取得する場合

こちらはPipelineAまたはPipelineBで2つ以上のデータを事前に出力済みである必要があります。
複数のvcfファイルから共通の変異を抽出します。

```
chmod +x get_intersection_variants.sh;
sudo bash get_intersection_variants.sh　\
<output directory name> \
<set name>
<sample name or group name 1> \
<sample name or group name 2> \
<sample name or group name 3> \
...
<sample name or group name n>;
```

それぞれの引数の説明はこちらになります。
```
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<set name>: 共通変異セットの名前です。<sample name><group name>とは異なる重複のない任意の名前を指定してください。
<sample name or group name n>: 解析対象とするサンプル名あるいはグループ名です。対象とできる数は無制限となっています。
```

完了したら次のコマンドを実行して、SNP検出とモチーフ抽出を行います。

```
chmod +x motif_estimation.sh;
sudo bash motif_estimation.sh　\
<set name> \
<output directory name> \
<VAF threshold>;
```

それぞれの引数の説明はこちらになります。
```
<set name>: 共通変異セットの名前です。get_intersection_variants.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<VAF threshold>: VAFを分類するときの閾値です。0.0-1.0までの数値を入力してください。
```

様々な閾値を試したい場合は、上記の閾値を変更して処理してください。閾値ごとに別々に出力ディレクトリを作成する必要はありません。

完了したら次のコマンドを実行して、結果ファイルをtar.gz形式にパッケージングします。

```
chmod +x get_result.sh;
sudo bash get_result.sh　\
<set name> \
<output directory name> \
<result_name>;
```

それぞれの引数の説明はこちらになります。
```
<set name>: 共通変異セットの名前です。get_intersection_variants.shと同じものを入力してください。
<output directory name>: 出力ディレクトリ名です。run.shと同じものを入力してください。
<result_name>: パッケージされたtar.gzデータの名前です。`<output directory name> /report`上に`<result_name>.tar.gz`というファイルが作成されます。
```

tar.gzファイルにはサンプル名を含む各閾値ごとのデータを全て含みます。

### 補助スクリプト

scriptsディレクトリの下位に配置されたそれぞれのREADME.mdをご参照ください。

</details>