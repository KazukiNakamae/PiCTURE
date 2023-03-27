# PiCTURE
Pipeline for CRISPR-induced Transcriptome Unintended RNA Editing Analysis

## Run PiCTURE with Docker

### Download PiCTURE pipeline using preparation.sh

最新のPiCTUREをダウンロードして、Dockerディレクトリをカレンとディレクトリに設定してください。

```
git clone https://github.com/KazukiNakamae/PiCTURE.git;
cd PiCTURE/Docker;
```

### Parepare PiCTURE pipeline using preparation.sh

PiCTUREパイプラインは出力ディレクトリ内で処理が完結するように設計されています。
まず次のコマンドを実行して出力ディレクトリとDockerイメージの作成を行なってください。

```
preparation.sh　<output directory name>
```

[GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)からhg38用のダウンロードページに入り、以下のファイルをダウンロードしてください。
- resources-broad-hg38-v0-Homo_sapiens_assembly38.dict 
- resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta
- resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai
- resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf
- resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx

ダウンロードしたファイルを下記のコマンドにしたがって出力ディレクトリ内に配置してください。
```
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.dict <output directory name>/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta <output directory name>/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai <output directory name>/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf <output directory name>/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx <output directory name>/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.dict <output directory name>/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai <output directory name>/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai;
cp resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta <output directory name>/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta;
```

完了です。この出力ファイルはサンプルごとに何度も作る必要はなく、一度作ってしまえばあとは下記のrun.shで使いまわせます。

### Run PiCTURE pipeline using run.sh

各サンプルごとに次のコマンドを実行してバリアントコールまでを実施してください。

```
run.sh　\
<raw forward.fastq> \
<raw reverse.fastq> \
<sample name> \
<output directory name> \
joint_preparation;
```

それぞれの引数の説明はこちらになります。
```
<raw forward.fastq> & <raw reverse.fastq>: RAN-seqのraw fastqデータのパスを入力します。gzipファイルは受け付けていません。
<sample name>: サンプル名です。ユニークな文字列であれば何度もよいですが、基本的にはSRAのRun IDを入力することを推奨します。
<output directory name>: 出力ディレクトリ名です。preparation.shで入力したものと同一である必要があります。
```

完了すると、出力ディレクトリ内に各ステップの処理結果が出力されます。

## 開発工程におけるgithubの扱い方

まず最初はリポジストリをクローンしてもらいます。それ以降は鈴木さんはmainからtestブランチを作成して、testブランチへプッシュを行うようにしましょう。
私がそれを確認して、追加の修正を指示やコーディングを行ったあと、マージしてmainのほうへ反映します。
具体的な流れとしては、

- (中前) issue欄に課題等を記述する
- (鈴木さん) リポジストリをクローン（操作A）
- (鈴木さん) mainからtestブランチを作成（操作B）
- (鈴木さん) (*)コーディング
- (鈴木さん) testブランチへのプッシュ（操作C）
- (鈴木さん) サイトへ行ってpull requestをだす
- (中前) 中前が許可して、ブランチの内容をmainへマージする
- (中前) 該当issueをクローズする
- (中前) コーディング・修正作業
- (中前) issue欄に次の課題等を記述する
- (鈴木さん) mainからのプル（操作D）
- (鈴木さん) mainからtestブランチを作成（操作B）
- (*)に戻る

とする予定です。

---

### （操作A）リポジストリをクローン

```bash
git clone https://github.com/KazukiNakamae/PiCTURE;
```

### （操作B）mainからtestブランチを作成

```bash
# ブランチの作成と切り替え
# ブランチ名はtest_20210825_nakamae、みたいな感じでお願いします。
git checkout -b test_(日付)_(名前);
# ファイルのダウンロード
git pull origin main;
```

### （操作C）testブランチへのプッシュ

```bash
# 更新ファイルの登録
git add --all;
# コミット
git commit -m "add new file";
# リモートリポジトリの設定
git remote add origin https://github.com//（自分のユーザー名）.git;
# ブランチへプッシュ
git push origin test_(日付)_(名前);
```

### （操作D）mainからのプル

```bash
# ブランチの確認
git branch;
# mainへ移行
git checkout main;
# mainの内容を反映
git pull origin main;
```

### mainへのプッシュ（中前のみ使用）

```bash
# 更新ファイルの登録
git add --all;
# コミット
git commit -m "add new file";
# リモートリポジトリの設定
git remote add origin https://github.com//KazukiNakamae.git;
# mainへプッシュ
git push origin main;
```
