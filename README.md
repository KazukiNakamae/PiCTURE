# PiCTURE
**Pi**peline for **C**RISPR-induced **T**ranscriptome-wide **U**nintended **R**NA **E**diting

## Run PiCTURE with Docker

### Download PiCTURE pipeline using preparation.sh

最新のPiCTUREをダウンロードして、Dockerディレクトリをカレンとディレクトリに設定してください。

```
git clone https://github.com/KazukiNakamae/PiCTURE.git;
cd PiCTURE/Docker;
```

### Parepare PiCTURE pipeline using preparation.sh

PiCTUREパイプラインは出力ディレクトリ内で処理が完結するように設計されています。
まず[GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)からhg38用のダウンロードページに入り、以下のファイルをダウンロードしてください。
- resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta

このファイルをDockerディクレクトリに配置後、
次のコマンドを実行して出力ディレクトリとDockerイメージの作成を行なってください。

```
chmod +x preparation.sh;
sudo preparation.sh　<output directory name>;
```

完了です。この出力ファイルはサンプルごとに何度も作る必要はなく、一度作ってしまえばあとは下記のrun.shで使いまわせます。

### Run PiCTURE pipeline

各サンプルごとに次のコマンドを実行してBQSRまでを実施してください。

```
chmod +x run.sh;
sudo run.sh　\
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

#### 単一のRNA-seqデータを取得する場合

バリアントコールして、単一のRNA-seqデータからvcfファイルを生成します。

```
chmod +x variant_identification_from_singledb.sh;
sudo variant_identification_from_singledb.sh　\
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
sudo motif_estimation.sh　\
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
chmod +x get_result_from_singledb.sh;
sudo get_result_from_singledb.sh　\
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
