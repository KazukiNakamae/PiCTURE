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
joint_preparation;
```

それぞれの引数の説明はこちらになります。
```
<raw forward.fastq> & <raw reverse.fastq>: RAN-seqのraw fastqデータのパスを入力します。gzipファイルは受け付けていません。
<sample name>: サンプル名です。ユニークな文字列であれば何度もよいですが、基本的にはSRAのRun IDを入力することを推奨します。
<output directory name>: 出力ディレクトリ名です。preparation.shで入力したものと同一である必要があります。
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
variant_identification_from_multidb.sh　\
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
get_intersection_variants.sh　\
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
