# PiCTURE
Pipeline for CRISPR-induced Transcriptome Unintended RNA Editing Analysis

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
git clone https://github.com/KazukiNakamae/PtTIDE;
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
