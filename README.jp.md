# HMMER 3.4 PEZY-SC 対応版

HMMER 3.4 の一部の処理を PEZY-SC3 で実行するように改変した版です。
このファイルでは、このリポジトリのソースからビルドしてインストールする手順と、実行時の動作・設定を説明します。
本家 HMMER のインストール手順、オプションなどの説明は `official_readme/` 配下のファイルと `Userguide.pdf` にあります。

English version: [README.md](README.md)

以下では、次のプレースホルダを使います。

| 記号 | 意味 |
|---|---|
| `<SRC>` | このリポジトリを clone したディレクトリ（`configure` がある） |
| `<PREFIX>` | インストール先（`<PREFIX>/bin` に実行ファイルが入る） |
| `<PZSDK>` | PZSDK のディレクトリ（既定は `/opt/pezy/pzsdk`） |

---

## 1. 概要

- デバイスを使うのは `hmmsearch`・`phmmer`・`jackhmmer`・`hmmscan` の 4 つです。
  `hmmbuild` などそれ以外のプログラムは本家 HMMER 3.4 と同じく CPU で動きます。
  デバイスを使えない条件のときは、自動で CPU の処理に切り替わります（§5.2）。
- **PEZY-SC3 が無い環境でも、CPU だけで動く版（CPU 版）をビルドして使えます。**
  検索結果は本家 HMMER 3.4 と同じで、`hmmscan` はクエリ配列が多いときに本家より速くなります。
  ビルド方法は §3.5 を見てください。
- ビルドすると、ソースの commit ID が実行ファイルに埋め込まれ、`-h` の表示と結果ファイルのヘッダに出ます（§4）。
  どの版で出した結果かを後から確認できます。
  - commit ID は git から読むので、`git clone` したツリーでビルドしてください。
    追跡しているファイルに未コミットの変更があると、ID の末尾に `-dirty` が付きます。
  - GitHub の「Download ZIP」などで取得した `.git` の無いツリーでは、ID は `unknown-dirty` になります。
    その場合は、取得した commit ID を 1 行だけ書いたファイル `GIT_COMMIT` を `<SRC>` 直下に置くと、
    git より優先してその内容が使われます。

---

## 2. 必要なもの

| 項目 | 要件 | 備考 |
|---|---|---|
| OS / CPU | Linux x86_64 | 動作確認は AlmaLinux 9.1 |
| C コンパイラ | **システムの gcc** | 動作確認は gcc 11.3.1。**conda の gcc は使えない**（§7） |
| make / binutils | GNU make、`objcopy` | `objcopy` でデバイスカーネルを実行ファイルに埋め込む |
| PZSDK | `<PZSDK>/bin/pzcc`、`<PZSDK>/lib/libpzcl.so`、`pzcl/pzcl.h` | ヘッダは SDK の版によって `include/` か `inc/` にあるが、どちらでも自動で見つかる |
| デバイス | PEZY-SC3（ドライバ込み） | 見えているカードは全部使う（動作確認は SC3 × 4） |
| pthread | 必須 | |
| git | 推奨 | commit ID の埋め込みに使う（§1）。無くてもビルドはできる |
| python3 | 任意 | `make check` だけで使う |
| autoconf | **不要** | `configure` は生成済み |

PZSDK が見つからない環境でもビルドは通りますが、その場合は **PEZY を使わない CPU 版**ができます。
エラーは出ないので、ビルド後に必ず §4 の確認をしてください。

---

## 3. ビルドとインストール

PEZY-SC3 を使う版（PEZY 版）のビルド手順です。PEZY-SC3 を使わない CPU 版のビルドは §3.5 を見てください。

```bash
git clone https://github.com/t-oda-ic/hmmer_pz.git <SRC>
cd <SRC>

# (1) もし conda などを使用している場合、PATH をデフォルトのものに戻す（conda の gcc や ld を使わせない）
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

# (2) PZSDK の場所と、カーネルを生成するデバイスの世代
#     既定値は PZSDK_PATH=/opt/pezy/pzsdk、PZ_ARCH=sc3。既定のままでよければ不要
#     make と make install で値をそろえるため、コマンドライン引数ではなく環境変数で渡す
export PZSDK_PATH=<PZSDK>
export PZ_ARCH=sc3

# (3) configure → make → make install
#     --prefix を省くとインストール先は /usr/local（書き込みに root 権限が必要）
#     make だけで <SRC>/src/ に実行ファイルができ、そのまま使える（make install はそれを <PREFIX>/bin にコピーするだけ）
./configure --prefix=<PREFIX>
make -j32
make install
```

`make install` が終わると、次の実行ファイルができます。

```
<PREFIX>/bin/jackhmmer
<PREFIX>/bin/hmmsearch
<PREFIX>/bin/hmmscan
<PREFIX>/bin/hmmbuild
（ほかに phmmer, nhmmer, hmmpress, hmmalign などの HMMER 一式）
```

### 3.1 インストール先と、make install をしない場合

- **`--prefix` を指定しない場合、インストール先は `/usr/local` になります**（autoconf の既定値）。
  実行ファイルは `/usr/local/bin/`、man ページは `/usr/local/share/man/man1/` に入ります。
  `/usr/local` に書き込むには通常 root 権限が必要なので、一般ユーザーが `make install` すると権限エラーになります。
  自分のディレクトリに入れる場合は `--prefix` を指定してください。
- **`make install` をしなくても、`make` だけで `<SRC>/src/` に実行ファイルができ、そのまま実行できます**
  （`<SRC>/src/hmmsearch`、`<SRC>/src/hmmscan`、`<SRC>/src/jackhmmer`、`<SRC>/src/hmmbuild` など）。
  `make install` は、ここにできたものを `<PREFIX>/bin` にコピーするだけです。
  - HMMER と Easel のライブラリは静的にリンクされ、デバイスカーネルも実行ファイルに埋め込まれています。
    `libpzcl.so` は rpath で `<PZSDK>/lib` から読み込みます。このため、`src/` に置いたままでも動きます。
  - `make clean` や `make distclean` を実行すると、`src/` の実行ファイルも消えます。
  - この場合、§4 の `<PREFIX>/bin` は `<SRC>/src` に読み替えてください。
- 次のものは `make` でビルドされますが、`make install` ではインストールされません。
  必要ならビルドした場所から直接使うか、手でコピーしてください。
  - `src/` の補助プログラム `hmmc2`、`hmmerfm-exactmatch`
  - Easel のツール `esl-*`（`easel/miniapps/` にできる）

### 3.2 注意

- **`configure` を作り直さないでください**（`autoconf` や `autoreconf` を実行しない）。
  同梱の `configure` は、この版で変更した `configure.ac` から autoconf 2.71 で生成したものです。
  別の版の autoconf で作り直すと内容が変わります。
- **一度ビルドしたツリーを別の場所へコピーした場合や、`--prefix` を変える場合は、
  先に `make distclean` を実行してください。**
  `Makefile` と `config.status` には configure を実行したときのインストール先が書き込まれています。
  また、古い `.o` は作り直されません。
- **SDK の場所や版、`PZ_ARCH`、`HMMER_NO_PEZY` を変えたときは、`make clean` からやり直してください。**
  これらは make 変数なので、値が変わっても既存の `.o` は作り直されません。

### 3.3 configure のオプション

`--prefix` 以外は指定しなくてかまいません。既定の設定は本家 HMMER 3.4 を普通にビルドしたときと同じです
（`CC=gcc`、`CFLAGS=-O3`、DP 実装 `sse`、スレッド有効）。
本家 HMMER と結果を比べる場合は、`CFLAGS` などを変えないでください。

### 3.4 make のオプション

環境変数でもコマンドライン引数でも指定できます。

| 変数 | 既定 | 効果 |
|---|---|---|
| `PZSDK_PATH` | `/opt/pezy/pzsdk` | PZSDK のディレクトリ |
| `PZ_ARCH` | `sc3` | `pzcc --arch` に渡す世代。動作確認済みは `sc3` だけ |
| `HMMER_NO_PEZY=1` | 未設定 | PEZY を組み込まない（CPU 版をビルドする） |
| `OBJCOPY` | `objcopy` | `objcopy` の名前が違う環境で指定する |

PEZY を組み込むのは、`<PZSDK>/bin/pzcc`、`pzcl/pzcl.h`（`include/` か `inc/` のどちらか）、
`<PZSDK>/lib/libpzcl.so` の 3 つがすべて見つかったときだけです。このとき次のようにビルドします。

- `src/` の `.c` を `-DHAVE_PEZY` 付きでコンパイルする
- `src/p7_pz_kernel.pzc` を `pzcc --arch $(PZ_ARCH) -O3 -ffp-contract=off` でコンパイルし、
  `objcopy` で `libhmmer.a` に埋め込む。実行時にカーネルのファイルは要らない
  （`-ffp-contract=off` は CPU 版と結果をビット単位で一致させるために必要なので、外さないこと）
- `-lpzcl -lstdc++` と `-Wl,-rpath,<PZSDK>/lib` を付けてリンクする

### 3.5 CPU 版のビルド（PEZY-SC3 を使わない場合）

PEZY-SC3 が無い環境でも、CPU だけで動く版（CPU 版）をビルドして使えます。

- 必要なものは本家 HMMER 3.4 と同じです。§2 のうち PZSDK、デバイス、`objcopy` は要りません。
- 検索結果は本家 HMMER 3.4 と同じです。
- `hmmsearch`・`phmmer`・`jackhmmer` などは、本家と同じ処理で動きます。
- `hmmscan` は、クエリ配列が `HMMER_HMMCACHE_MINQ`（既定 4）本以上あると、プロファイル DB を一度だけメモリに
  読み込んで、すべてのクエリを処理します（本家はクエリごとにプロファイル DB を読み直します）。
  このため、クエリ配列が多いほど本家より速くなります。プロファイル DB が `HMMER_HMMCACHE_MAX`（既定 64 GiB）を
  超える場合は、本家と同じ処理になります。

**PZSDK が無い環境**では、PEZY 版と同じ手順でビルドすると自動的に CPU 版になります。
`PZSDK_PATH` と `PZ_ARCH` の設定は要りません。

```bash
git clone https://github.com/t-oda-ic/hmmer_pz.git <SRC>
cd <SRC>
# export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin   # 必要であるならば PATH をデフォルトのものに戻す（§3 (1) と同じ）
./configure --prefix=<PREFIX>
make -j32
make install
```

**PZSDK がある環境**で CPU 版をビルドするときは、`HMMER_NO_PEZY=1` を設定します。
`make` と `make install` の両方で必要なので、環境変数で渡してください。
`make install` のときに設定を忘れると、PEZY 用の部分がビルドされて PEZY 版と CPU 版が混ざった実行ファイルになります。
PEZY 版をビルドしたことのあるツリーでは、先に `make clean` を実行してください（§3.2）。

```bash
export HMMER_NO_PEZY=1
./configure --prefix=<PREFIX>
make -j32
make install
```

ビルドしたら、§4 の (b) で何も表示されず、(c) が 1 になる（バナーに `CPU-only binary, built without PZSDK` と出る）ことを確認してください。
`make check`（§6）は CPU 版でも実行できます。§6 の PEZY 単体テストと、`HMMER_PEZY` で始まる環境変数（§5.3）は CPU 版では使いません。

---

## 4. ビルド後の確認（必須）

PZSDK の検出に失敗すると、エラーを出さずに CPU 版ができます。CPU 版でも検索結果は正しいため、
結果を比較するテストでは気づけません。次の 3 点を必ず確認してください。
`make install` をしていない場合は、`<PREFIX>/bin` を `<SRC>/src` に読み替えてください。

```bash
# (a) 4 つの実行ファイルがあること
ls -l <PREFIX>/bin/{jackhmmer,hmmsearch,hmmscan,hmmbuild}

# (b) libpzcl.so がリンクされていること（何も出なければ CPU 版）
ldd <PREFIX>/bin/hmmsearch | grep pzcl

# (c) バナーに "CPU-only binary" が付いていないこと（0 なら PEZY 版）
<PREFIX>/bin/hmmsearch -h | head -6
<PREFIX>/bin/hmmsearch -h | grep -c 'CPU-only binary'
```

(c) は次のように表示されます。

```
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# PEZY-SC3 accelerated build (PEZY Computing K.K.)                                          <- PEZY 版
# PEZY-SC3 accelerated build (PEZY Computing K.K.); CPU-only binary, built without PZSDK    <- CPU 版
# git commit <commit ID>
# contact: https://github.com/t-oda-ic/hmmer_pz/issues
```

- `# git commit` の行には、ビルドしたソースの commit ID が出ます（§1）。結果ファイルのヘッダにも同じ行が入るので、
  どの版で出した結果かを後から確認できます。`# contact` の行は `-h` の表示にだけ出ます。
- `PEZY-SC3` は固定の文字列なので、`PZ_ARCH` を変えても表示は変わりません。

デバイスを使えるかどうかは、同梱のチュートリアル用ファイルを使った小さな検索で確認できます。

```bash
cd <作業ディレクトリ>       # pzcl.log がここに出る
HMMER_PEZY_VERBOSE=1 <PREFIX>/bin/hmmsearch --cpu 4 \
    <SRC>/tutorial/globins4.hmm <SRC>/tutorial/globins45.fa > /dev/null
#  -> stderr に "# pezy: 4 devices, ..." のように、開けたカード数が出る
#     "# pezy: <n> devices" の行が出なければ、デバイスを使わずに CPU で処理している
#     （"# pezy: [0.000s] start" のような時刻の行は、CPU で処理するときにも出る）
grep -E -i 'error|MMU' pzcl.log    # 何も出ないこと
```

---

## 5. 実行時の動作と設定

### 5.1 複数のカードの使い方

既定では、見えているカードをすべて使います。開けなかったカード（他のプロセスが使用中のカードなど）は
飛ばして、残りのカードで実行します。1 枚も開けなかったときだけ CPU で処理します。

- `hmmsearch` / `phmmer` / `jackhmmer`: 検索対象の配列 DB をカードの枚数で分け、各カードが自分の担当分を処理します。
- `hmmscan`: プロファイル DB を全カードに置き、空いたカードから順にモデルを割り当てます。

`PZCL` の実行ログ `pzcl.log` は、プログラムを実行したディレクトリに出ます。

### 5.2 CPU で処理する条件

次の場合は、デバイスを使わずに本家と同じ CPU の処理で実行します（結果は変わりません）。

| 対象 | 条件 |
|---|---|
| 全体 | デバイスが無い、全カードが他のプロセスで使用中、`HMMER_PEZY=0` を指定した |
| `hmmsearch` / `phmmer` / `jackhmmer` | 検索対象の配列 DB が次のどれかに当たる: 通常ファイルの FASTA でない（パイプ、gzip、Stockholm など）、ギャップ `-` を含む、100,000 残基を超える配列がある、ファイルが `HMMER_PEZY_MAXDB`（既定 8 GiB）を超える |
| `hmmsearch` / `phmmer` | `--restrictdb_stkey` か `--restrictdb_n` を指定した |
| 上記 3 つ（クエリごと） | クエリのモデル長が 2040 を超える（そのクエリだけ CPU で処理する） |
| `hmmscan` | クエリ配列が `HMMER_HMMCACHE_MINQ`（既定 4）本より少ない、プロファイル DB の大きさが `HMMER_HMMCACHE_MAX`（既定 64 GiB）を超える |

カーネルが制限時間内に終わらなかった場合も、CPU で処理し直します（stderr に警告が出ます）。
`hmmsearch` / `phmmer` / `jackhmmer` ではそのクエリだけを、`hmmscan` ではそれ以降のすべてのクエリを CPU で処理します。

### 5.3 環境変数

| 変数 | 既定 | 効果 |
|---|---|---|
| `HMMER_PEZY=0` | — | デバイスを使わず、本家と同じ CPU の処理で実行する（結果の比較用）。`0`・`n`・`N` で始まる値で無効になる |
| `HMMER_PEZY_NDEV=<n>` | 見えているカード全部 | 使うカードの数の上限。`1` で 1 枚だけ使う |
| `HMMER_PEZY_VERBOSE=1` | — | 開いたカードの数、カードごとの処理量の割合、処理ごとの所要時間を stderr に出す（診断用） |
| `HMMER_PEZY_MAXDB=<bytes>` | 8 GiB | `hmmsearch` などで、デバイス用にメモリへ読み込む配列 DB の大きさの上限 |
| `HMMER_PEZY_RTHREADS=<n>` | `--cpu` + 1 | 配列 DB の読み込みに使うスレッド数 |
| `HMMER_HMMCACHE_MINQ=<n>` | 4 | `hmmscan`: クエリ配列がこの本数以上あれば、プロファイル DB をメモリに読み込んでデバイスを使う |
| `HMMER_HMMCACHE_MAX=<bytes>` | 64 GiB | `hmmscan`: メモリに読み込むプロファイル DB の大きさの上限 |
| `HMMER_HMMCACHE_QBATCH=<n>` | 512 | `hmmscan`: デバイスに一度に渡すクエリ配列の数 |

通常は `HMMER_PEZY`・`HMMER_PEZY_NDEV`・`HMMER_PEZY_VERBOSE` 以外を変える必要はありません。

---

## 6. テスト

```bash
cd <SRC>
make check        # HMMER 本体のテストスイート（python3 が必要）
```

PEZY 部分の単体テストは次のようにビルドして実行します。どれも、デバイス（または高速化した読み込み）の結果が
本家の CPU 実装の結果と完全に一致するかを調べます。

```bash
cd <SRC>/src
make pz_ssv_utest pz_ssvmulti_utest pz_fwdparse_utest pz_fpprobe_utest pz_pzdb_utest

./pz_ssv_utest      [-n <nseq>] <hmmfile> <seqfile>                 # SSV フィルタ（hmmsearch 用）
./pz_ssvmulti_utest [-n <nmodel>] [-t <threads>] <hmmfile> <seqfile> # 複数モデルの SSV フィルタ（hmmscan 用。<hmmfile> は hmmpress 済み）
./pz_fwdparse_utest [-n <nseq>] <hmmfile> <seqfile>                 # parser 版 Forward / Backward（impl_sse とビット単位で比較）
./pz_fpprobe_utest  [-n <ntrial>]                                   # デバイスの浮動小数点演算が x86 SSE とビット単位で一致するか
./pz_pzdb_utest     [-n <nsample>] [-t <threads>] <seqfile>         # 配列 DB の高速読み込み（esl_sqio と比較）
```

カードの枚数を変えて試すと、DB の分け方の誤りも確認できます。

```bash
for nd in 1 2 3 4; do
  HMMER_PEZY_NDEV=$nd ./pz_ssv_utest <hmmfile> <seqfile> | tail -1
done
```

---

## 7. トラブルシュート

| 症状 | 原因と対処 |
|---|---|
| (c) の行に `CPU-only binary` が付く / `ldd` に `pzcl` が出ない | PZSDK が見つかっていない。`PZSDK_PATH` を確認し、`<PZSDK>/bin/pzcc`、`<PZSDK>/{include,inc}/pzcl/pzcl.h`、`<PZSDK>/lib/libpzcl.so` があるか見る。直したら **`make clean` からやり直す** |
| リンク時に `libstdc++.so.6` 関係のエラーや `undefined reference to std::...` が出る | conda の gcc / ld を使っている。`libpzcl.so` はシステムの libstdc++ を必要とする。§3 (1) のように PATH から conda を外す |
| 実行時に `libpzcl.so: cannot open shared object file` が出る | ビルドしたときと実行するときで SDK の場所が違う。実行ファイルは `<PZSDK>/lib` を rpath に持つ。同じパスに SDK を置くか、`LD_LIBRARY_PATH=<SDK の lib>` を設定する |
| `make install` が意図しない場所に書き込もうとする | 前にビルドしたときの `Makefile` が残っている。`make distclean` してから configure し直す |
| `# git commit` が `unknown-dirty` になる / 別の ID が出る | `.git` の無いツリー（ZIP で取得したものなど）でビルドした。§1 のとおり `GIT_COMMIT` を置くか、`git clone` したツリーでビルドする。ツリーを別の git リポジトリの中に置くと、そのリポジトリの ID が出る |
| `# git commit` の ID に `-dirty` が付く | ビルドしたツリーの、git が追跡しているファイルに未コミットの変更がある |
| デバイスを使わずに CPU で処理する | §5.2 の条件に当たっている。`HMMER_PEZY_VERBOSE=1` を付けても `# pezy: <n> devices` の行が出なければ CPU で処理している（理由は表示されない） |
| 一部のカードを使わない | 開けなかったカードは飛ばす。`HMMER_PEZY_VERBOSE=1` で `skipping a device` と理由が出る |
| 実行が止まる・異常終了する | 実行したディレクトリの `pzcl.log` を `grep "MMU ERROR"` などで確認する |

---

## 8. ライセンスと問い合わせ

- ライセンスは本家 HMMER と同じ BSD 3-Clause です。本家 HMMER の著作権表示を含め、`LICENSE` を参照してください。
  Easel のライセンスは `easel/LICENSE` にあります。
- 不具合の報告や問い合わせは [GitHub Issues](https://github.com/t-oda-ic/hmmer_pz/issues) にお願いします。
  この改変版についての問い合わせは、本家 HMMER の開発者ではなくこちらにお願いします。
