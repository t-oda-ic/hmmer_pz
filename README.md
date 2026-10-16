# HMMER 3.4 for PEZY-SC

A modified version of HMMER 3.4 that runs part of its work on PEZY-SC3.
This file describes how to build and install it from the source in this repository, and how it behaves at run time.
For the original HMMER installation instructions and option reference, see the files under `official_readme/` and `Userguide.pdf`.

日本語版: [README.jp.md](README.jp.md)

The following placeholders are used below.

| Placeholder | Meaning |
|---|---|
| `<SRC>` | The directory you cloned this repository into (it contains `configure`) |
| `<PREFIX>` | The install location (executables go to `<PREFIX>/bin`) |
| `<PZSDK>` | The PZSDK directory (default `/opt/pezy/pzsdk`) |

---

## 1. Overview

- Four programs use the device: `hmmsearch`, `phmmer`, `jackhmmer` and `hmmscan`.
  All other programs, such as `hmmbuild`, run on the CPU exactly as in the original HMMER 3.4.
  When the device cannot be used, the programs fall back to the CPU automatically (§5.2).
- **Even without PEZY-SC3, you can build and use a version that runs on the CPU only (the CPU-only build).**
  Its search results are the same as the original HMMER 3.4, and `hmmscan` is faster than the original when there are many query sequences.
  See §3.5 for how to build it.
- The build embeds the commit ID of the source in the executables; it appears in the `-h` output and in the header of result files (§4).
  This lets you tell later which version produced a given result.
  - The commit ID is read from git, so build from a tree obtained with `git clone`.
    If tracked files have uncommitted changes, `-dirty` is appended to the ID.
  - In a tree without `.git`, such as one obtained with GitHub's "Download ZIP", the ID becomes `unknown-dirty`.
    In that case, put a file named `GIT_COMMIT` containing just the commit ID on one line directly under `<SRC>`;
    its contents take precedence over git.

---

## 2. Requirements

| Item | Requirement | Notes |
|---|---|---|
| OS / CPU | Linux x86_64 | Tested on AlmaLinux 9.1 |
| C compiler | **The system gcc** | Tested with gcc 11.3.1. **conda's gcc does not work** (§7) |
| make / binutils | GNU make, `objcopy` | `objcopy` embeds the device kernel in the executables |
| PZSDK | `<PZSDK>/bin/pzcc`, `<PZSDK>/lib/libpzcl.so`, `pzcl/pzcl.h` | Depending on the SDK version the headers are under `include/` or `inc/`; either is found automatically |
| Device | PEZY-SC3 (with driver) | Every visible card is used (tested with 4 x SC3) |
| pthread | Required | |
| git | Recommended | Used to embed the commit ID (§1). The build works without it |
| python3 | Optional | Only used by `make check` |
| autoconf | **Not needed** | `configure` is pre-generated |

The build also succeeds where PZSDK cannot be found, but then it produces a **CPU-only build that does not use PEZY**.
No error is reported, so always run the checks in §4 after building.

---

## 3. Build and install

These are the steps for the build that uses PEZY-SC3 (the PEZY build). For the CPU-only build, which does not use PEZY-SC3, see §3.5.

```bash
git clone https://github.com/t-oda-ic/hmmer_pz.git <SRC>
cd <SRC>

# (1) If you are using conda, remove conda from PATH (so that conda's gcc and ld are not used)
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

# (2) PZSDK location and the device generation to build the kernel for
#     Defaults are PZSDK_PATH=/opt/pezy/pzsdk and PZ_ARCH=sc3; not needed if the defaults are fine
#     Pass them as environment variables, not command-line arguments, so make and make install see the same values
export PZSDK_PATH=<PZSDK>
export PZ_ARCH=sc3

# (3) configure -> make -> make install
#     Without --prefix the install location is /usr/local (writing there needs root)
#     make alone produces usable executables in <SRC>/src/ (make install just copies them to <PREFIX>/bin)
./configure --prefix=<PREFIX>
make -j32
make install
```

After `make install`, you have the following executables.

```
<PREFIX>/bin/jackhmmer
<PREFIX>/bin/hmmsearch
<PREFIX>/bin/hmmscan
<PREFIX>/bin/hmmbuild
(plus the rest of the HMMER suite: phmmer, nhmmer, hmmpress, hmmalign, ...)
```

### 3.1 Install location, and skipping make install

- **Without `--prefix`, the install location is `/usr/local`** (the autoconf default).
  Executables go to `/usr/local/bin/` and man pages to `/usr/local/share/man/man1/`.
  Writing to `/usr/local` usually needs root, so `make install` as an ordinary user fails with a permission error.
  Use `--prefix` to install into your own directory.
- **You do not need `make install`: `make` alone builds executables in `<SRC>/src/` that can be run in place**
  (`<SRC>/src/hmmsearch`, `<SRC>/src/hmmscan`, `<SRC>/src/jackhmmer`, `<SRC>/src/hmmbuild`, ...).
  `make install` only copies them to `<PREFIX>/bin`.
  - The HMMER and Easel libraries are linked statically and the device kernel is embedded in the executables.
    `libpzcl.so` is loaded from `<PZSDK>/lib` via rpath, so the executables also work where they are in `src/`.
  - `make clean` and `make distclean` delete the executables in `src/` as well.
  - In this case, read `<PREFIX>/bin` in §4 as `<SRC>/src`.
- The following are built by `make` but not installed by `make install`.
  Use them where they were built, or copy them by hand if needed.
  - The auxiliary programs `hmmc2` and `hmmerfm-exactmatch` in `src/`
  - The Easel tools `esl-*` (built in `easel/miniapps/`)

### 3.2 Cautions

- **Do not regenerate `configure`** (do not run `autoconf` or `autoreconf`).
  The bundled `configure` was generated with autoconf 2.71 from this version's modified `configure.ac`.
  Regenerating it with a different autoconf version changes its contents.
- **If you copy a built tree to another location, or change `--prefix`,
  run `make distclean` first.**
  `Makefile` and `config.status` record the install location given to configure,
  and stale `.o` files are not rebuilt.
- **After changing the SDK location or version, `PZ_ARCH`, or `HMMER_NO_PEZY`, start again from `make clean`.**
  These are make variables, so changing them does not cause existing `.o` files to be rebuilt.

### 3.3 configure options

Nothing other than `--prefix` needs to be specified. The defaults are the same as a plain build of the original HMMER 3.4
(`CC=gcc`, `CFLAGS=-O3`, DP implementation `sse`, threads enabled).
If you compare results with the original HMMER, do not change `CFLAGS` and the like.

### 3.4 make options

These can be given as environment variables or as command-line arguments.

| Variable | Default | Effect |
|---|---|---|
| `PZSDK_PATH` | `/opt/pezy/pzsdk` | PZSDK directory |
| `PZ_ARCH` | `sc3` | Device generation passed to `pzcc --arch`. Only `sc3` has been tested |
| `HMMER_NO_PEZY=1` | unset | Build without PEZY (CPU-only build) |
| `OBJCOPY` | `objcopy` | Set this where `objcopy` has a different name |

PEZY support is built in only when all three of `<PZSDK>/bin/pzcc`, `pzcl/pzcl.h` (under either `include/` or `inc/`),
and `<PZSDK>/lib/libpzcl.so` are found. In that case the build:

- compiles the `.c` files in `src/` with `-DHAVE_PEZY`
- compiles `src/p7_pz_kernel.pzc` with `pzcc --arch $(PZ_ARCH) -O3 -ffp-contract=off` and
  embeds it in `libhmmer.a` with `objcopy`, so no kernel file is needed at run time
  (`-ffp-contract=off` is required for bit-exact agreement with the CPU build; do not remove it)
- links with `-lpzcl -lstdc++` and `-Wl,-rpath,<PZSDK>/lib`

### 3.5 CPU-only build (without PEZY-SC3)

Even without PEZY-SC3, you can build and use a version that runs on the CPU only (the CPU-only build).

- The requirements are the same as for the original HMMER 3.4. Of the items in §2, PZSDK, the device and `objcopy` are not needed.
- Search results are the same as the original HMMER 3.4.
- `hmmsearch`, `phmmer`, `jackhmmer` and the others run the same code as the original.
- When there are at least `HMMER_HMMCACHE_MINQ` (default 4) query sequences, `hmmscan` reads the profile database into memory
  once and processes every query from there (the original re-reads the profile database for each query).
  So the more query sequences there are, the faster it is than the original. If the profile database is larger than
  `HMMER_HMMCACHE_MAX` (default 64 GiB), it runs the same code as the original.

**Where PZSDK is not installed**, following the same steps as for the PEZY build automatically produces the CPU-only build.
There is no need to set `PZSDK_PATH` or `PZ_ARCH`.

```bash
git clone https://github.com/t-oda-ic/hmmer_pz.git <SRC>
cd <SRC>
# export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin   # remove conda from PATH if necessary (same as §3 (1))
./configure --prefix=<PREFIX>
make -j32
make install
```

**Where PZSDK is installed**, set `HMMER_NO_PEZY=1` to build the CPU-only version.
It is needed by both `make` and `make install`, so pass it as an environment variable.
If you forget it for `make install`, the PEZY parts get built and the executables end up as a mix of the PEZY and CPU-only builds.
In a tree where the PEZY build was built before, run `make clean` first (§3.2).

```bash
export HMMER_NO_PEZY=1
./configure --prefix=<PREFIX>
make -j32
make install
```

After building, check that (b) in §4 prints nothing and (c) prints 1 (the banner says `CPU-only binary, built without PZSDK`).
`make check` (§6) also works with the CPU-only build. The PEZY unit tests in §6 and the environment variables starting with `HMMER_PEZY` (§5.3) do not apply to the CPU-only build.

---

## 4. Post-build checks (required)

If PZSDK detection fails, a CPU-only build is produced without any error. Because the CPU-only build also produces correct results,
tests that compare results will not catch it. Always check the following three points.
If you did not run `make install`, read `<PREFIX>/bin` as `<SRC>/src`.

```bash
# (a) the four executables exist
ls -l <PREFIX>/bin/{jackhmmer,hmmsearch,hmmscan,hmmbuild}

# (b) libpzcl.so is linked (no output means a CPU-only build)
ldd <PREFIX>/bin/hmmsearch | grep pzcl

# (c) the banner does not say "CPU-only binary" (0 means a PEZY build)
<PREFIX>/bin/hmmsearch -h | head -6
<PREFIX>/bin/hmmsearch -h | grep -c 'CPU-only binary'
```

The output of (c) looks like this.

```
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# PEZY-SC3 accelerated build (PEZY Computing K.K.)                                          <- PEZY build
# PEZY-SC3 accelerated build (PEZY Computing K.K.); CPU-only binary, built without PZSDK    <- CPU-only build
# git commit <commit ID>
# contact: https://github.com/t-oda-ic/hmmer_pz/issues
```

- The `# git commit` line shows the commit ID of the source that was built (§1). The same line goes into the header of result files,
  so you can tell later which version produced a result. The `# contact` line appears only in the `-h` output.
- `PEZY-SC3` is a fixed string, so it does not change with `PZ_ARCH`.

To check that the device can be used, run a small search with the bundled tutorial files.

```bash
cd <working directory>       # pzcl.log is written here
HMMER_PEZY_VERBOSE=1 <PREFIX>/bin/hmmsearch --cpu 4 \
    <SRC>/tutorial/globins4.hmm <SRC>/tutorial/globins45.fa > /dev/null
#  -> stderr shows the number of cards opened, e.g. "# pezy: 4 devices, ..."
#     If no "# pezy: <n> devices" line appears, the search ran on the CPU without the device
#     (timing lines such as "# pezy: [0.000s] start" appear even when the search runs on the CPU)
grep -E -i 'error|MMU' pzcl.log    # should print nothing
```

---

## 5. Run-time behavior and settings

### 5.1 Using multiple cards

By default every visible card is used. Cards that cannot be opened (for example, ones in use by another process)
are skipped and the run continues on the rest. Only when no card at all can be opened does the run use the CPU.

- `hmmsearch` / `phmmer` / `jackhmmer`: the target sequence database is split across the cards, and each card processes its share.
- `hmmscan`: the profile database is loaded onto every card, and models are handed to whichever card becomes free first.

The PZCL runtime log `pzcl.log` is written to the directory the program is run from.

### 5.2 When the CPU is used

In the following cases the programs run the same CPU code as the original instead of using the device (results do not change).

| Applies to | Condition |
|---|---|
| All | No device, every card is in use by other processes, or `HMMER_PEZY=0` is set |
| `hmmsearch` / `phmmer` / `jackhmmer` | The target sequence database is any of: not a regular FASTA file (pipe, gzip, Stockholm, ...), contains gap characters `-`, has a sequence longer than 100,000 residues, or the file is larger than `HMMER_PEZY_MAXDB` (default 8 GiB) |
| `hmmsearch` / `phmmer` | `--restrictdb_stkey` or `--restrictdb_n` is given |
| The three above (per query) | The query model length exceeds 2040 (only that query runs on the CPU) |
| `hmmscan` | Fewer query sequences than `HMMER_HMMCACHE_MINQ` (default 4), or the profile database is larger than `HMMER_HMMCACHE_MAX` (default 64 GiB) |

If a kernel does not finish within the time limit, the work is also redone on the CPU (a warning is printed to stderr):
`hmmsearch` / `phmmer` / `jackhmmer` redo only that query on the CPU, while `hmmscan` uses the CPU for all remaining queries.

### 5.3 Environment variables

| Variable | Default | Effect |
|---|---|---|
| `HMMER_PEZY=0` | — | Do not use the device; run the same CPU code as the original (for comparing results). Any value starting with `0`, `n` or `N` disables the device |
| `HMMER_PEZY_NDEV=<n>` | all visible cards | Upper limit on the number of cards used. `1` uses a single card |
| `HMMER_PEZY_VERBOSE=1` | — | Print the number of cards opened, each card's share of the work, and per-step timings to stderr (diagnostic) |
| `HMMER_PEZY_MAXDB=<bytes>` | 8 GiB | For `hmmsearch` etc., the largest target sequence database that is loaded into memory for the device |
| `HMMER_PEZY_RTHREADS=<n>` | `--cpu` + 1 | Number of threads used to read the target sequence database |
| `HMMER_HMMCACHE_MINQ=<n>` | 4 | `hmmscan`: with at least this many query sequences, the profile database is loaded into memory and the device is used |
| `HMMER_HMMCACHE_MAX=<bytes>` | 64 GiB | `hmmscan`: the largest profile database that is loaded into memory |
| `HMMER_HMMCACHE_QBATCH=<n>` | 512 | `hmmscan`: number of query sequences passed to the device at once |

Normally only `HMMER_PEZY`, `HMMER_PEZY_NDEV` and `HMMER_PEZY_VERBOSE` need to be changed.

---

## 6. Tests

```bash
cd <SRC>
make check        # the HMMER test suite (requires python3)
```

Build and run the unit tests for the PEZY parts as follows. Each checks that the result on the device (or of the accelerated reader)
matches the original CPU implementation exactly.

```bash
cd <SRC>/src
make pz_ssv_utest pz_ssvmulti_utest pz_fwdparse_utest pz_fpprobe_utest pz_pzdb_utest

./pz_ssv_utest      [-n <nseq>] <hmmfile> <seqfile>                 # SSV filter (for hmmsearch)
./pz_ssvmulti_utest [-n <nmodel>] [-t <threads>] <hmmfile> <seqfile> # many-model SSV filter (for hmmscan; <hmmfile> must be hmmpress'ed)
./pz_fwdparse_utest [-n <nseq>] <hmmfile> <seqfile>                 # parser-mode Forward / Backward (bit-for-bit against impl_sse)
./pz_fpprobe_utest  [-n <ntrial>]                                   # whether the device's float arithmetic matches x86 SSE bit for bit
./pz_pzdb_utest     [-n <nsample>] [-t <threads>] <seqfile>         # fast sequence database reader (compared with esl_sqio)
```

Running with different numbers of cards also checks how the database is split.

```bash
for nd in 1 2 3 4; do
  HMMER_PEZY_NDEV=$nd ./pz_ssv_utest <hmmfile> <seqfile> | tail -1
done
```

---

## 7. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| The line in (c) says `CPU-only binary` / `ldd` does not show `pzcl` | PZSDK was not found. Check `PZSDK_PATH`, and that `<PZSDK>/bin/pzcc`, `<PZSDK>/{include,inc}/pzcl/pzcl.h` and `<PZSDK>/lib/libpzcl.so` exist. After fixing, **start again from `make clean`** |
| Errors about `libstdc++.so.6` or `undefined reference to std::...` at link time | conda's gcc / ld is being used. `libpzcl.so` needs the system libstdc++. Remove conda from PATH as in §3 (1) |
| `libpzcl.so: cannot open shared object file` at run time | The SDK is in a different place at run time than at build time. The executables have `<PZSDK>/lib` as their rpath. Put the SDK at the same path, or set `LD_LIBRARY_PATH=<the SDK's lib>` |
| `make install` tries to write somewhere unexpected | A `Makefile` from an earlier build is left over. Run `make distclean` and configure again |
| `# git commit` shows `unknown-dirty` / a different ID | The tree was built without `.git` (for example, one obtained as a ZIP). Put a `GIT_COMMIT` file as described in §1, or build from a tree obtained with `git clone`. If the tree sits inside another git repository, that repository's ID is shown |
| `-dirty` is appended to the `# git commit` ID | Files tracked by git in the built tree have uncommitted changes |
| The CPU is used instead of the device | One of the conditions in §5.2 applies. If no `# pezy: <n> devices` line appears even with `HMMER_PEZY_VERBOSE=1`, the run is on the CPU (the reason is not shown) |
| Some cards are not used | Cards that cannot be opened are skipped. With `HMMER_PEZY_VERBOSE=1`, `skipping a device` is printed with the reason |
| A run hangs or terminates abnormally | Check `pzcl.log` in the directory the program was run from, e.g. with `grep "MMU ERROR"` |

---

## 8. License and contact

- The license is BSD 3-Clause, the same as the original HMMER. See `LICENSE`, which also contains the original HMMER copyright notices.
  The Easel license is in `easel/LICENSE`.
- Please report bugs and send questions via [GitHub Issues](https://github.com/t-oda-ic/hmmer_pz/issues).
  Please direct questions about this modified version here rather than to the original HMMER developers.
