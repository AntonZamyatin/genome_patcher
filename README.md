# genome_patcher

Patches a genome assembly using curated patch sequences. Assembly contigs are aligned to the patch sequences with minimap2, then the tool builds patch-centric chains and emits a corrected FASTA where patched regions are replaced or augmented with the appropriate assembly or patch sequence.

## How it works

1. **Align** — minimap2 aligns assembly contigs (query) against patch sequences (reference) and produces a sorted, indexed BAM.
2. **Patch** — alignments are parsed into hits, filtered by identity and length, and used to build chains along each patch sequence. Each chain is a recipe of segments (from the patch or the assembly) that together form the corrected output contig.
3. **Visualize** — per-patch hit maps and per-chain segment diagrams are written as PNG files.

### Chain linking

When two patch chains share an open endpoint on the same assembly contig (one chain's right edge and another's left edge land on the same contig), a bridge segment of assembly sequence connects them into a single output contig.

### Preserve mode (`-p`)

By default the chain body is built from assembly segments projected back from patch coordinates. With `-p`, the original patch sequence is kept in the chain body instead; only the flanks and bridges come from the assembly. Use this when you trust the patch sequence itself and want the assembly only to provide flanking and linking context.

---

## Requirements

### External tools

| Tool | Version tested | Notes |
|---|---|---|
| [minimap2](https://github.com/lh3/minimap2) | 2.30 | must be available in `PATH` |

### Python

Python >= 3.10

### Python packages

| Package | Version tested |
|---|---|
| [pysam](https://github.com/pysam-developers/pysam) | 0.23.3 |
| [biopython](https://biopython.org) | 1.87 |
| [matplotlib](https://matplotlib.org) | 3.10.8 |

---

## Installation

```bash
# 1. Clone the repository
git clone <repo-url>
cd genome_patcher

# 2. Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 3. Install Python dependencies
pip install pysam biopython matplotlib

# 4. Ensure minimap2 is available
minimap2 --version
```

---

## Usage

The tool has three subcommands: `align`, `patch`, and `run`.

```
python genome_patcher.py <command> [options]
```

---

### `run` — align and patch in one step

```bash
python genome_patcher.py run \
  --assembly <assembly.fasta> \
  --patch <patch.fasta> \
  -o <out_dir> \
  -l <label> \
  [options]
```

---

### `align` — alignment stage only

Runs minimap2 and writes a sorted, indexed BAM. Use this if you want to inspect the BAM before patching, or to reuse an existing alignment.

```bash
python genome_patcher.py align \
  --assembly <assembly.fasta> \
  --patch <patch.fasta> \
  -o <out_dir> \
  -l <label> \
  [options]
```

---

### `patch` — patching stage only

Reads an existing BAM and produces the patched FASTA, alignment table, and plots.

```bash
python genome_patcher.py patch \
  --assembly <assembly.fasta> \
  --patch <patch.fasta> \
  -o <out_dir> \
  -l <label> \
  [--bam <path/to/sorted.bam>] \
  [options]
```

---

## Options

### Input / output

| Option | Commands | Default | Description |
|---|---|---|---|
| `--assembly` | all | required | Assembly FASTA (`.fa` / `.fa.gz`), used as minimap2 query |
| `--patch` | all | required | Patch FASTA (`.fa` / `.fa.gz`), used as minimap2 reference |
| `-o / --out-dir` | all | required | Output directory (created if absent) |
| `-l / --label` | all | required | Prefix for all output files |
| `--bam` | `patch` | `<out-dir>/<label>.sorted.bam` | Pre-existing BAM to use instead of the default path |

### Alignment

| Option | Commands | Default | Description |
|---|---|---|---|
| `-t / --threads` | `align`, `run` | `16` | Threads for minimap2 and BAM sorting |
| `--tmp-dir` | `align`, `run` | `<out-dir>/tmp` | Directory for intermediate SAM/BAM files |
| `--keep-unmapped` | `align`, `run` | off | Keep unmapped contigs in the output BAM (default: drop them) |

### Filtering

| Option | Commands | Default | Description |
|---|---|---|---|
| `--min-pid` | `patch`, `run` | `0.999` | Minimum alignment identity. Identity is computed as `matches / (matches + mismatches + insertions + deletions)`, penalising all forms of sequence divergence including structural indels |
| `--min-aln-len` | `patch`, `run` | `1000` | Minimum aligned length in bases (matches + mismatches, indels excluded) |
| `--border-tolerance` | `patch`, `run` | `10` | Bases from each patch end within which a terminal hit is considered to reach the border and triggers chain linking |

### Patching mode

| Option | Commands | Default | Description |
|---|---|---|---|
| `-p / --preserve-mode` | `patch`, `run` | off | Keep the patch sequence in the chain body; only flanks and bridge segments come from the assembly |

### General

| Option | Commands | Description |
|---|---|---|
| `-v / --verbose` | all | Print progress details to stdout |

---

## Outputs

All files are written to `<out-dir>/` with the prefix `<label>`.

| File | Description |
|---|---|
| `<label>.sorted.bam` | Sorted, indexed minimap2 alignment (`align` / `run`) |
| `<label>.patched.fasta` | Patched assembly: emitted chains + untouched contigs |
| `<label>.changed.txt` | One assembly contig name per line — contigs consumed by emitted chains |
| `<label>.alignments.tsv` | Tab-separated hit table (see below) |
| `patch_plots/<label>.<patch_id>.all_hits.png` | All filtered hits on a given patch sequence |
| `patch_plots/<label>.<chain_id>.patched_chain.png` | Segment layout of each emitted chain |

### Alignment table columns

| Column | Description |
|---|---|
| `patch_seq_name` | Patch sequence name |
| `ctg_name` | Assembly contig name |
| `p_start` / `p_end` | Hit coordinates on the patch sequence |
| `q_start` / `q_end` | Hit coordinates on the assembly contig |
| `strand` | `+` forward, `-` reverse-complement |
| `pid` | Alignment identity (indel-aware) |
| `plen` | Percent of the assembly contig covered by this hit |

---

## Example

```bash
python genome_patcher.py run \
  --assembly reference.hpc.fasta \
  --patch patch_seqs.fasta \
  -o out \
  -l mPleAur1 \
  --min-pid 0.999 \
  --min-aln-len 1000 \
  --border-tolerance 10 \
  -t 16 \
  -p \
  --verbose
```