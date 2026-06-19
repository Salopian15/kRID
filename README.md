# kRID
kmer based redundancy and identity detection of germplasm samples.
# kRID
Alignment-free pipeline for detecting duplicate and redundant accessions across plant g
ermplasm collections from genome sequencing data.

kRID compares each query sample against a reference panel using shared
k-mer content (via [KMC](https://github.com/refresh-bio/KMC)). It was developed to reconcile the multiple
different *Eragrostis tef* germplasm databases currently available to researchers, 
namely the USDA, Ethiopian (EIAR) and IBERS sample panels, but the
method is collection- and species-agnostic.

---

## Why it exists

Germplasm banks containing duplicate samples and use different naming schemes
compared to other banks is a known issue and makes it difficult to assess whether
a sample being used in one location with one name attached to it is the same sample
labelled as such in another database.


## Repository contents

| File | Purpose |
|------|---------|
| `v6.sh` | Main pipeline. Self-contained SLURM orchestrator (panel build → union + variable mask → per-query scoring → aggregation). |
| `aggregate.sh` | Cohort-level aggregation of per-query results into a single report. Supports `--exclude-self` for all-vs-all runs. |

---

## Requirements

- **KMC ≥ 3.0** (`kmc`, `kmc_tools`) on `PATH`, slurm modules should work fine too
- **SLURM** (`sbatch`) — the pipeline submits a dependent 4-phase job chain, currently hard-coded to the JIC slurm node names but easy to change
- **bash ≥ 4**, **gawk**, coreutils

> **Cluster note:** the `#SBATCH -p jic-medium,nbi-medium` partition lines in
> `v6.sh` are specific and won't work outside of the JIC cluster. Edit them (and the
> time/memory/core requests) to match your scheduler before running elsewhere.

---

## Quick start

Put `v6.sh` and `aggregate.sh` in the same directory (the pipeline calls
`aggregate.sh` automatically as its final phase).

```bash
# panel = precompiled KMC databases or FASTQ/FASTA; query = FASTQ/FASTA/KMC
bash v6.sh --ignore-corrupt  panel_dir/  query_dir/  my_workspace
```

Three positional arguments: `<panel_input_dir> <query_input_dir> <project_workspace>`.
Use a fresh `<project_workspace>` per run — results are written there and
nothing outside it is touched, older results may be overwritten too.

When the job chain finishes, the headline result is in
`project_workspace/cohort_report.txt`.

---

## Input formats

- **Panel** (`panel_dir/`): a directory of **precompiled KMC databases**
  (`*.kmc_pre` / `*.kmc_suf`, symlinked in), and/or **FASTQ/FASTA** that the
  pipeline will build into KMC databases.
- **Query** (`query_dir/`): **FASTQ** (paired or single), **FASTA**, or
  precompiled **KMC** databases.

Sample names are taken from file basenames.

---

## Parameters (environment variables)

All optional; defaults match the low-coverage reconciliation use case.

| Variable | Default | Meaning |
|----------|---------|---------|
| `QUERY_CI` | `2` | k-mer min-count for FASTQ **queries** (drop sequencing-error k-mers). Use 2 for low coverage, higher for deep data. |
| `PANEL_CI` | `4` | k-mer min-count for FASTQ **panel** builds (ignored for precompiled KMC panels). |
| `MASK_CORE` | `1` | Build the variable-fraction mask and report `*_var` metrics (`0` = legacy/whole-set only). |
| `CORE_FREQ_FRAC` | `0.95` | A k-mer present in ≥ this fraction of panel accessions is "core" and masked out; the rest form the variable set. |
| `MASK_PAR` | `12` | Parallel `kmc_tools` workers used to build the mask (Phase 2). |
| `JACCARD_THRESHOLD` | `0.85` | Jaccard cutoff for the `REDUNDANT` flag. |
| `CONTAINMENT_THRESHOLD` | `0.95` | Containment cutoff for the `REDUNDANT`/`SUBSET` flag. |
| `MIN_PANEL_KMERS` | `1000000` | Panel accessions with fewer k-mers are excluded by a sanity check. |


---

## Output layout

```
project_workspace/
├── panel_union/
│   ├── unified_teff_kmers.*        # union of all panel k-mers
│   ├── accession_sizes.tsv         # |k-mers| per panel accession
│   ├── panel_variable.*            # the variable-fraction mask (if MASK_CORE=1)
│   ├── accession_sizes_var.tsv     # |variable k-mers| per accession
│   ├── panel_union_size.txt        # |U|
│   └── panel_excluded.tsv          # accessions failing the sanity check
├── query_runs/<sample>/
│   ├── summary.txt / summary.tsv           # ranked by whole-set Jaccard
│   ├── summary_var.txt / summary_var.tsv   # ranked by variable-fraction Jaccard
│   └── rates.tsv / rates_var.tsv           # full per-panel scores for that query
├── cohort_top_hits.tsv             # best hit per query (whole-set ranking)
├── cohort_all_hits.tsv             # all per-query top-N hits concatenated
├── cohort_panel_frequency.tsv      # how often each panel accession is the top hit
├── cohort_top_hits_var.tsv         # best hit per query (variable-fraction ranking)
└── cohort_report.txt               # human-readable cohort summary (LEGACY + VARIABLE sections)
```

### Metrics

For a query *B* against a panel accession *A* (`ni = |A∩B|`, `na = |A|`, `nb = |B|`):

- **`jaccard`** = `ni / (na + nb − ni)` — symmetric overlap (depth-sensitive).
- **`containment`** = `ni / min(na, nb)` — directional; for a shallow query
  this is "fraction of the query's k-mers found in the reference."
- **`jaccard_var` / `containment_var`** — the same, computed only on the
  variable (core-masked) k-mer set.
- **`ssr`** (shared-state rate) — a simple-matching coefficient over the panel
  universe; diagnostic only (it counts joint absences and is not used for flags).
- **`flag` / `flag_var`** — `REDUNDANT` (`jaccard ≥ JACCARD_THRESHOLD` **and**
  `containment ≥ CONTAINMENT_THRESHOLD`), `SUBSET` (`containment ≥` threshold,
  `jaccard <` threshold), or `.`.

> In a low-coverage cross-collection comparison, `containment` (and especially
> `containment_var`) is the reliable duplicate signal; `jaccard` mostly reports
> how completely the query was sequenced and rarely reaches the `REDUNDANT` bar.

---

## Aggregation

`aggregate.sh` runs automatically as Phase 4, but you can re-run it any time:

```bash
bash aggregate.sh my_workspace
```

Will combine results across the per accession results into easier to read files.

---
