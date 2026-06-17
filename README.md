# kRID
kmer based redundancy and identity detection of germplasm samples.
# DEDUKT

**DEDUplication by K-mer Testing** — an alignment-free pipeline for detecting
duplicate and redundant accessions across plant germplasm collections from
whole-genome sequencing data.

DEDUKT compares each query sample against a reference panel using shared
k-mer content (via [KMC](https://github.com/refresh-bio/KMC)), and is built
specifically to stay robust when the two collections were sequenced at very
different depths — the situation where naïve whole-genome similarity breaks
down. It was developed to reconcile a low-coverage USDA teff (*Eragrostis tef*)
collection against a high-coverage Ethiopian (EIAR) reference panel, but the
method is collection- and species-agnostic.

---

## Why it exists

When one collection is sequenced at ~15× and another at ~4×, the shallow
samples recover only a fraction of their k-mers. This caps symmetric metrics
(Jaccard) and makes the conserved, species-wide "core" dominate any
containment score, so genuine duplicates and unrelated samples look alike.
DEDUKT addresses this with two ideas:

1. **Directional containment** — score what fraction of the *smaller* (shallow)
   sample's k-mers are present in the reference accession, rather than a
   symmetric overlap that the depth gap caps.
2. **Variable-fraction core-masking** — restrict scoring to k-mers that are
   *not* shared by almost every accession in the panel (the "variable" or
   polymorphic fraction). This is the k-mer analogue of dropping monomorphic
   sites and genotyping on polymorphic markers only, and it removes the
   conserved core that otherwise inflates similarity.

Both legacy whole-set metrics and the core-masked metrics are reported side by
side, so you can see the effect of masking on identical input.

---

## Repository contents

| File | Purpose |
|------|---------|
| `v6.sh` | Main pipeline. Self-contained SLURM orchestrator (panel build → union + variable mask → per-query scoring → aggregation). |
| `aggregate.sh` | Cohort-level aggregation of per-query results into a single report. Supports `--exclude-self` for all-vs-all runs. |
| `plot_score_distributions.py` | Jaccard & containment distributions (whole-set vs core-masked). |
| `cluster_usda_eiar.py` | Clustered heatmap + MDS embedding of query↔panel similarity. |
| `analyze_low_scorers.py` | Diagnose low-scoring samples by sequencing depth and (with metadata) by species / QC fields. |

---

## Requirements

- **KMC ≥ 3.0** (`kmc`, `kmc_tools`) on `PATH`
- **SLURM** (`sbatch`) — the pipeline submits a dependent 4-phase job chain
- **bash ≥ 4**, **gawk**, coreutils
- **Python ≥ 3.8** with `pandas`, `numpy`, `scipy`, `matplotlib` (for the
  analysis/plotting scripts only)

> **Cluster note:** the `#SBATCH -p jic-medium,nbi-medium` partition lines in
> `v6.sh` are specific to the John Innes Centre cluster. Edit them (and the
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
nothing outside it is touched.

When the job chain finishes, the headline result is in
`my_workspace/cohort_report.txt`.

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

Example — sweep the masking threshold into separate workspaces:

```bash
for f in 0.80 0.90 0.95 0.99; do
  CORE_FREQ_FRAC=$f bash v6.sh --ignore-corrupt panel_dir/ query_dir/ ws_f${f/./}
done
```

---

## Output layout

```
my_workspace/
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

It resolves columns by header name, so it works on both legacy and
variable-fraction workspaces.

### All-vs-all (within-collection duplicates)

To find duplicates *within* one collection, use it as both panel and query, and
exclude the trivial self-matches:

```bash
PANEL_CI=2 bash v6.sh --ignore-corrupt collection_dir/ collection_dir/ self_run
bash aggregate.sh --exclude-self self_run
```

`--exclude-self` drops rows where the panel accession equals the query sample,
so each sample's reported top hit is its best *non-self* match.

---

## Analysis & plotting

```bash
# Score distributions: whole-set vs core-masked Jaccard and containment
python3 plot_score_distributions.py cohort_top_hits.tsv score_distributions.png

# Clustered heatmap + MDS embedding of query↔panel similarity
python3 cluster_usda_eiar.py cohort_all_hits.tsv jaccard_var cluster

# Diagnose low scorers by depth, and by species/QC if metadata is supplied
python3 analyze_low_scorers.py cohort_top_hits.tsv --meta metadata.tsv
```

`metadata.tsv` is a tab-separated table keyed by a `sample` column, with any of
`species`, `coverage`, `n_reads`, `pct_dups`, … (e.g. assembled from GRIN
passport data and/or MultiQC). When a `species` column is present, the script
splits containment by species — useful for confirming that wild/progenitor
accessions, not bad data, drive low scores.

---

## Typical workflows

1. **Cross-collection reconciliation** — low-coverage queries vs a high-coverage
   panel: `bash v6.sh --ignore-corrupt panel/ query/ recon` then read the
   VARIABLE-FRACTION section of `recon/cohort_report.txt`.
2. **Positive control** — run a collection with *known* duplicates against the
   panel and confirm those pairs score high on `containment_var` (calibrates
   the masking threshold).
3. **Within-collection deduplication** — all-vs-all with `--exclude-self`.
4. **Masking sensitivity** — sweep `CORE_FREQ_FRAC`.

---

## Notes & caveats

- The variable-mask build (Phase 2) binarises every panel accession and sums
  occurrence across the panel; it needs transient scratch roughly the size of
  the panel and is the most compute-heavy step. It parallelises over
  `MASK_PAR` workers and degrades cleanly to legacy-only metrics if it fails.
- Counter caps: the per-k-mer occurrence sum uses `-cs` sized to the panel; for
  panels larger than ~250 accessions check that KMC counter limits are adequate.
- `aggregate.sh --exclude-self` matches self by exact `sample == panel` string
  equality; ensure panel and query naming are identical in all-vs-all runs.

---



## License

Add a license of your choice (e.g. MIT) before publishing.
