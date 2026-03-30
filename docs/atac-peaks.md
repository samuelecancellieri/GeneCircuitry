---
layout: default
title: ATAC Peaks Processing
nav_order: 3
parent: Reference
---

# ATAC Peaks Processing

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

The `genecircuitry/atac_peaks_processing.py` module converts **ATAC-seq peak BED files** into a TF motif matrix (TF info PKL) that can be used as a **chromatin-informed base GRN** for CellOracle. This replaces the default promoter-only base GRN with open-chromatin-derived TF binding sites.

{: .note }

> This module requires the `atac` optional dependency group: `pip install -e ".[atac]"`. Dependencies: `celloracle`, `genomepy`, `gimmemotifs`.

---

## When to use ATAC peaks

By default, CellOracle builds GRNs using a **promoter-only base GRN** derived from publicly available ChIP-seq data. This is appropriate for most scRNA-seq datasets.

If you have **matched scATAC-seq data** (or bulk ATAC-seq from the same cell type), you can build a **chromatin-aware base GRN** by scanning open chromatin peaks for TF binding motifs. This approach:

- Captures cell-type-specific regulatory elements beyond promoters
- Uses enhancers and other distal regulatory regions
- Improves GRN specificity for rare populations with unique open chromatin landscapes

---

## `process_atac_peaks()`

The main function — runs the full BED → TF motif matrix pipeline.

```python
from genecircuitry.atac_peaks_processing import process_atac_peaks

tf_info_path = process_atac_peaks(
    bed_path="data/atac_peaks.bed",   # path to BED file with ATAC peaks
    species="human",                   # "human" (hg38) or "mouse" (mm10)
    output_dir=None,                   # defaults to config.OUTPUT_DIR
    fpr=None,                          # false positive rate for motif scanning
                                       # defaults to config.ATAC_MOTIF_SCAN_FPR
    motif_score_threshold=None,        # minimum motif score for filtering
                                       # defaults to config.ATAC_MOTIF_SCORE_THRESHOLD
)
# tf_info_path: path to the saved PKL file
```

### Parameters

| Parameter               | Type            | Default                             | Description                                 |
| ----------------------- | --------------- | ----------------------------------- | ------------------------------------------- |
| `bed_path`              | `str`           | required                            | Path to BED file with pre-called ATAC peaks |
| `species`               | `str`           | `"human"`                           | Species: `"human"` → hg38, `"mouse"` → mm10 |
| `output_dir`            | `str \| None`   | `config.OUTPUT_DIR`                 | Directory for output files                  |
| `fpr`                   | `float \| None` | `config.ATAC_MOTIF_SCAN_FPR`        | False positive rate for motif hit calling   |
| `motif_score_threshold` | `int \| None`   | `config.ATAC_MOTIF_SCORE_THRESHOLD` | Minimum motif score to keep                 |

### Return value

Path (str) to the saved TF info pickle file (e.g. `output/atac/enriched_motifs.pkl`).

---

## Full pipeline workflow

```
BED file (ATAC peaks)
    │
    ▼ (1) _annotate_bed_peaks()
TSS annotation CSV (peak → nearest gene)
    │
    ▼ (2) _ensure_genome_installed()
Reference genome (hg38 / mm10)
    │
    ▼ (3) CellOracle motif_analysis.TFinfo_from_bed()
Raw motif matrix
    │
    ▼ (4) Filter by motif_score_threshold
Enriched TF motif matrix
    │
    ▼ (5) Save as PKL
tf_info.pkl  →  passed to create_oracle_object(TG_to_TF_dictionary=...)
```

---

## Usage with CellOracle

After generating the TF motif matrix, pass it to [`create_oracle_object()`](../user-guide/celloracle/#step-1--create-oracle-object-create_oracle_object):

```python
from genecircuitry.atac_peaks_processing import process_atac_peaks
from genecircuitry.celloracle_processing import (
    perform_grn_pre_processing,
    create_oracle_object,
    run_PCA, run_KNN, run_links,
)

# Step 1: convert ATAC BED → TF motif matrix
tf_info_path = process_atac_peaks(
    bed_path="data/scatac_peaks.bed",
    species="human",
    output_dir="results/atac/",
)
print(f"TF motif matrix saved to: {tf_info_path}")

# Step 2: use in CellOracle (instead of default promoter base GRN)
adata_grn = perform_grn_pre_processing(adata, cluster_key="leiden")
oracle = create_oracle_object(
    adata_grn,
    cluster_column_name="leiden",
    embedding_name="X_umap",
    TG_to_TF_dictionary=tf_info_path,   # ← ATAC-informed base GRN
    raw_count_layer="raw_count",
)
# Continue with run_PCA, run_KNN, run_links as normal...
```

---

## BED file format

The input BED file must have at least 3 columns (chrom, start, end):

```
chr1    10000    10500
chr1    15000    15800
chr2    23000    23450
...
```

Standard 6-column BED format is also supported (with strand, score, name fields). The file should contain pre-called peaks — GeneCircuitry does not call peaks from raw ATAC data.

---

## Reference genome installation

The first time you run `process_atac_peaks()` for a given species, GeneCircuitry checks whether the reference genome is installed and downloads it if not:

```
Checking hg38 installation...
  hg38 installation: False
  Installing genome hg38...
  ✓ Genome hg38 installed
```

Genome files are stored in the default `genomepy` directory (`~/.local/share/genomes/`). Subsequent runs reuse the cached genome.

---

## Configuration reference

| Config parameter             | Default | Description                                 |
| ---------------------------- | ------- | ------------------------------------------- |
| `ATAC_MOTIF_SCAN_FPR`        | `0.02`  | False positive rate for TF motif scanning   |
| `ATAC_MOTIF_SCORE_THRESHOLD` | `10`    | Minimum motif score to retain in the matrix |

Adjust these globally before processing:

```python
from genecircuitry import config
config.update_config(
    ATAC_MOTIF_SCAN_FPR=0.01,        # stricter FPR → fewer false positives
    ATAC_MOTIF_SCORE_THRESHOLD=15,   # higher score → higher-confidence motifs
)
```

---

## Pipeline integration

ATAC peak processing is integrated as an optional step in the main pipeline:

```bash
python -m genecircuitry.pipeline \
    --input data.h5ad \
    --output results/ \
    --steps load preprocessing clustering atac_peaks celloracle
```

The pipeline expects the BED file path to be configured in `config.ATAC_BED_PATH` (set it via `config.update_config(ATAC_BED_PATH="path/to/peaks.bed")` before running) or passed at runtime.

---

## Troubleshooting

**`ImportError: genomepy`** — install with `pip install -e ".[atac]"`.

**`ImportError: gimmemotifs`** — same: `pip install -e ".[atac]"`. Note: `gimmemotifs` may require additional system libraries (see [gimmemotifs docs](https://gimmemotifs.readthedocs.io/)).

**Genome download fails** — check internet connection and disk space. Genomes are several GB. Set a custom genome directory via `genomepy.config["genome_dir"] = "/path/to/genomes"`.

**No TF motifs found after filtering** — lower `ATAC_MOTIF_SCORE_THRESHOLD` or increase `ATAC_MOTIF_SCAN_FPR`. Very narrow peak sets may have few high-confidence motif hits.
