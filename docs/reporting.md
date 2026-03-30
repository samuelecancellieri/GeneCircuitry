---
layout: default
title: Reporting
nav_order: 6
parent: User Guide
---

# Reporting

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

GeneCircuitry generates **HTML and PDF reports** that consolidate figures, statistics, and pipeline metadata into a single shareable document. Reports are produced automatically at the end of each pipeline run and can also be generated programmatically.

The reporting engine lives in `genecircuitry/reporting/`. The main entry point is `generate_report()` from `genecircuitry/reporting/__init__.py`.

---

## `generate_report()`

```python
from genecircuitry.reporting import generate_report

outputs = generate_report(
    output_dir="results/",            # directory for all report outputs
    title="My Analysis Report",       # main title
    subtitle="",                      # optional subtitle line
    adata=adata,                      # processed AnnData (for statistics section)
    celloracle_result=(oracle, links), # tuple (oracle, links), or None to skip GRN section
    hotspot_result=hs,                # Hotspot object, or None to skip module section
    log_file="results/logs/pipeline.log",  # path to pipeline.log for step summary
    embed_images=True,                # True = inline base64 images (portable HTML)
                                     # False = relative paths (smaller file, requires figures/)
    formats=["html", "pdf"],         # list of output formats
)

# outputs: dict with paths
# {"html": "results/report.html", "pdf": "results/report.pdf"}
print(outputs)
```

### Parameters

| Parameter           | Type              | Default                           | Description                                                        |
| ------------------- | ----------------- | --------------------------------- | ------------------------------------------------------------------ |
| `output_dir`        | `str`             | required                          | Root directory — report is saved here                              |
| `title`             | `str`             | `"GeneCircuitry Analysis Report"` | Report title                                                       |
| `subtitle`          | `str`             | `""`                              | Subtitle (e.g. sample name or condition)                           |
| `adata`             | `AnnData`         | `None`                            | Used to populate the data summary section                          |
| `celloracle_result` | `tuple \| None`   | `None`                            | `(oracle, links)` for GRN section                                  |
| `hotspot_result`    | `Hotspot \| None` | `None`                            | Hotspot object for modules section                                 |
| `log_file`          | `str \| None`     | `None`                            | Path to `pipeline.log` for step timeline                           |
| `embed_images`      | `bool`            | `True`                            | Embed figures as base64 (portable) or use relative paths (smaller) |
| `formats`           | `list[str]`       | `["html"]`                        | `"html"` and/or `"pdf"`                                            |

### Return value

A `dict` mapping format name to output path:

```python
{"html": "results/report.html", "pdf": "results/report.pdf"}
```

---

## HTML vs PDF reports

| Feature                | HTML                                  | PDF                           |
| ---------------------- | ------------------------------------- | ----------------------------- |
| Interactive navigation | ✅ collapsible sections, anchor links | ❌ static                     |
| Portable (single file) | ✅ when `embed_images=True`           | ✅                            |
| File size              | Larger with embedded images           | Smaller                       |
| Requires extra dep     | No                                    | Yes: `pip install weasyprint` |
| Best for               | Sharing with collaborators, archiving | Printing, submission          |

{: .note }

> PDF generation requires [WeasyPrint](https://weasyprint.org/). Install it with `pip install weasyprint`. If not available, PDF is skipped and a warning is logged.

---

## Report sections

A GeneCircuitry report automatically includes:

| Section               | Contents                                                  | Requires             |
| --------------------- | --------------------------------------------------------- | -------------------- |
| **Analysis settings** | Key config parameters, pipeline version, timestamp        | always               |
| **Data summary**      | Cell count, gene count, cluster breakdown, QC metrics     | `adata`              |
| **QC**                | Before/after filtering metrics, embedded QC plots         | `adata` + QC figures |
| **Clustering**        | UMAP plot, cluster sizes, leiden resolution used          | `adata`              |
| **GRN inference**     | Per-cluster TF network statistics, top TFs, GRN plots     | `celloracle_result`  |
| **Gene modules**      | Module count, genes per module, local correlation heatmap | `hotspot_result`     |
| **Pipeline log**      | Step-by-step timeline with durations                      | `log_file`           |

---

## Using `embed_images`

**`embed_images=True`** (default):

- All figures are read from disk and encoded as base64 strings inside the HTML
- The resulting `.html` file is fully self-contained — open it anywhere without needing the `figures/` directory alongside it
- File size increases proportionally to the number and resolution of figures

**`embed_images=False`**:

- The HTML references figures with relative paths (e.g. `figures/qc/violin_pre_filter.png`)
- The HTML must be kept in `output_dir/` alongside the `figures/` directory to render correctly
- Much smaller file size — suitable for web hosting or sharing via a directory link

---

## `ReportGenerator` class

For programmatic control over report sections, use `ReportGenerator` directly:

```python
from genecircuitry.reporting.generator import ReportGenerator
from genecircuitry.reporting.sections import (
    create_data_summary_section,
    create_settings_section,
    create_celloracle_section,
    create_hotspot_section,
)

# Initialise generator
gen = ReportGenerator(
    output_dir="results/",
    title="Custom Report",
    subtitle="TypeA cells only",
)

# Add sections
gen.add_section(create_settings_section())
gen.add_section(create_data_summary_section(adata))
gen.add_section(create_celloracle_section(oracle, links))
gen.add_section(create_hotspot_section(hs))

# Generate outputs
gen.generate_html("report.html", embed_images=True)
gen.generate_pdf("report.pdf")
```

### `ReportSection` dataclass

Each section is a `ReportSection` with:

| Field         | Type                  | Description                         |
| ------------- | --------------------- | ----------------------------------- |
| `title`       | `str`                 | Section heading                     |
| `content`     | `str`                 | HTML or Markdown content            |
| `section_id`  | `str`                 | Anchor ID (auto-derived from title) |
| `figures`     | `list[str]`           | File paths to embed                 |
| `tables`      | `list[dict]`          | Tabular data                        |
| `metrics`     | `dict`                | Key-value statistics                |
| `subsections` | `list[ReportSection]` | Nested subsections                  |

---

## Example: report for stratified analysis

```python
from genecircuitry.reporting import generate_report

for cluster_name, cluster_data in stratification_results.items():
    generate_report(
        output_dir=f"results/stratified_analysis/{cluster_name}/",
        title=f"GRN Analysis — {cluster_name}",
        adata=cluster_data["adata"],
        celloracle_result=cluster_data.get("celloracle"),
        hotspot_result=cluster_data.get("hotspot"),
        log_file=f"results/stratified_analysis/{cluster_name}/logs/pipeline.log",
        formats=["html"],
    )
```

---

## Report in the pipeline

When running via the CLI, a report is automatically generated as the final step:

```bash
python -m genecircuitry.pipeline --input data.h5ad --output results/
# → results/report.html (and report.pdf if weasyprint installed)
```

To generate a report without running the full pipeline:

```bash
python -m genecircuitry.pipeline \
    --input results/clustered_adata.h5ad \
    --output results/ \
    --steps report
```
