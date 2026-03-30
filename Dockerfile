# ============================================================
# GeneCircuitry Docker Image
# Base: condaforge/miniforge3 (conda-forge + bioconda ready)
# Usage: docker run samuelecancellieri/genecircuitry --help
# ============================================================
FROM condaforge/miniforge3:latest

LABEL maintainer="Samuele Cancellieri <samuelc@uio.no>" \
      description="GeneCircuitry – transcriptional regulatory network analysis" \
      version="0.1.6"

# ── Step 1: all conda-installable dependencies ───────────────────────────────
# Mirrors pixi.toml [dependencies]. No strict channel priority to avoid
# solver conflicts; conda-forge + bioconda covers everything listed here.
RUN conda install -y -n base \
        -c conda-forge \
        -c bioconda \
        "python>=3.9,<3.11" \
        "numpy>=1.20.0,<2.0.0" \
        "pandas>=1.0.3,<=1.5.3" \
        "scanpy>=1.9.0" \
        "anndata>=0.8.0" \
        "matplotlib>=3.6.3" \
        "seaborn>=0.11.0" \
        "scipy>=1.7.0" \
        "networkx>=2.6.0" \
        "leidenalg" \
        "adjusttext>=0.7.3" \
        "pip" \
    && conda clean -afy

# ── Step 2: pip-only packages (PyPI only — not on conda-forge/bioconda) ───────
RUN pip install --no-cache-dir \
        "fa2-modified" \
        "celloracle==0.18.0" \
        "hotspotsc==1.1.3"

# ── Step 3: optional reporting dependencies ───────────────────────────────────
RUN pip install --no-cache-dir \
        weasyprint

# ── Step 4: copy source and install genecircuitry itself ──────────────────────
WORKDIR /opt/genecircuitry

COPY . .

RUN pip install --no-cache-dir --no-deps .

# ── Runtime directories (bind-mount your data and output here) ────────────────
RUN mkdir -p /data /output
ENV GENECIRCUITRY_OUTPUT_DIR=/output

# ── Entry point ───────────────────────────────────────────────────────────────
ENTRYPOINT ["genecircuitry"]
CMD ["--help"]

