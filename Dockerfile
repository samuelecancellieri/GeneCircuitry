# ============================================================
# GeneCircuitry Docker Image
# Base: condaforge/miniforge3 (conda-forge + bioconda ready)
# Usage: docker run samuelecancellieri/genecircuitry --help
# ============================================================
FROM condaforge/miniforge3:latest

LABEL maintainer="Samuele Cancellieri <samuelc@uio.no>" \
      description="GeneCircuitry – transcriptional regulatory network analysis" \
      version="0.1.6"

# ── Activate strict channel priority for reproducibility ─────────────────────
# RUN conda config --system --set channel_priority strict

# ── Step 1: conda-available core dependencies ────────────────────────────────
# All packages listed here are available on conda-forge or bioconda.
# celloracle and hotspotsc are pip-only and installed separately below.
RUN conda install -y -n base \
        -c bioconda \
        -c conda-forge \
        "python>=3.9,<3.11" \
        "numpy>=1.20.0,<2.0.0" \
        "pandas>=1.3.0" \
        "scanpy>=1.9.0" \
        "anndata>=0.8.0" \
        "matplotlib-base>=3.6.3" \
        "seaborn>=0.11.0" \
        "scipy>=1.7.0" \
        "networkx>=2.6.0" \
        leidenalg \
        pip \
    && conda clean -afy

# ── Step 2: pip-only dependencies (not on conda-forge / bioconda) ─────────────
RUN pip install --no-cache-dir \
        "adjustText>=0.7.3" \
        "fa2-modified" \
        "celloracle==0.18.0" \
        "hotspotsc==1.1.3"

# ── Step 3: optional reporting dependencies (pip) ─────────────────────────────
RUN pip install --no-cache-dir \
        weasyprint \
        pdfkit

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

