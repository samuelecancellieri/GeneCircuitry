# ============================================================
# GeneCircuitry Docker Image
# Base: Ubuntu 22.04 (ships Python 3.10, within >=3.9,<3.11)
# Usage: docker run genecircuitry-docker genecircuitry --help
# ============================================================
FROM ubuntu:22.04

LABEL maintainer="Samuele Cancellieri <samuelc@uio.no>" \
      description="GeneCircuitry – transcriptional regulatory network analysis" \
      version="0.1.4"

# ── prevent interactive prompts during package installation ──────────────────
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ── system packages ──────────────────────────────────────────────────────────
# build tools, Python 3.10, HDF5, BLAS/LAPACK, Cairo/Pango (weasyprint),
# wkhtmltopdf (pdfkit fallback), git (some pip packages clone at build time)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        gfortran \
        make \
        cmake \
        git \
        curl \
        wget \
        ca-certificates \
        # Python 3.10
        python3.10 \
        python3.10-dev \
        python3.10-distutils \
        python3-pip \
        # HDF5 (h5py / anndata / celloracle)
        libhdf5-dev \
        pkg-config \
        # BLAS / LAPACK (numpy, scipy)
        libopenblas-dev \
        liblapack-dev \
        # igraph C library (python-igraph / leidenalg)
        libigraph-dev \
        # Cairo + Pango (weasyprint PDF rendering)
        libcairo2-dev \
        libpango1.0-dev \
        libpangocairo-1.0-0 \
        libgdk-pixbuf2.0-dev \
        libffi-dev \
        shared-mime-info \
        # wkhtmltopdf (pdfkit PDF rendering fallback)
        wkhtmltopdf \
        # misc
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── make python3.10 the default python/pip ───────────────────────────────────
RUN update-alternatives --install /usr/bin/python  python  /usr/bin/python3.10 1 \
 && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1 \
 && update-alternatives --install /usr/bin/pip     pip     /usr/bin/pip3       1

# ── upgrade pip / setuptools / wheel ─────────────────────────────────────────
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# ── Step 1: Cython + numpy (must precede packages that build C extensions) ───
RUN pip install --no-cache-dir \
        Cython \
        "numpy==1.26.4" 

RUN pip install --no-cache-dir \
        "velocyto==0.17.17"

# ── Step 2: celloracle (pinned; requires numpy<2 and Cython already present) ─
RUN pip install --no-cache-dir \
        "celloracle==0.18.0"

# ── Step 3: hotspotsc ────────────────────────────────────────────────────────
RUN pip install --no-cache-dir \
        "hotspotsc==1.1.3"

# ── Step 4: remaining runtime dependencies from requirements.txt ─────────────
RUN pip install --no-cache-dir \
        leidenalg \
        adjustText \
        gseapy \
        fa2-modified \
        weasyprint \
        pdfkit

# ── Step 5: copy source and install the package itself ───────────────────────
WORKDIR /opt/genecircuitry

COPY . .

RUN pip install --no-cache-dir --no-deps .

# ── runtime environment ───────────────────────────────────────────────────────
# Default output/data directories (users can bind-mount over these)
RUN mkdir -p /data /output
ENV GENECIRCUITRY_OUTPUT_DIR=/output

# ── entry point ───────────────────────────────────────────────────────────────
# docker run genecircuitry-docker genecircuitry --help
ENTRYPOINT ["genecircuitry"]
CMD ["--help"]
