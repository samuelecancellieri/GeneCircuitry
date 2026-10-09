"""
ATAC Peaks Processing Module for GeneCircuitry
==========================================

Processes ATAC-seq peak BED files through CellOracle's motif analysis
pipeline to generate enriched motif matrices (TF info) that can be used
as custom base GRN input for CellOracle GRN inference.

Typical workflow:
    1. Read a BED file with pre-called ATAC peaks
    2. Annotate peaks with TSS / gene info
    3. Save annotated peaks as CSV
    4. Ensure reference genome is installed
    5. Scan annotated peaks for TF motifs
    6. Filter motifs by score
    7. Export enriched motif matrix as a PKL file

Usage within GeneCircuitry pipeline:
    from genecircuitry.atac_peaks_processing import process_atac_peaks
    tf_info_path = process_atac_peaks("peaks.bed", species="human")
"""

import os
import pickle
from typing import Optional, Union, Tuple, Any

import pandas as pd
from anndata import AnnData

from genecircuitry import config


def _extract_raw_counts(
    adata: AnnData,
    raw_layer: Optional[str] = None,
) -> AnnData:
    """
    Ensure the scATAC AnnData has raw counts in .X.

    If raw_layer is specified:
    - If found in adata.layers, copies adata.layers[raw_layer] into adata.X.
    - If raw_layer in ('raw', '.raw') and adata.raw is not None, extracts counts from adata.raw.
    - Otherwise raises ValueError.

    If raw_layer is not specified (default None):
    - Automatically checks adata.layers for standard raw count names
      ('raw_counts', 'raw_count', 'raw', 'counts', 'count').
    - If no candidate layer matches, checks if adata.raw is present.
    - If any other layer contains 'raw' or 'count', falls back to that layer.
    - If none found, keeps adata.X as is with an informational note.

    Parameters
    ----------
    adata : AnnData
        Loaded scATAC AnnData object.
    raw_layer : Optional[str], default None
        Explicit layer name or 'raw' / '.raw'.

    Returns
    -------
    AnnData
        AnnData object with raw counts in .X.
    """
    if raw_layer is not None:
        raw_layer_clean = raw_layer.strip()
        if raw_layer_clean in adata.layers:
            print(f"  ✓ Using raw counts from specified layer '{raw_layer_clean}'")
            adata.X = adata.layers[raw_layer_clean].copy()
            return adata

        lower_layers = {k.lower(): k for k in adata.layers.keys()}
        if raw_layer_clean.lower() in lower_layers:
            matched = lower_layers[raw_layer_clean.lower()]
            print(f"  ✓ Using raw counts from specified layer '{matched}'")
            adata.X = adata.layers[matched].copy()
            return adata

        if raw_layer_clean.lower() in ("raw", ".raw"):
            if adata.raw is None:
                raise ValueError(
                    f"Specified raw layer '{raw_layer_clean}', but scATAC AnnData has no .raw."
                )
            if set(adata.var_names).issubset(set(adata.raw.var_names)):
                adata.X = adata.raw[:, adata.var_names].X.copy()
            else:
                raw_ad = adata.raw.to_adata()
                for k, v in adata.obsm.items():
                    if k not in raw_ad.obsm:
                        raw_ad.obsm[k] = v
                adata = raw_ad
            print(f"  ✓ Using raw counts from .raw as requested ('{raw_layer_clean}')")
            return adata

        raise ValueError(
            f"Specified raw layer '{raw_layer_clean}' not found in scATAC data. "
            f"Available layers: {list(adata.layers.keys())}, has .raw: {adata.raw is not None}."
        )

    # Auto-detection mode (raw_layer is None)
    candidate_layers = ["raw_counts", "raw_count", "raw", "counts", "count"]
    lower_layers = {k.lower(): k for k in adata.layers.keys()}
    for cand in candidate_layers:
        if cand in lower_layers:
            matched = lower_layers[cand]
            print(f"  ✓ Auto-detected and using raw counts from layer '{matched}'")
            adata.X = adata.layers[matched].copy()
            return adata

    if adata.raw is not None:
        if set(adata.var_names).issubset(set(adata.raw.var_names)):
            adata.X = adata.raw[:, adata.var_names].X.copy()
        else:
            raw_ad = adata.raw.to_adata()
            for k, v in adata.obsm.items():
                if k not in raw_ad.obsm:
                    raw_ad.obsm[k] = v
            adata = raw_ad
        print("  ✓ Auto-detected and using raw counts from .raw")
        return adata

    partial_matches = [
        k for k in adata.layers.keys()
        if "raw" in k.lower() or "count" in k.lower()
    ]
    if partial_matches:
        matched = partial_matches[0]
        print(f"  ✓ Auto-detected and using raw counts from layer '{matched}'")
        adata.X = adata.layers[matched].copy()
        return adata

    print("  ℹ No raw counts layer or .raw detected; using .X as counts")
    return adata


def _load_scatac_data(
    scatac_input: Union[str, AnnData],
    raw_layer_atac: Optional[str] = None,
) -> AnnData:
    """
    Load single-cell ATAC-seq data from file path or AnnData object.

    Supports:
    - In-memory AnnData object
    - AnnData .h5ad file
    - MuData .h5mu file (extracting the 'atac' modality)

    Checks for and extracts raw counts from a layer or .raw, or from
    the explicitly specified `raw_layer_atac`.

    Parameters
    ----------
    scatac_input : Union[str, AnnData]
        Path to .h5ad or .h5mu file, or an in-memory AnnData object.
    raw_layer_atac : Optional[str], default None
        Name of layer in AnnData to use for raw counts (or 'raw'/'.raw').
        If None, automatically checks for layers ('raw_counts', 'counts',
        'raw') or .raw.

    Returns
    -------
    AnnData
        The loaded scATAC AnnData object with raw counts in .X.
    """
    if isinstance(scatac_input, AnnData) or hasattr(scatac_input, "var_names"):
        print("  ✓ Using in-memory scATAC AnnData object")
        adata = scatac_input.copy()
        return _extract_raw_counts(adata, raw_layer=raw_layer_atac)

    if not isinstance(scatac_input, str):
        raise TypeError(
            f"Unsupported type for scatac_input: {type(scatac_input)}. Expected str or AnnData."
        )

    if not os.path.exists(scatac_input):
        raise FileNotFoundError(f"scATAC file not found: {scatac_input}")

    if scatac_input.endswith(".h5mu"):
        try:
            import muon as mu
        except ImportError as err:
            raise ImportError(
                "Package 'muon' is required to load .h5mu files. "
                "Install via `pip install muon` or `pip install -e '.[atac]'`."
            ) from err

        print(f"  Loading MuData file: {scatac_input}")
        mdata = mu.read_h5mu(scatac_input)
        if "atac" in mdata.mod:
            print("  ✓ Extracted 'atac' modality from MuData")
            adata = mdata.mod["atac"].copy()
            return _extract_raw_counts(adata, raw_layer=raw_layer_atac)
        elif len(mdata.mod) == 1:
            mod_key = list(mdata.mod.keys())[0]
            print(f"  Note: using sole modality '{mod_key}' from {scatac_input}")
            adata = mdata.mod[mod_key].copy()
            return _extract_raw_counts(adata, raw_layer=raw_layer_atac)
        else:
            raise ValueError(
                f"Modalities in {scatac_input}: {list(mdata.mod.keys())}, "
                "but no 'atac' modality found."
            )

    # Standard AnnData (.h5ad)
    import scanpy as sc

    print(f"  Loading AnnData file: {scatac_input}")
    adata = sc.read_h5ad(scatac_input)
    return _extract_raw_counts(adata, raw_layer=raw_layer_atac)


def _parse_threshold(
    coaccess_series: pd.Series,
    threshold_spec: Optional[Union[float, int, str]] = None,
) -> float:
    """
    Parse a co-accessibility threshold specification into an absolute numerical cutoff.

    Supports:
    - Quantile string, e.g. '0.95q' or '95%'
    - Absolute score float/int, e.g. 0.8

    Parameters
    ----------
    coaccess_series : pd.Series
        Series of co-accessibility scores.
    threshold_spec : Optional[Union[float, int, str]], default None
        Threshold specification. Defaults to config.SCATAC_COACCESS_THRESHOLD.

    Returns
    -------
    float
        Calculated score cutoff.
    """
    if threshold_spec is None:
        threshold_spec = config.SCATAC_COACCESS_THRESHOLD

    if isinstance(threshold_spec, str):
        threshold_spec_clean = threshold_spec.strip().lower()
        if threshold_spec_clean.endswith("q"):
            q_val = float(threshold_spec_clean[:-1])
            return float(coaccess_series.quantile(q_val))
        elif threshold_spec_clean.endswith("%"):
            q_val = float(threshold_spec_clean[:-1]) / 100.0
            return float(coaccess_series.quantile(q_val))
        else:
            try:
                val = float(threshold_spec_clean)
                return val
            except ValueError as err:
                raise ValueError(
                    f"Invalid co-accessibility threshold specification: {threshold_spec}. "
                    "Use e.g. '0.95q' or a float score like '0.8'."
                ) from err

    return float(threshold_spec)



def _get_ref_genome(species: str) -> str:
    """
    Map species name to reference genome identifier.

    Parameters
    ----------
    species : str
        Species name (e.g., 'human', 'mouse').

    Returns
    -------
    str
        Reference genome identifier (e.g., 'hg38', 'mm10').

    Raises
    ------
    ValueError
        If species is not supported.
    """
    species_map = {
        "human": "hg38",
        "mouse": "mm10",
    }
    ref_genome = species_map.get(species.lower())
    if ref_genome is None:
        raise ValueError(
            f"Unsupported species '{species}' for ATAC peak processing. "
            f"Supported: {list(species_map.keys())}"
        )
    return ref_genome


def _ensure_genome_installed(ref_genome: str) -> None:
    """
    Check and install reference genome if not already present.

    Parameters
    ----------
    ref_genome : str
        Reference genome identifier (e.g., 'hg38').
    """
    from celloracle import motif_analysis as ma

    genome_installed = ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=None)
    print(f"  {ref_genome} installation: {genome_installed}")

    if not genome_installed:
        import genomepy

        print(f"  Installing genome {ref_genome}...")
        genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=None)
        print(f"  ✓ Genome {ref_genome} installed")
    else:
        print(f"  ✓ Genome {ref_genome} is already installed")


def _annotate_bed_peaks(bed_path: str, ref_genome: str) -> pd.DataFrame:
    """
    Read a BED file and annotate peaks with TSS / gene information.

    Uses CellOracle's motif_analysis utilities to:
    1. Read the BED file
    2. Convert peaks to string representation
    3. Annotate each peak with the nearest TSS and gene name

    Parameters
    ----------
    bed_path : str
        Path to the BED file with ATAC peaks.
    ref_genome : str
        Reference genome identifier (e.g., 'hg38').

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``['peak_id', 'gene_short_name']``.
    """
    from celloracle import motif_analysis as ma

    # Read the BED file using CellOracle
    bed = ma.read_bed(bed_path)
    print(f"  ✓ Read BED file: {len(bed)} peaks")

    # Convert to peak string list
    peaks = ma.process_bed_file.df_to_list_peakstr(bed)

    # Annotate peaks with TSS / gene information
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)
    print(f"  ✓ TSS annotation complete: {len(tss_annotated)} peaks")

    # Build final DataFrame with peak_id and gene_short_name
    peak_id_tss = ma.process_bed_file.df_to_list_peakstr(tss_annotated)
    tss_df = pd.DataFrame(
        {
            "peak_id": peak_id_tss,
            "gene_short_name": tss_annotated.gene_short_name.values,
        }
    )
    tss_df = tss_df.reset_index(drop=True)

    return tss_df


def process_atac_peaks(
    bed_path: str,
    species: str = "human",
    output_dir: Optional[str] = None,
    fpr: Optional[float] = None,
    motif_score_threshold: Optional[int] = None,
    log_dir: Optional[str] = None,
    force: bool = False,
    **kwargs,
) -> str:
    """
    Process ATAC-seq peaks BED file to generate enriched TF motif matrix.

    The full workflow is:
    1. Read the BED file and annotate peaks with TSS / gene info (CSV)
    2. Load the annotated CSV, validate peak format
    3. Scan peaks for TF binding motifs (or load existing TFinfo HDF5)
    4. Filter motifs by score
    5. Save the resulting TF info matrix DataFrame and dictionary PKL files
    6. Save checkpoint to log_dir to avoid recomputing on subsequent runs

    Parameters
    ----------
    bed_path : str
        Path to the BED file containing ATAC peaks.
    species : str, default 'human'
        Species name. Determines reference genome ('human' -> hg38,
        'mouse' -> mm10).
    output_dir : str, optional
        Directory to save output files. Defaults to config.OUTPUT_DIR.
    fpr : float, optional
        False positive rate for motif scanning.
        Defaults to config.ATAC_MOTIF_SCAN_FPR.
    motif_score_threshold : int, optional
        Minimum motif score for filtering.
        Defaults to config.ATAC_MOTIF_SCORE_THRESHOLD.
    log_dir : str, optional
        Path to log directory for checkpoint tracking.
    force : bool, default False
        If True, re-run motif scanning and processing even if checkpoint or
        output files exist.

    Returns
    -------
    str
        Path to the saved enriched ATAC peaks dictionary pickle file.

    Raises
    ------
    FileNotFoundError
        If the BED file does not exist.
    ValueError
        If species is not supported.
    """
    import pickle
    from celloracle import motif_analysis as ma

    # Resolve defaults from config
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    if fpr is None:
        fpr = config.ATAC_MOTIF_SCAN_FPR
    if motif_score_threshold is None:
        motif_score_threshold = config.ATAC_MOTIF_SCORE_THRESHOLD

    # Validate input
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"ATAC peaks BED file not found: {bed_path}")

    # Create output subdirectory
    atac_output_dir = os.path.join(output_dir, "celloracle")
    os.makedirs(atac_output_dir, exist_ok=True)

    csv_path = os.path.join(atac_output_dir, "tss_annotated_peaks.csv")
    tfi_path = os.path.join(atac_output_dir, "motif_enriched_tfi.celloracle.tfinfo")
    pkl_path_df = os.path.join(atac_output_dir, "enriched_atac_peaks_df.pkl")
    pkl_path_dict = os.path.join(atac_output_dir, "enriched_atac_peaks_dict.pkl")

    # ------------------------------------------------------------------
    # Step 0: Checkpoint check (avoid recomputing motif analysis)
    # ------------------------------------------------------------------
    step_hash = None
    if log_dir:
        from genecircuitry.pipeline.controller import (
            check_checkpoint,
            compute_input_hash,
        )

        step_hash = compute_input_hash(
            bed_path,
            species=species,
            fpr=fpr,
            threshold=motif_score_threshold,
        )

    if not force:
        has_checkpoint = bool(
            log_dir and check_checkpoint(log_dir, "atac_peaks", step_hash)
        )
        files_exist = os.path.exists(pkl_path_dict) and os.path.exists(pkl_path_df)

        if has_checkpoint or files_exist:
            print(f"\n  ✓ Found existing enriched ATAC peaks (checkpoint hit):")
            print(f"    DataFrame:  {pkl_path_df}")
            print(f"    Dictionary: {pkl_path_dict}")
            print("  ⏭ Skipping motif analysis (already computed).")
            return pkl_path_dict

    # ------------------------------------------------------------------
    # Stage 1: BED → TSS-annotated CSV
    # ------------------------------------------------------------------
    ref_genome = _get_ref_genome(species)
    _ensure_genome_installed(ref_genome)

    if not force and os.path.exists(csv_path):
        print(f"\n  [3.5.1] Found existing TSS-annotated peaks CSV: {csv_path}")
        peaks = pd.read_csv(csv_path, index_col=0)
        peaks.reset_index(drop=True, inplace=True)
        print(f"  ✓ Loaded {len(peaks)} annotated peaks from CSV")
    else:
        print(f"\n  [3.5.1] Annotating BED peaks with TSS info...")
        tss_df = _annotate_bed_peaks(bed_path, ref_genome)
        tss_df.to_csv(csv_path)
        print(f"  ✓ TSS-annotated peaks saved to: {csv_path}")
        peaks = tss_df

    # ------------------------------------------------------------------
    # Stage 2: Annotated CSV → motif scan → TFinfo HDF5
    # ------------------------------------------------------------------
    if not force and os.path.exists(tfi_path):
        print(f"\n  [3.5.2] Found existing TFinfo object: {tfi_path}")
        print("  Loading TFinfo object without re-running motif scanning...")
        tfi = ma.load_tfinfo(file_path=tfi_path)
        print("  ✓ Loaded TFinfo object")
    else:
        print(f"\n  [3.5.2] Loading annotated peaks for motif scanning...")
        peaks = pd.read_csv(csv_path, index_col=0)
        peaks.reset_index(drop=True, inplace=True)
        print(f"  ✓ Loaded {len(peaks)} annotated peaks")

        # Validate peak format
        peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

        # Create TFinfo object
        tfi = ma.TFinfo(
            peak_data_frame=peaks,
            ref_genome=ref_genome,
            genomes_dir=None,
        )

        # Load motifs
        from gimmemotifs.motif import default_motifs

        if species.lower() in ["human", "mouse"]:
            motifs = default_motifs()
        else:
            print("  WARNING - Species not recognized, using default motifs.")
            motifs = None

        # Scan for motifs
        print(f"  Scanning peaks for TF motifs (FPR={fpr})...")
        tfi.scan(fpr=fpr, motifs=motifs, verbose=False, n_cpus=config.N_JOBS)
        print("  ✓ Motif scanning complete")

        # Save TFinfo object as HDF5
        tfi.to_hdf5(file_path=tfi_path)
        print(f"  ✓ TFinfo object saved to: {tfi_path}")

    # ------------------------------------------------------------------
    # Stage 3: Filter motifs, generate DF & Dict PKL, save checkpoint
    # ------------------------------------------------------------------
    print(f"\n  [3.5.3] Filtering motifs and generating TF info matrix...")
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=motif_score_threshold)
    print(f"  ✓ Filtered motifs (score threshold={motif_score_threshold})")

    # Generate TF info dataframe and dictionary
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()
    df_dict = tfi.to_dictionary()

    # Save DataFrame result as pickle
    df.to_pickle(pkl_path_df)
    print(f"  ✓ Enriched ATAC peaks DataFrame saved to: {pkl_path_df}")
    print(f"  TF info matrix shape: {df.shape}")

    # Save Dictionary result as pickle
    with open(pkl_path_dict, "wb") as f:
        pickle.dump(df_dict, f)
    print(f"  ✓ Enriched ATAC peaks dictionary saved to: {pkl_path_dict}")

    # Save checkpoint if log_dir is provided
    if log_dir:
        from genecircuitry.pipeline.controller import (
            compute_input_hash,
            write_checkpoint,
        )

        if step_hash is None:
            step_hash = compute_input_hash(
                bed_path,
                species=species,
                fpr=fpr,
                threshold=motif_score_threshold,
            )
        write_checkpoint(
            log_dir,
            "atac_peaks",
            step_hash,
            bed_path=bed_path,
            pkl_path=pkl_path_dict,
            df_pkl_path=pkl_path_df,
            n_peaks=len(peaks),
            tf_info_shape=list(df.shape),
        )
        print(f"  ✓ Checkpoint saved to: {log_dir}/atac_peaks.checkpoint")

    return pkl_path_dict


def process_scatac_data(
    scatac_data: Union[str, AnnData],
    species: str = "human",
    output_dir: Optional[str] = None,
    raw_layer_atac: Optional[str] = None,
    coaccess_threshold: Optional[Union[float, str]] = None,
    compute_metacells: Optional[bool] = None,
    fpr: Optional[float] = None,
    motif_score_threshold: Optional[int] = None,
    log_dir: Optional[str] = None,
    force: bool = False,
    n_jobs: Optional[int] = None,
    **kwargs,
) -> str:
    """
    Process single-cell ATAC-seq data through CIRCE co-accessibility inference
    and CellOracle motif analysis to generate an enriched TF info matrix.

    Pipeline workflow:
    1. Load scATAC AnnData (.h5ad, .h5mu, or AnnData object) and extract raw counts
       from specified layer, auto-detected layer, or .raw.
    2. Standardize peak format (chr_start_end) and add CIRCE region annotations
    3. (Optional) Compute CIRCE metacells
    4. Compute ATAC co-accessibility network and extract links
    5. Annotate peaks with TSS and integrate with CIRCE connections
    6. Filter peak connections using co-accessibility threshold (default 95th percentile)
    7. Scan peaks for TF binding motifs and filter by score
    8. Save resulting TF info matrix (DF & dict PKL) and checkpoint

    Parameters
    ----------
    scatac_data : Union[str, AnnData]
        Path to .h5ad or .h5mu file, or in-memory AnnData object.
    species : str, default 'human'
        Species name ('human' -> hg38, 'mouse' -> mm10).
    output_dir : str, optional
        Directory to save output files. Defaults to config.OUTPUT_DIR.
    raw_layer_atac : str, optional
        Layer name in scATAC AnnData containing raw counts (or 'raw' / '.raw').
        If None, automatically detects raw counts from candidate layers or .raw.
    coaccess_threshold : Union[float, str], optional
        Threshold for filtering co-accessible peak connections.
        Defaults to config.SCATAC_COACCESS_THRESHOLD (e.g. '0.95q').
    compute_metacells : bool, optional
        Whether to compute CIRCE metacells before network inference.
        Defaults to config.SCATAC_COMPUTE_METACELLS.
    fpr : float, optional
        False positive rate for motif scanning. Defaults to config.ATAC_MOTIF_SCAN_FPR.
    motif_score_threshold : int, optional
        Minimum motif score for filtering. Defaults to config.ATAC_MOTIF_SCORE_THRESHOLD.
    log_dir : str, optional
        Path to log directory for checkpoint tracking.
    force : bool, default False
        If True, re-run analysis even if checkpoint or output files exist.
    n_jobs : int, optional
        Number of CPUs for parallel processing. Defaults to config.N_JOBS.

    Returns
    -------
    str
        Path to the saved enriched scATAC peaks dictionary pickle file.
    """
    # Resolve defaults from config
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    if raw_layer_atac is None:
        raw_layer_atac = kwargs.get("raw_layer", config.SCATAC_RAW_LAYER)
    if coaccess_threshold is None:
        coaccess_threshold = config.SCATAC_COACCESS_THRESHOLD
    if compute_metacells is None:
        compute_metacells = config.SCATAC_COMPUTE_METACELLS
    if fpr is None:
        fpr = config.ATAC_MOTIF_SCAN_FPR
    if motif_score_threshold is None:
        motif_score_threshold = config.ATAC_MOTIF_SCORE_THRESHOLD
    if n_jobs is None:
        n_jobs = config.N_JOBS

    # Create output subdirectory
    atac_output_dir = os.path.join(output_dir, "celloracle")
    os.makedirs(atac_output_dir, exist_ok=True)

    coaccess_csv_path = os.path.join(atac_output_dir, "scatac_coaccessibility_links.csv")
    csv_path = os.path.join(atac_output_dir, "scatac_tss_annotated_peaks.csv")
    tfi_path = os.path.join(atac_output_dir, "scatac_motif_enriched_tfi.celloracle.tfinfo")
    pkl_path_df = os.path.join(atac_output_dir, "enriched_scatac_peaks_df.pkl")
    pkl_path_dict = os.path.join(atac_output_dir, "enriched_scatac_peaks_dict.pkl")

    # ------------------------------------------------------------------
    # Step 0: Checkpoint check (avoid recomputing motif analysis)
    # ------------------------------------------------------------------
    step_hash = None
    input_repr = scatac_data if isinstance(scatac_data, str) else str(scatac_data.shape)
    if log_dir:
        from genecircuitry.pipeline.controller import (
            check_checkpoint,
            compute_input_hash,
        )

        step_hash = compute_input_hash(
            input_repr,
            species=species,
            fpr=fpr,
            threshold=motif_score_threshold,
            coaccess_threshold=str(coaccess_threshold),
            compute_metacells=compute_metacells,
            raw_layer_atac=str(raw_layer_atac),
        )

    if not force:
        has_checkpoint = bool(
            log_dir and check_checkpoint(log_dir, "scatac_data", step_hash)
        )
        files_exist = os.path.exists(pkl_path_dict) and os.path.exists(pkl_path_df)

        if has_checkpoint or files_exist:
            print(f"\n  ✓ Found existing enriched scATAC peaks (checkpoint hit):")
            print(f"    DataFrame:  {pkl_path_df}")
            print(f"    Dictionary: {pkl_path_dict}")
            print("  ⏭ Skipping scATAC CIRCE analysis (already computed).")
            return pkl_path_dict

    # ------------------------------------------------------------------
    # Stage 1: Load scATAC data & CIRCE co-accessibility inference
    # ------------------------------------------------------------------
    ref_genome = _get_ref_genome(species)
    _ensure_genome_installed(ref_genome)

    atac = _load_scatac_data(scatac_data, raw_layer_atac=raw_layer_atac)

    if not force and os.path.exists(coaccess_csv_path):
        print(f"\n  [scATAC 1/3] Found existing CIRCE co-accessibility links: {coaccess_csv_path}")
        circe_network = pd.read_csv(coaccess_csv_path)
        print(f"  ✓ Loaded {len(circe_network)} links from CSV")
    else:
        try:
            import circe as ci
        except ImportError as err:
            raise ImportError(
                "Package 'circe' (circe-py) is required for scATAC co-accessibility analysis. "
                "Install via `pip install circe-py` or `pip install -e '.[atac]'`."
            ) from err

        print("\n  [scATAC 1/3] Standardizing peak identifiers and adding region annotations...")
        atac.var_names = (
            atac.var_names.str.replace(":", "_", regex=False)
            .str.replace("-", "_", regex=False)
        )
        atac = ci.add_region_infos(atac, sep=("_", "_"))

        if compute_metacells:
            print("  Computing CIRCE metacells...")
            atac = ci.metacells.compute_metacells(atac)
        else:
            print("  Skipping metacell aggregation (using single cells)")

        print("  Inferring cis-coaccessibility network with CIRCE...")
        ci.compute_atac_network(atac,n_jobs=config.GRN_N_JOBS)
        circe_network = ci.extract_atac_links(atac)
        circe_network = circe_network.rename(columns={"score": "coaccess"})
        circe_network.to_csv(coaccess_csv_path, index=False)
        print(
            f"  ✓ CIRCE co-accessibility links saved ({len(circe_network)} links) to: {coaccess_csv_path}"
        )

    # ------------------------------------------------------------------
    # Stage 2: TSS annotation & CIRCE connections integration
    # ------------------------------------------------------------------
    if not force and os.path.exists(csv_path):
        print(f"\n  [scATAC 2/3] Found existing TSS-annotated peaks CSV: {csv_path}")
        peaks = pd.read_csv(csv_path, index_col=0)
        peaks.reset_index(drop=True, inplace=True)
        print(f"  ✓ Loaded {len(peaks)} annotated peaks from CSV")
    else:
        from celloracle import motif_analysis as ma

        print("\n  [scATAC 2/3] Annotating peaks with TSS info...")
        peak_list = (
            atac.var_names.values
            if hasattr(atac, "var_names")
            else list(set(circe_network["Peak1"]).union(set(circe_network["Peak2"])))
        )
        tss_annotated = ma.get_tss_info(peak_str_list=peak_list, ref_genome=ref_genome)
        print(f"  ✓ Found {len(tss_annotated)} TSS annotations")

        print("  Integrating TSS peaks with CIRCE co-accessibility network...")
        integrated = ma.integrate_tss_peak_with_cicero(
            tss_peak=tss_annotated,
            cicero_connections=circe_network,
        )

        cutoff = _parse_threshold(integrated["coaccess"], coaccess_threshold)
        peaks = integrated[integrated["coaccess"] >= cutoff]
        print(
            f"  ✓ Filtered peak connections with coaccess >= {cutoff:.4f}: {len(peaks)} links kept"
        )
        peaks = peaks[["peak_id", "gene_short_name"]].drop_duplicates().reset_index(drop=True)
        peaks.to_csv(csv_path)
        print(f"  ✓ Saved TSS-annotated filtered peaks to: {csv_path}")

    # ------------------------------------------------------------------
    # Stage 3: TF motif scanning & filtering
    # ------------------------------------------------------------------
    from celloracle import motif_analysis as ma

    if not force and os.path.exists(tfi_path):
        print(f"\n  [scATAC 3/3] Found existing TFinfo object: {tfi_path}")
        print("  Loading TFinfo object without re-running motif scanning...")
        tfi = ma.load_tfinfo(file_path=tfi_path)
        print("  ✓ Loaded TFinfo object")
    else:
        print(f"\n  [scATAC 3/3] Scanning peaks for TF motifs...")
        peaks = pd.read_csv(csv_path, index_col=0)
        peaks.reset_index(drop=True, inplace=True)
        peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

        tfi = ma.TFinfo(
            peak_data_frame=peaks,
            ref_genome=ref_genome,
            genomes_dir=None,
        )

        from gimmemotifs.motif import default_motifs

        if species.lower() in ["human", "mouse"]:
            motifs = default_motifs()
        else:
            print("  WARNING - Species not recognized, using default motifs.")
            motifs = None

        print(f"  Scanning peaks for TF motifs (FPR={fpr}, n_jobs={n_jobs})...")
        tfi.scan(fpr=fpr, motifs=motifs, verbose=False, n_cpus=n_jobs)
        tfi.to_hdf5(file_path=tfi_path)
        print(f"  ✓ TFinfo object saved to: {tfi_path}")

    # Post-filtering
    print(f"\n  Filtering motifs and generating TF info matrix...")
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=motif_score_threshold)
    print(f"  ✓ Filtered motifs (score threshold={motif_score_threshold})")

    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()
    df_dict = tfi.to_dictionary()

    df.to_pickle(pkl_path_df)
    print(f"  ✓ Enriched scATAC peaks DataFrame saved to: {pkl_path_df}")
    print(f"  TF info matrix shape: {df.shape}")

    with open(pkl_path_dict, "wb") as f:
        pickle.dump(df_dict, f)
    print(f"  ✓ Enriched scATAC peaks dictionary saved to: {pkl_path_dict}")

    # Checkpointing
    if log_dir:
        from genecircuitry.pipeline.controller import write_checkpoint

        if step_hash is None:
            step_hash = compute_input_hash(
                input_repr,
                species=species,
                fpr=fpr,
                threshold=motif_score_threshold,
                coaccess_threshold=str(coaccess_threshold),
                compute_metacells=compute_metacells,
            )
        write_checkpoint(
            log_dir,
            "scatac_data",
            step_hash,
            scatac_data=str(scatac_data) if isinstance(scatac_data, str) else "in_memory_adata",
            pkl_path=pkl_path_dict,
            df_pkl_path=pkl_path_df,
            n_peaks=len(peaks),
            tf_info_shape=list(df.shape),
        )
        print(f"  ✓ Checkpoint saved to: {log_dir}/scatac_data.checkpoint")

    return pkl_path_dict

