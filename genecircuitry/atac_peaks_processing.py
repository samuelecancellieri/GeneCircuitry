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
from typing import Optional

import pandas as pd

from genecircuitry import config


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
) -> str:
    """
    Process ATAC-seq peaks BED file to generate enriched TF motif matrix.

    The full workflow is:
    1. Read the BED file and annotate peaks with TSS / gene info (CSV)
    2. Load the annotated CSV, validate peak format
    3. Scan peaks for TF binding motifs
    4. Filter motifs by score
    5. Save the resulting TF info matrix as a pickle file

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

    Returns
    -------
    str
        Path to the saved enriched ATAC peaks pickle file.

    Raises
    ------
    FileNotFoundError
        If the BED file does not exist.
    ValueError
        If species is not supported.
    """
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

    # ------------------------------------------------------------------
    # Stage 1: BED → TSS-annotated CSV
    # ------------------------------------------------------------------
    ref_genome = _get_ref_genome(species)
    _ensure_genome_installed(ref_genome)

    print(f"\n  [3.5.1] Annotating BED peaks with TSS info...")
    tss_df = _annotate_bed_peaks(bed_path, ref_genome)

    # Save annotated peaks as CSV
    csv_path = os.path.join(atac_output_dir, "tss_annotated_peaks.csv")
    tss_df.to_csv(csv_path)
    print(f"  ✓ TSS-annotated peaks saved to: {csv_path}")

    # ------------------------------------------------------------------
    # Stage 2: Annotated CSV → motif scan → filtered PKL
    # ------------------------------------------------------------------
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
    tfi_path = os.path.join(atac_output_dir, "motif_enriched_tfi.celloracle.tfinfo")
    tfi.to_hdf5(file_path=tfi_path)
    print(f"  ✓ TFinfo object saved to: {tfi_path}")

    # Filter motifs
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=motif_score_threshold)
    print(f"  ✓ Filtered motifs (score threshold={motif_score_threshold})")

    # Generate TF info dataframe
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()

    # Save result as pickle
    pkl_path = os.path.join(atac_output_dir, "enriched_atac_peaks.pkl")
    df.to_pickle(pkl_path)
    print(f"  ✓ Enriched ATAC peaks saved to: {pkl_path}")
    print(f"  TF info matrix shape: {df.shape}")

    return pkl_path
