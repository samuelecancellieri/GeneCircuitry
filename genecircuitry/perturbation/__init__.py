"""
In-Silico Perturbation Module for GeneCircuitry
===============================================

Provides end-to-end in-silico gene perturbation analysis based on CellOracle GRN inference:
- Automatic selection of top driver TFs via differential centrality and Hotspot modules
- Directional signal propagation simulation (KO / OE)
- Cell identity shift vector fields and grid arrows
- Trajectory flow and inner product perturbation scores (PS)
- Cluster-to-cluster transition analysis and displacement metrics
- Publication-ready comparative figures
"""

from .tf_selection import (
    select_top_perturbed_tfs,
    compute_cluster_pair_centrality_diff,
    compute_tf_hotspot_enrichment,
)
from .simulation import (
    setup_oracle_for_simulation,
    simulate_single_tf_perturbation,
    extract_target_gene_shifts,
)
from .development_analysis import (
    compute_cluster_displacement_score,
    compute_perturbation_cluster_transitions,
    compute_developmental_flow_inner_product,
)
from .pipeline import run_perturbation_analysis

__all__ = [
    "select_top_perturbed_tfs",
    "compute_cluster_pair_centrality_diff",
    "compute_tf_hotspot_enrichment",
    "setup_oracle_for_simulation",
    "simulate_single_tf_perturbation",
    "extract_target_gene_shifts",
    "compute_cluster_displacement_score",
    "compute_perturbation_cluster_transitions",
    "compute_developmental_flow_inner_product",
    "run_perturbation_analysis",
]
