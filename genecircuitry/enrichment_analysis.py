import gseapy

# from gseapy import barplot

from . import config


def gseapy_ora_enrichment_analysis(
    gene_list: list,
    gene_sets: list = config.ENRICHMENT_GENE_SETS,
    pval_cutoff: float = 0.05,
    species: str = "human",
):
    """
    Perform ORA enrichment analysis using gseapy.

    Parameters:
    gene_list (list): A list of gene symbols ranked by some metric (e.g., log fold change).
    gene_sets (str or dict): Path to a gene set file in GMT format or a dictionary of gene sets.
    outdir (str): Directory to save the results.

    Returns:
    gseapy.enrichr: The result object containing enrichment analysis results.
    """
    # Perform GSEA prerank analysis

    colors = dict()
    for i, gene_set in enumerate(gene_sets):
        colors[gene_set] = f"C{i}"

    import time

    sleep_time = 2
    max_sleep = 5
    enr = None

    while sleep_time <= max_sleep:
        time.sleep(sleep_time)
        try:
            enr = gseapy.enrichr(
                gene_list=gene_list[:100], gene_sets=gene_sets, organism=species
            )
            break
        except Exception as e:
            if "429" in str(e):
                if sleep_time == max_sleep:
                    raise e
                sleep_time += 1
            else:
                raise e
    enr.results = enr.results[enr.results["Adjusted P-value"] < pval_cutoff]

    return enr
