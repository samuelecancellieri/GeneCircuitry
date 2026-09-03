"""
Unit Tests for Local Enrichment Analysis Engine
===============================================

Tests for genecircuitry.enrichment_analysis module:
- Bundled gene set loading
- Local ORA hypergeometric testing
- Custom dictionary and GMT file support
- Background universe options
- Offline execution (zero network requests)
- Integration with comparative analysis and module annotations
"""

import os
import socket
import tempfile
import shutil
from pathlib import Path
import pytest
import pandas as pd
import numpy as np

from genecircuitry import config
from genecircuitry.enrichment_analysis import (
    EnrichmentResult,
    get_bundled_gene_sets_dir,
    get_cache_gene_sets_dir,
    read_gmt,
    write_gmt,
    list_available_gene_sets,
    load_gene_set,
    cache_gene_sets,
    perform_ora_enrichment,
    gseapy_ora_enrichment_analysis,
    run_enrichr_online,
)
from genecircuitry.comparative_analysis import compute_module_pathway_enrichments
from genecircuitry.pipeline.controller import create_parser


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    d = tempfile.mkdtemp(prefix="test_enrichment_")
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def sample_gene_set_dict():
    """Create a sample dictionary of gene sets."""
    return {
        "Pathway_Apoptosis": [
            "TP53", "BAX", "CASP3", "CASP8", "CASP9", "FAS", "BCL2", "BID", "APAF1", "CYCS"
        ],
        "Pathway_Cell_Cycle": [
            "CDK1", "CDK2", "CDK4", "CDK6", "CCNA1", "CCNB1", "CCND1", "CCNE1", "E2F1", "RB1"
        ],
        "Pathway_DNA_Repair": [
            "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "RAD51", "PARP1", "TP53BP1"
        ],
    }


@pytest.fixture
def sample_gmt_file(temp_dir, sample_gene_set_dict):
    """Create a temporary GMT file."""
    gmt_path = os.path.join(temp_dir, "custom_pathways.gmt")
    write_gmt(sample_gene_set_dict, gmt_path, description="Custom test pathways")
    return gmt_path


# =============================================================================
# Bundled & GMT File I/O Tests
# =============================================================================


def test_bundled_gene_sets_exist():
    """Verify that bundled gene sets exist in package data."""
    bundled_dir = get_bundled_gene_sets_dir()
    assert os.path.isdir(bundled_dir), f"Bundled dir not found: {bundled_dir}"

    hallmark_gmt = os.path.join(bundled_dir, "MSigDB_Hallmark_2020.gmt")
    assert os.path.isfile(hallmark_gmt), "MSigDB_Hallmark_2020.gmt missing"

    available = list_available_gene_sets(include_cached=False)
    assert "MSigDB_Hallmark_2020" in available


def test_read_and_write_gmt(temp_dir, sample_gene_set_dict):
    """Test reading and writing GMT files."""
    dest = os.path.join(temp_dir, "test.gmt")
    write_gmt(sample_gene_set_dict, dest, description="Test description")

    assert os.path.isfile(dest)
    loaded = read_gmt(dest)

    assert len(loaded) == len(sample_gene_set_dict)
    for term, genes in sample_gene_set_dict.items():
        assert term in loaded
        assert set(loaded[term]) == set(genes)


def test_read_gmt_nonexistent():
    """Test reading a non-existent GMT file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_gmt("/non/existent/path.gmt")


# =============================================================================
# Gene Set Resolver Tests
# =============================================================================


def test_load_gene_set_dict(sample_gene_set_dict):
    """Test loading a gene set from a Python dictionary."""
    name, gs = load_gene_set(sample_gene_set_dict)
    assert name == "custom"
    assert gs == sample_gene_set_dict


def test_load_gene_set_file(sample_gmt_file, sample_gene_set_dict):
    """Test loading a gene set from a GMT file path."""
    name, gs = load_gene_set(sample_gmt_file)
    assert name == "custom_pathways"
    assert len(gs) == len(sample_gene_set_dict)


def test_load_gene_set_bundled():
    """Test loading default bundled hallmark gene sets."""
    name, gs = load_gene_set("MSigDB_Hallmark_2020", organism="human")
    assert name == "MSigDB_Hallmark_2020"
    assert isinstance(gs, dict)
    assert len(gs) == 50
    assert "Apoptosis" in gs or "HALLMARK_APOPTOSIS" in gs or any("Apoptosis" in k for k in gs.keys())


def test_load_gene_set_not_found(temp_dir):
    """Test requesting a non-existent gene set raises ValueError."""
    # Use empty cache dir so it doesn't find external files
    with pytest.raises(ValueError) as excinfo:
        load_gene_set("NON_EXISTENT_LIBRARY_XYZ", cache_dir=temp_dir)
    assert "could not be found locally" in str(excinfo.value)


# =============================================================================
# Local ORA Statistical Testing
# =============================================================================


def test_perform_ora_enrichment_bundled():
    """Test local ORA enrichment with bundled MSigDB Hallmark."""
    # Query genes for apoptosis and p53
    query = ["TP53", "MDM2", "BAX", "CASP3", "CASP8", "FAS", "ATM", "CHEK2", "CDKN1A"]
    res = perform_ora_enrichment(query, gene_sets=["MSigDB_Hallmark_2020"], pval_cutoff=0.05)

    assert isinstance(res, EnrichmentResult)
    assert not res.empty
    assert isinstance(res.results, pd.DataFrame)
    assert "Term" in res.results.columns
    assert "Adjusted P-value" in res.results.columns
    assert "Combined Score" in res.results.columns
    assert "Overlap" in res.results.columns
    assert "Genes" in res.results.columns

    # Check that top hits include Apoptosis or p53 Pathway
    terms = res.results["Term"].tolist()
    assert any("Apoptosis" in t or "p53" in t for t in terms)

    # Verify P-values are significant and sorted
    assert res.results["Adjusted P-value"].iloc[0] <= 0.05
    assert res.results["P-value"].is_monotonic_increasing


def test_perform_ora_enrichment_custom_dict(sample_gene_set_dict):
    """Test local ORA enrichment with custom dictionary."""
    query = ["TP53", "BAX", "CASP3", "CASP8", "CASP9"]
    res = perform_ora_enrichment(query, gene_sets=sample_gene_set_dict, pval_cutoff=0.05)

    assert not res.empty
    top_term = res.results.iloc[0]["Term"]
    assert top_term == "Pathway_Apoptosis"
    assert res.results.iloc[0]["Adjusted P-value"] < 0.05


def test_perform_ora_background_options(sample_gene_set_dict):
    """Test ORA with different background universe specifications."""
    query = ["TP53", "BAX", "CASP3", "CASP8", "CASP9"]

    # 1. Background = None (library union)
    res_none = perform_ora_enrichment(query, gene_sets=sample_gene_set_dict, background=None, pval_cutoff=0.1)
    # 2. Background = int
    res_int = perform_ora_enrichment(query, gene_sets=sample_gene_set_dict, background=20000, pval_cutoff=0.1)
    # 3. Background = list of gene symbols
    custom_bg = list(set.union(*[set(v) for v in sample_gene_set_dict.values()])) + [f"OTHER_{i}" for i in range(1000)]
    res_list = perform_ora_enrichment(query, gene_sets=sample_gene_set_dict, background=custom_bg, pval_cutoff=0.1)

    assert not res_none.empty
    assert not res_int.empty
    assert not res_list.empty

    # With larger genome background, P-value is smaller
    p_none = res_none.results.iloc[0]["P-value"]
    p_int = res_int.results.iloc[0]["P-value"]
    assert p_int < p_none


def test_perform_ora_empty_and_no_hits(sample_gene_set_dict):
    """Test edge cases: empty gene list, non-overlapping gene list."""
    # Empty query
    res_empty = perform_ora_enrichment([], gene_sets=sample_gene_set_dict)
    assert res_empty.empty
    assert list(res_empty.results.columns) == [
        "Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"
    ]

    # No overlap
    res_no_overlap = perform_ora_enrichment(["UNKNOWN1", "UNKNOWN2"], gene_sets=sample_gene_set_dict)
    assert res_no_overlap.empty
    assert isinstance(res_no_overlap.results, pd.DataFrame)


def test_gseapy_ora_enrichment_analysis_compatibility(sample_gene_set_dict):
    """Test that legacy gseapy_ora_enrichment_analysis function behaves identically."""
    query = ["TP53", "BAX", "CASP3", "CASP8", "CASP9"]
    res = gseapy_ora_enrichment_analysis(query, gene_sets=sample_gene_set_dict, pval_cutoff=0.05)

    assert isinstance(res, EnrichmentResult)
    assert hasattr(res, "results")
    assert hasattr(res, "res2d")
    assert not res.results.empty
    assert res.results.iloc[0]["Term"] == "Pathway_Apoptosis"


# =============================================================================
# Offline Verification (Network Blocking)
# =============================================================================


def test_offline_zero_network_calls(monkeypatch, sample_gene_set_dict):
    """Verify that ORA execution makes zero socket connections."""
    def block_socket(*args, **kwargs):
        raise RuntimeError("Network socket connection attempted during local enrichment!")

    monkeypatch.setattr(socket, "socket", block_socket)

    query = ["TP53", "MDM2", "BAX", "CASP3", "CASP8"]

    # 1. Bundled library
    res_bundled = perform_ora_enrichment(query, gene_sets=["MSigDB_Hallmark_2020"])
    assert not res_bundled.empty

    # 2. Custom dict
    res_dict = perform_ora_enrichment(["TP53", "BAX", "CASP3", "CASP8", "CASP9"], gene_sets=sample_gene_set_dict)
    assert not res_dict.empty


# =============================================================================
# Integration Tests
# =============================================================================


def test_integration_comparative_pathway_enrichments(sample_gene_set_dict):
    """Test compute_module_pathway_enrichments with local enrichment."""
    module_dict = {
        "1": ["TP53", "BAX", "CASP3", "CASP8", "CASP9"],
        "2": ["CDK1", "CDK2", "CDK4", "CCNA1", "CCNB1"],
    }

    df_enr = compute_module_pathway_enrichments(
        module_dict,
        gene_sets=sample_gene_set_dict,
        top_n_terms=2,
        pval_cutoff=0.05,
    )

    assert isinstance(df_enr, pd.DataFrame)
    assert not df_enr.empty
    assert "module" in df_enr.columns
    assert "term" in df_enr.columns
    assert "adjusted_p_value" in df_enr.columns

    # Module 1 should be enriched for Apoptosis
    m1_rows = df_enr[df_enr["module"] == "Module 1"]
    assert not m1_rows.empty
    assert "Apoptosis" in m1_rows["term"].values[0]


# =============================================================================
# Online Enrichr API Tests
# =============================================================================


class DummyEnrichrObj:
    """Mock object returned by gseapy.enrichr."""

    def __init__(self, df: pd.DataFrame):
        self.results = df
        self.res2d = df


def test_perform_ora_enrichment_online_mocked(monkeypatch):
    """Test online ORA enrichment workflow with mocked gseapy.enrichr."""
    mock_results = pd.DataFrame(
        [
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "Apoptosis",
                "Overlap": "4/161",
                "P-value": 0.0001,
                "Adjusted P-value": 0.001,
                "Odds Ratio": 50.0,
                "Combined Score": 350.0,
                "Genes": "CASP3;CASP8;BAX;TP53",
            },
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "p53 Pathway",
                "Overlap": "3/200",
                "P-value": 0.0005,
                "Adjusted P-value": 0.002,
                "Odds Ratio": 30.0,
                "Combined Score": 200.0,
                "Genes": "TP53;MDM2;BAX",
            },
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "Insignificant Pathway",
                "Overlap": "1/200",
                "P-value": 0.2,
                "Adjusted P-value": 0.4,
                "Odds Ratio": 1.5,
                "Combined Score": 2.0,
                "Genes": "BAX",
            },
        ]
    )

    calls = []

    def mock_enrichr(gene_list, gene_sets, organism, **kwargs):
        calls.append(
            {"gene_list": gene_list, "gene_sets": gene_sets, "organism": organism}
        )
        return DummyEnrichrObj(mock_results)

    import gseapy

    monkeypatch.setattr(gseapy, "enrichr", mock_enrichr)

    query = ["TP53", "MDM2", "BAX", "CASP3", "CASP8"]
    res = perform_ora_enrichment(
        query,
        gene_sets=["MSigDB_Hallmark_2020"],
        pval_cutoff=0.05,
        online=True,
    )

    assert len(calls) == 1
    assert calls[0]["organism"] == "human"
    assert isinstance(res, EnrichmentResult)
    assert not res.empty
    # Insignificant pathway should be filtered out by pval_cutoff=0.05
    assert len(res.results) == 2
    assert "Apoptosis" in res.results["Term"].values
    assert "p53 Pathway" in res.results["Term"].values
    assert "Insignificant Pathway" not in res.results["Term"].values


def test_perform_ora_enrichment_online_via_config(monkeypatch):
    """Test that setting config.ENRICHMENT_ONLINE=True activates online Enrichr by default."""
    mock_results = pd.DataFrame(
        [
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "Apoptosis",
                "Overlap": "3/161",
                "P-value": 0.001,
                "Adjusted P-value": 0.01,
                "Odds Ratio": 40.0,
                "Combined Score": 180.0,
                "Genes": "CASP3;BAX;TP53",
            }
        ]
    )

    calls = []

    def mock_enrichr(gene_list, gene_sets, organism, **kwargs):
        calls.append(True)
        return DummyEnrichrObj(mock_results)

    import gseapy

    monkeypatch.setattr(gseapy, "enrichr", mock_enrichr)

    # Save original setting
    original = config.ENRICHMENT_ONLINE
    try:
        config.ENRICHMENT_ONLINE = True
        query = ["TP53", "BAX", "CASP3"]
        res = perform_ora_enrichment(query)
        assert len(calls) == 1
        assert not res.empty
        assert res.results.iloc[0]["Term"] == "Apoptosis"

        # Also test legacy wrapper gseapy_ora_enrichment_analysis
        res2 = gseapy_ora_enrichment_analysis(query)
        assert len(calls) == 2
        assert not res2.empty
    finally:
        config.ENRICHMENT_ONLINE = original


def test_online_enrichr_retry_on_429(monkeypatch):
    """Test that _run_enrichr_online retries on HTTP 429 and succeeds."""
    mock_results = pd.DataFrame(
        [
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "Apoptosis",
                "Overlap": "3/161",
                "P-value": 0.001,
                "Adjusted P-value": 0.01,
                "Odds Ratio": 40.0,
                "Combined Score": 180.0,
                "Genes": "CASP3;BAX;TP53",
            }
        ]
    )

    call_count = [0]

    def mock_enrichr_rate_limited(*args, **kwargs):
        call_count[0] += 1
        if call_count[0] == 1:
            raise RuntimeError("HTTP 429: Too Many Requests")
        return DummyEnrichrObj(mock_results)

    import gseapy

    monkeypatch.setattr(gseapy, "enrichr", mock_enrichr_rate_limited)

    query = ["TP53", "BAX", "CASP3"]
    res = run_enrichr_online(
        query,
        gene_sets=["MSigDB_Hallmark_2020"],
        max_retries=3,
        retry_delay=0.01,
    )

    assert call_count[0] == 2
    assert not res.empty
    assert res.results.iloc[0]["Term"] == "Apoptosis"


def test_online_enrichr_retry_exhausted(monkeypatch):
    """Test that _run_enrichr_online raises RuntimeError after exhausting retries."""

    def mock_enrichr_fail(*args, **kwargs):
        raise RuntimeError("HTTP 429: Rate limit exceeded permanently")

    import gseapy

    monkeypatch.setattr(gseapy, "enrichr", mock_enrichr_fail)

    with pytest.raises(RuntimeError) as excinfo:
        run_enrichr_online(
            ["TP53", "BAX"],
            gene_sets=["MSigDB_Hallmark_2020"],
            max_retries=2,
            retry_delay=0.01,
        )

    assert "Online Enrichr API analysis failed" in str(excinfo.value)
    assert "config.ENRICHMENT_ONLINE = False" in str(excinfo.value)


def test_online_enrichr_cli_flags():
    """Test that controller argument parser recognizes --enrichment-online and its aliases."""
    parser = create_parser()

    args1 = parser.parse_args(["--enrichment-online"])
    assert args1.enrichment_online is True

    args2 = parser.parse_args(["--use-enrichr"])
    assert args2.enrichment_online is True

    args3 = parser.parse_args(["--online-enrichment"])
    assert args3.enrichment_online is True

    args_default = parser.parse_args([])
    assert args_default.enrichment_online is False


def test_online_enrichr_comparative_integration(monkeypatch):
    """Test compute_module_pathway_enrichments with online=True."""
    mock_results = pd.DataFrame(
        [
            {
                "Gene_set": "MSigDB_Hallmark_2020",
                "Term": "Apoptosis",
                "Overlap": "3/161",
                "P-value": 0.001,
                "Adjusted P-value": 0.01,
                "Odds Ratio": 40.0,
                "Combined Score": 180.0,
                "Genes": "CASP3;BAX;TP53",
            }
        ]
    )

    calls = []

    def mock_enrichr(*args, **kwargs):
        calls.append(True)
        return DummyEnrichrObj(mock_results)

    import gseapy

    monkeypatch.setattr(gseapy, "enrichr", mock_enrichr)

    module_dict = {"1": ["TP53", "BAX", "CASP3", "CASP8", "CASP9"]}
    df_enr = compute_module_pathway_enrichments(
        module_dict,
        gene_sets=["MSigDB_Hallmark_2020"],
        online=True,
    )

    assert len(calls) == 1
    assert not df_enr.empty
    assert "Apoptosis" in df_enr["term"].values[0]


def test_online_enrichr_live_call():
    """Live test of online Enrichr API (skipped gracefully if no network connection)."""
    import urllib.request

    try:
        urllib.request.urlopen("https://maayanlab.cloud/Enrichr/", timeout=3)
    except Exception:
        pytest.skip("Online Enrichr API is unreachable, skipping live test.")

    res = run_enrichr_online(
        ["TP53", "MDM2", "BAX", "CASP3", "CASP8", "FAS", "ATM"],
        gene_sets=["MSigDB_Hallmark_2020"],
        pval_cutoff=0.05,
    )
    assert isinstance(res, EnrichmentResult)
    assert not res.empty
    assert "Term" in res.results.columns
    assert any("Apoptosis" in t or "p53" in t for t in res.results["Term"].values)

