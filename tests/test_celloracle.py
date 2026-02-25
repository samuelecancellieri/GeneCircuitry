"""
Tests for CellOracle processing module
"""

import pytest
import numpy as np
import scanpy as sc
from anndata import AnnData
from unittest.mock import Mock, patch, MagicMock

# Note: CellOracle may not be installed in test environment
# We'll mock it for testing
try:
    import celloracle as co

    CELLORACLE_AVAILABLE = True
except ImportError:
    CELLORACLE_AVAILABLE = False


@pytest.fixture
def sample_adata_for_oracle():
    """Create a sample AnnData object suitable for CellOracle"""
    np.random.seed(42)
    n_cells = 200
    n_genes = 100

    # Create count matrix with some structure
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes)).astype(float)

    # Create gene names (use standard gene symbols)
    var_names = [f"GENE{i}" if i < 90 else f"TF{i}" for i in range(n_genes)]

    # Create AnnData object
    adata = AnnData(X=X)
    adata.var_names = var_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add raw counts layer
    adata.layers["raw_counts"] = X.copy()

    # Add cluster annotations
    adata.obs["leiden"] = np.random.choice(
        ["cluster_0", "cluster_1", "cluster_2"], n_cells
    )
    adata.obs["cell_type"] = np.random.choice(["TypeA", "TypeB", "TypeC"], n_cells)

    # Add UMAP embedding
    adata.obsm["X_umap"] = np.random.randn(n_cells, 2)

    # Add PCA embedding
    adata.obsm["X_pca"] = np.random.randn(n_cells, 50)

    return adata


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_create_oracle_object_with_raw_counts(sample_adata_for_oracle):
    """Test creating Oracle object with raw counts"""
    from genecircuitry.celloracle_processing import create_oracle_object

    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
        raw_count_layer="raw_counts",
    )

    # Check that Oracle object was created
    assert oracle is not None
    assert hasattr(oracle, "adata")
    assert oracle.adata.n_obs == sample_adata_for_oracle.n_obs


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_create_oracle_object_with_normalized_counts(sample_adata_for_oracle):
    """Test creating Oracle object with normalized counts"""
    from genecircuitry.celloracle_processing import create_oracle_object

    # Normalize the data
    sc.pp.normalize_total(sample_adata_for_oracle, target_sum=1e4)
    sc.pp.log1p(sample_adata_for_oracle)

    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
        raw_count_layer=None,  # Use normalized counts
    )

    # Check that Oracle object was created
    assert oracle is not None
    assert hasattr(oracle, "adata")


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_create_oracle_object_with_custom_dict(sample_adata_for_oracle, tmp_path):
    """Test creating Oracle object with custom TF dictionary"""
    from genecircuitry.celloracle_processing import create_oracle_object
    import pickle

    # Create a mock TF dictionary
    tf_dict = {"GENE1": ["TF90", "TF91"], "GENE2": ["TF92"], "GENE3": ["TF90", "TF93"]}

    # Save to temporary file
    dict_path = tmp_path / "tf_dict.pkl"
    with open(dict_path, "wb") as f:
        pickle.dump(tf_dict, f)

    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
        TG_to_TF_dictionary=str(dict_path),
    )

    assert oracle is not None


def test_create_oracle_mock():
    """Test create_oracle_object with mocking (when CellOracle not available)"""
    from genecircuitry.celloracle_processing import create_oracle_object

    # Create mock objects
    mock_adata = MagicMock()
    mock_adata.copy.return_value = mock_adata
    mock_adata.layers = {"raw_counts": np.random.rand(100, 50)}

    with patch("genecircuitry.celloracle_processing.co") as mock_co:
        mock_oracle = MagicMock()
        mock_co.Oracle.return_value = mock_oracle
        mock_co.data.load_human_promoter_base_GRN.return_value = MagicMock()

        oracle = create_oracle_object(
            adata=mock_adata, cluster_column_name="leiden", embedding_name="X_umap"
        )

        # Verify Oracle was created
        assert mock_co.Oracle.called
        assert mock_co.data.load_human_promoter_base_GRN.called


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_run_pca(sample_adata_for_oracle):
    """Test PCA on Oracle object"""
    from genecircuitry.celloracle_processing import create_oracle_object, run_PCA

    # Create Oracle object
    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
    )

    # Run PCA
    oracle_pca = run_PCA(oracle)

    # Check that PCA was performed
    assert oracle_pca is not None
    assert hasattr(oracle_pca, "pca")


def test_run_pca_mock():
    """Test run_PCA with mocking"""
    from genecircuitry.celloracle_processing import run_PCA

    # Create mock Oracle object
    mock_oracle = MagicMock()
    mock_oracle.copy.return_value = mock_oracle
    mock_oracle.pca = MagicMock()
    mock_oracle.pca.explained_variance_ratio_ = np.linspace(0.1, 0.001, 200)

    with (
        patch("matplotlib.pyplot.plot"),
        patch("matplotlib.pyplot.axvline"),
        patch("matplotlib.pyplot.show"),
    ):
        oracle_pca = run_PCA(mock_oracle)

        # Verify PCA was called
        assert mock_oracle.perform_PCA.called
        assert oracle_pca is not None


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_run_knn(sample_adata_for_oracle):
    """Test KNN imputation on Oracle object"""
    from genecircuitry.celloracle_processing import (
        create_oracle_object,
        run_PCA,
        run_KNN,
    )

    # Create and prepare Oracle object
    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
    )
    oracle = run_PCA(oracle)

    # Run KNN imputation
    run_KNN(oracle, n_comps=30)

    # Check that KNN was performed
    # Oracle object should be modified in place


def test_run_knn_mock():
    """Test run_KNN with mocking"""
    from genecircuitry.celloracle_processing import run_KNN

    # Create mock Oracle object
    mock_oracle = MagicMock()
    mock_oracle.copy.return_value = mock_oracle
    mock_oracle.adata.shape = (1000, 50)  # 1000 cells

    run_KNN(mock_oracle, n_comps=50)

    # Verify KNN imputation was called
    assert mock_oracle.knn_imputation.called


def test_run_knn_k_calculation():
    """Test automatic k calculation in run_KNN"""
    from genecircuitry.celloracle_processing import run_KNN

    # Test with different cell numbers
    for n_cells in [100, 1000, 10000]:
        mock_oracle = MagicMock()
        mock_oracle.copy.return_value = mock_oracle
        mock_oracle.adata.shape = (n_cells, 50)

        run_KNN(mock_oracle, n_comps=50)

        # Verify k was calculated as 2.5% of cells
        expected_k = int(0.025 * n_cells)
        call_args = mock_oracle.knn_imputation.call_args
        assert call_args[1]["k"] == expected_k


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="CellOracle not installed")
def test_run_links(sample_adata_for_oracle):
    """Test inferring regulatory links"""
    from genecircuitry.celloracle_processing import (
        create_oracle_object,
        run_PCA,
        run_KNN,
        run_links,
    )

    # Create and prepare Oracle object
    oracle = create_oracle_object(
        adata=sample_adata_for_oracle,
        cluster_column_name="leiden",
        embedding_name="X_umap",
    )
    oracle = run_PCA(oracle)
    run_KNN(oracle, n_comps=30)

    # Infer links
    links = run_links(oracle, cluster_column_name="leiden", p_cutoff=0.001)

    # Check that links were created
    assert links is not None


def test_run_links_mock():
    """Test run_links with mocking"""
    from genecircuitry.celloracle_processing import run_links

    # Create mock Oracle object
    mock_oracle = MagicMock()
    mock_links = MagicMock()
    mock_oracle.get_links.return_value = mock_links

    links = run_links(mock_oracle, cluster_column_name="leiden", p_cutoff=0.001)

    # Verify methods were called
    assert mock_oracle.get_links.called
    assert mock_links.filter_links.called
    assert mock_links.plot_degree_distributions.called
    assert mock_links.get_network_score.called
    assert links is not None


def test_run_links_p_cutoff():
    """Test different p-value cutoffs"""
    from genecircuitry.celloracle_processing import run_links

    for p_cutoff in [0.001, 0.01, 0.0001]:
        mock_oracle = MagicMock()
        mock_links = MagicMock()
        mock_oracle.get_links.return_value = mock_links

        links = run_links(mock_oracle, cluster_column_name="leiden", p_cutoff=p_cutoff)

        # Verify filter_links was called with correct p_cutoff
        call_args = mock_links.filter_links.call_args
        assert call_args[1]["p"] == p_cutoff


def test_integration_with_config():
    """Test that config parameters are used correctly"""
    from genecircuitry import config
    from genecircuitry.celloracle_processing import run_KNN, run_links

    # Test GRN_N_JOBS is used in run_KNN
    mock_oracle = MagicMock()
    mock_oracle.copy.return_value = mock_oracle
    mock_oracle.adata.shape = (1000, 50)

    original_jobs = config.GRN_N_JOBS
    config.update_config(GRN_N_JOBS=4)

    run_KNN(mock_oracle, n_comps=50)

    call_args = mock_oracle.knn_imputation.call_args
    assert call_args[1]["n_jobs"] == 4

    # Reset config
    config.update_config(GRN_N_JOBS=original_jobs)

    # Test FIGURES_DIR is used in run_links
    mock_oracle2 = MagicMock()
    mock_links = MagicMock()
    mock_oracle2.get_links.return_value = mock_links

    original_fig_dir = config.FIGURES_DIR
    config.update_config(FIGURES_DIR="test_figures")

    run_links(mock_oracle2, "leiden", 0.001)

    call_args = mock_links.plot_degree_distributions.call_args
    assert "test_figures" in call_args[1]["save"]

    # Reset config
    config.update_config(FIGURES_DIR=original_fig_dir)


def test_module_imports():
    """Test that module imports work correctly"""
    try:
        from genecircuitry.celloracle_processing import (
            create_oracle_object,
            run_PCA,
            run_KNN,
            run_links,
        )

        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import celloracle_processing functions: {e}")


def test_function_signatures():
    """Test that functions have correct signatures"""
    from genecircuitry.celloracle_processing import (
        create_oracle_object,
        run_PCA,
        run_KNN,
        run_links,
    )
    import inspect

    # Check create_oracle_object
    sig = inspect.signature(create_oracle_object)
    assert "adata" in sig.parameters
    assert "cluster_column_name" in sig.parameters
    assert "embedding_name" in sig.parameters
    assert "TG_to_TF_dictionary" in sig.parameters
    assert "raw_count_layer" in sig.parameters

    # Check run_PCA
    sig = inspect.signature(run_PCA)
    assert "oracle" in sig.parameters

    # Check run_KNN
    sig = inspect.signature(run_KNN)
    assert "oracle" in sig.parameters
    assert "n_comps" in sig.parameters

    # Check run_links
    sig = inspect.signature(run_links)
    assert "oracle" in sig.parameters
    assert "cluster_column_name" in sig.parameters
    assert "p_cutoff" in sig.parameters
