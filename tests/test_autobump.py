"""
Tests for autobump bioconda script functions.
"""

from pathlib import Path
import sys

# Ensure scripts directory is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from autobump_bioconda import submit_bioconda_pr, update_meta_yaml_content


def test_update_meta_yaml_content():
    """Test regex replacements in meta.yaml content."""
    sample_yaml = """{% set name = "genecircuitry" %}
{% set version = "0.2.1" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/g/genecircuitry/genecircuitry-0.2.1.tar.gz
  sha256: oldhash123

build:
  number: 3
  noarch: python
"""
    updated = update_meta_yaml_content(sample_yaml, "0.2.2", "newhash456")
    assert '{% set version = "0.2.2" %}' in updated
    assert "sha256: newhash456" in updated
    assert "number: 0" in updated


def test_submit_bioconda_pr_dry_run(tmp_path=None):
    """Test submit_bioconda_pr in dry-run mode."""
    recipe_path = Path(__file__).parent.parent / "conda-recipe" / "meta.yaml"
    result = submit_bioconda_pr(
        package_name="genecircuitry",
        version="0.2.2",
        sha256="testsha256",
        meta_yaml_source_path=recipe_path,
        dry_run=True,
    )
    assert result is None


def test_submit_bioconda_pr_no_token():
    """Test submit_bioconda_pr when no token is provided."""
    import os
    # Temporarily clear token env vars
    orig_bioconda = os.environ.pop("BIOCONDA_TOKEN", None)
    orig_gh = os.environ.pop("GH_TOKEN", None)
    try:
        recipe_path = Path(__file__).parent.parent / "conda-recipe" / "meta.yaml"
        result = submit_bioconda_pr(
            package_name="genecircuitry",
            version="0.2.2",
            sha256="testsha256",
            meta_yaml_source_path=recipe_path,
            token=None,
            dry_run=False,
        )
        assert result is None
    finally:
        if orig_bioconda is not None:
            os.environ["BIOCONDA_TOKEN"] = orig_bioconda
        if orig_gh is not None:
            os.environ["GH_TOKEN"] = orig_gh


if __name__ == "__main__":
    test_update_meta_yaml_content()
    test_submit_bioconda_pr_dry_run()
    test_submit_bioconda_pr_no_token()
    print("All autobump tests passed successfully!")

