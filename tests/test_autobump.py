"""
Tests for autobump bioconda script functions.
"""

from pathlib import Path
import sys

# Ensure scripts directory is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from autobump_bioconda import update_meta_yaml_content


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
