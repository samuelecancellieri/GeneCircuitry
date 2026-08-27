"""
Tests for package version consistency.
"""

import genecircuitry


def test_version_string():
    """Ensure __version__ is defined and equals 0.2.2."""
    assert hasattr(genecircuitry, "__version__")
    assert genecircuitry.__version__ == "0.2.2"
