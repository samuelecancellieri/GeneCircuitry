"""
Tests for warning suppression and propagation across processes.
"""

import os
import subprocess
import sys
import warnings
import genecircuitry


def test_pythonwarnings_env_propagated():
    """Ensure PYTHONWARNINGS environment variable includes pkg_resources filter."""
    pythonwarnings = os.environ.get("PYTHONWARNINGS", "")
    assert "pkg_resources is deprecated as an API" in pythonwarnings


def test_pkg_resources_warning_ignored_in_process():
    """Ensure pkg_resources warning is filtered in current process."""
    with warnings.catch_warnings(record=True) as recorded:
        warnings.warn(
            "pkg_resources is deprecated as an API. See https://setuptools.pypa.io/",
            UserWarning,
        )
    assert len(recorded) == 0


def test_pkg_resources_warning_ignored_in_subprocess():
    """Ensure child processes inherit warning suppression via PYTHONWARNINGS."""
    cmd = [
        sys.executable,
        "-c",
        "import warnings; warnings.warn('pkg_resources is deprecated as an API.', UserWarning)",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert "pkg_resources is deprecated as an API" not in res.stderr
