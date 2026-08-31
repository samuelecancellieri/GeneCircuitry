"""
GeneCircuitry CLI entry point for `python -m genecircuitry.pipeline`
"""

import sys
from .controller import main

if __name__ == "__main__":
    sys.exit(main())

