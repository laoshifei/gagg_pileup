"""Append verified results to local notebooks, preserving existing cells."""
from pathlib import Path
import os
import sys

FOLDER = Path(__file__).resolve().parent
for key in ('GAGG_RESIDUAL_PROJECT', 'GAGG_CONDITIONAL_PROJECT',
            'GAGG_RESOLUTION_PROJECT', 'GAGG_PAIR_PROJECT', 'GAGG_COMPONENT_PROJECT'):
    os.environ[key] = str(FOLDER.parents[2])
sys.path.insert(0, str(FOLDER / 'inference-source'))

if __name__ == '__main__':
    from pileup.analysis.alpha_pair_residual_notebook import deliver, catalog
    deliver()
    catalog()
