"""Evaluate saved outputs with the frozen analysis source."""
from pathlib import Path
import os
import sys

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
for key in ('GAGG_RESIDUAL_PROJECT', 'GAGG_CONDITIONAL_PROJECT',
            'GAGG_RESOLUTION_PROJECT', 'GAGG_PAIR_PROJECT', 'GAGG_COMPONENT_PROJECT'):
    os.environ[key] = str(PROJECT)
source = FOLDER / 'inference-source'
sys.path.insert(0, str(source if source.exists() else PROJECT))

if __name__ == '__main__':
    from pileup.analysis.alpha_pair_residual_report import report, figures
    if sys.argv[1] == 'figures':
        figures()
    else:
        report(sys.argv[1])
