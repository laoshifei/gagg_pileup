"""Append this experiment's report using its frozen reporting source."""
from pathlib import Path
import os
import sys

here = Path(__file__).resolve().parent
os.environ['GAGG_CURRENT_PAIR_ENERGY_ROOT'] = str(here.parents[2])
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, str(here / 'source'))

if __name__ == '__main__':
    from pileup.analysis.alpha_current_pair_energy_report import main
    main()
