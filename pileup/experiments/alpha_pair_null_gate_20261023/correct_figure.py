from pathlib import Path
import importlib.util
import os
import sys

folder=Path(__file__).resolve().parent
os.environ['GAGG_NULL_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    spec=importlib.util.spec_from_file_location('pileup.analysis.frozen_null_rank_figure',
        folder/'postprocessing-source/alpha_pair_null_calibration_figure.py')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    module.draw()
