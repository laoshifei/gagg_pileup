from pathlib import Path
import importlib.util
import os
import sys

folder=Path(__file__).resolve().parent
for key in ('GAGG_NULL_PROJECT','GAGG_RESIDUAL_PROJECT','GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT',
            'GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    spec=importlib.util.spec_from_file_location('pileup.analysis.frozen_null_validation',
        folder/'postprocessing-source/alpha_pair_null_validate.py')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    module.validate()
