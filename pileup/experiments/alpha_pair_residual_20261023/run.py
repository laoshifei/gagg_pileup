from pathlib import Path
import os
import sys

folder=Path(__file__).resolve().parent
project=folder.parents[2]
for key in ('GAGG_RESIDUAL_PROJECT','GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(project)
source=folder/'inference-source'
sys.path.insert(0,str(source if source.exists() else project))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_residual_study import main
    main()
