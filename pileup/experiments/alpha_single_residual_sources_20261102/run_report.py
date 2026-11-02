from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_RESIDUAL_SOURCES_ROOT']=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_single_residual_sources_report_v2 import main
    main()
