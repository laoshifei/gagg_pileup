from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
for name in ('GAGG_CURRENT_CONFIRMATION_ROOT','GAGG_CURRENT_RESPONSE_ROOT'):
    os.environ[name]=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_current_confirmation_report import main
    main()
