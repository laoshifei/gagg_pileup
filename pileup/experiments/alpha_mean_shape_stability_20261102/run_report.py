from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_SHAPE_STABILITY_ROOT']=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_mean_shape_stability_report import main
    main()
