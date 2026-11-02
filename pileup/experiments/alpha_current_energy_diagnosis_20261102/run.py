from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_CURRENT_ENERGY_DIAGNOSIS_ROOT']=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_current_energy_diagnosis_study import main
    main()
