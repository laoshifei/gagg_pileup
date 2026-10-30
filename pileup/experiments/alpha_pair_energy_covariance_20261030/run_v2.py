from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_COVARIANCE_PROJECT']=str(f.parents[2])
sys.path.insert(0,str(f/'source-v2'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_energy_covariance_study import main
    main()
