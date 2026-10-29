from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_GUARD_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'analysis-source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_guard_study import main
    main()
