from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_HISTORY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_history_study import main
    main()
