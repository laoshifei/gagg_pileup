from pathlib import Path
import os, sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_RECOVERY_ENERGY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_recovery_energy_report_v3 import main
    main()
