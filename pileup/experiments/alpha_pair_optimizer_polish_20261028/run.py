from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
root=str(folder.parents[2])
os.environ['GAGG_POLISH_PROJECT']=root
os.environ['GAGG_UNCERTAINTY_PROJECT']=root
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_optimizer_polish_study import main
    main()
