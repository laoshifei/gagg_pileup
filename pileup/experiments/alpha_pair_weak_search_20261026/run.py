from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('SEARCH','JOINT','CONFIRM','NULL','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_weak_search_study import main
    main()
