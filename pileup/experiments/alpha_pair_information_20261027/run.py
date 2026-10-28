from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_INFORMATION_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_information_study import main
    main()
