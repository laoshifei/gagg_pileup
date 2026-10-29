from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_UNCERTAINTY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_blind_uncertainty_delivery import check_extra
    check_extra()
    from pileup.analysis.alpha_pair_blind_uncertainty_delivery import main
    main()
