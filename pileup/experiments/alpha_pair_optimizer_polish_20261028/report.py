from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
root=str(folder.parents[2])
os.environ['GAGG_POLISH_PROJECT']=root
os.environ['GAGG_UNCERTAINTY_PROJECT']=root
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_optimizer_polish_report import main,check_extra
    check_extra();main()
