from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('MEANSCORE','SELECTION','CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_mean_score_delivery import notebook,docs,validate
    if len(sys.argv)>1 and sys.argv[1]=='validate':validate()
    else:notebook();docs();validate()
