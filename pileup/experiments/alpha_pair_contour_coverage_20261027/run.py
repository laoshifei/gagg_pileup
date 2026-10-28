from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_contour_study import main
    main()
