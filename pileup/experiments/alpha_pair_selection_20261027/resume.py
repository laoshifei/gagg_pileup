from pathlib import Path
import os,sys,runpy
folder=Path(__file__).resolve().parent
for key in ('SELECTION','CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'resume-source'))
if __name__=='__main__':runpy.run_module('pileup.analysis.alpha_pair_selection_resume',run_name='__main__')
