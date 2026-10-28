from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'coverage-source'))
pileup.analysis.__path__.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    from pileup.analysis import alpha_pair_contour_report as implementation
    assert len(sys.argv)==2 and sys.argv[1] in ['analyze', 'figures']
    getattr(implementation,sys.argv[1])()
