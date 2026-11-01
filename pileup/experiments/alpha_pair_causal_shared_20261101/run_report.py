from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_CAUSAL_SHARED_ROOT']=str(f.parents[2])
sys.path.insert(0,str(f/'report-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_causal_shared_study import compatibility,report
    compatibility();report()
