from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_FORWARD_PROJECT']=str(f.parents[2])
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_recovery_forward_report import main
    main()
