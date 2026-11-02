from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_TRANSPORT_ENV_ROOT']=str(f.parents[2])
sys.path.insert(0,str(f/'analysis-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_transport_environment_report import main
    main()
