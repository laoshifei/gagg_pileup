from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('JOINT','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
os.environ['GAGG_JOINT_FREEZE']='development-freeze-v2.json'
sys.path.insert(0,str(folder/'development-source-v2'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_joint_study import main
    main()
