from pathlib import Path
import importlib.util,os,sys
folder=Path(__file__).resolve().parent
for key in ('JOINT','NULL','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'development-source-v2'))
if __name__=='__main__':
    command=sys.argv[1]
    if command not in ('figures','validate','notebook'):raise SystemExit('figures | validate | notebook')
    stem='alpha_pair_joint_notebook' if command=='notebook' else 'alpha_pair_joint_report'
    spec=importlib.util.spec_from_file_location('pileup.analysis.'+stem,folder/'report-source'/(stem+'.py'))
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    getattr(module,'deliver' if command=='notebook' else command)()
