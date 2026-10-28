from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys
folder=Path(__file__).resolve().parent
for name in ('RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'evaluation-f64-source'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('analyze','figures','notebook','documents','validate'));args=parser.parse_args()
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((folder/'report-v3-freeze.json').read_text())
    for name,value in manifest['files'].items():assert digest(folder/name)==value,name
    choices={'analyze':('report','analyze'),'figures':('report','figures'),'notebook':('notebook','deliver'),'documents':('documents','deliver'),'validate':('validate','validate')}
    suffix,function=choices[args.command];name='alpha_pair_risk_'+suffix
    spec=importlib.util.spec_from_file_location('pileup.analysis.'+name,folder/'report-source-v3'/(name+'.py'))
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    getattr(module,function)()
