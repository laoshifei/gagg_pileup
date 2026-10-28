from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys
folder=Path(__file__).resolve().parent
for name in ('GATE','RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('notebook','documents','validate'));args=parser.parse_args()
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for manifest in ('freeze','report-freeze','delivery-freeze'):
        for name,value in json.loads((folder/(manifest+'.json')).read_text())['files'].items():assert digest(folder/name)==value,name
    name='alpha_pair_physical_gate_delivery'
    spec=importlib.util.spec_from_file_location('pileup.analysis.'+name,folder/'delivery-source'/(name+'.py'))
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    getattr(module,args.command)()
