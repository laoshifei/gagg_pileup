from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys
folder=Path(__file__).resolve().parent
for key in ('STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('notebook', 'documents', 'validate'));args=parser.parse_args()
    for name,value in json.loads((folder/'delivery-freeze.json').read_text())['files'].items():assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==value,name
    for name in ('alpha_pair_stability_delivery',):
        spec=importlib.util.spec_from_file_location('pileup.analysis.'+name,folder/'delivery-source'/(name+'.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    getattr(module,args.command)()
