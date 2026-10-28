from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys
folder=Path(__file__).resolve().parent
for name in ('EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('sweep','analyze','figures'));args=parser.parse_args()
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name,value in json.loads((folder/'report-freeze.json').read_text())['files'].items():assert digest(folder/name)==value,name
    for name in ('alpha_pair_risk_report','alpha_pair_evidence_report'):
        spec=importlib.util.spec_from_file_location('pileup.analysis.'+name,folder/'report-source'/(name+'.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    getattr(module,args.command)()
