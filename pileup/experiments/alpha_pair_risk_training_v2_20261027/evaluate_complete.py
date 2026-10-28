from pathlib import Path
import argparse,hashlib,json,os,sys
folder=Path(__file__).resolve().parent
for name in ('RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'evaluation-source'))
import pileup.analysis
pileup.analysis.__path__.append(str(folder/'dependency-source/pileup/analysis'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=('native-features','baseline','predict'));p.add_argument('--domain',choices=('native','continuous'),default='native');p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    manifest=json.loads((folder/'evaluation-freeze.json').read_text());digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert digest(folder/'freeze.json')==manifest['training_source_freeze_sha256']
    for name,value in manifest['files'].items():assert digest(folder/name)==value,name
    extension=json.loads((folder/'dependency-freeze.json').read_text())
    assert digest(folder/'evaluation-freeze.json')==extension['evaluation_freeze_sha256']
    for name,value in extension['files'].items():assert digest(folder/name)==value,name
    from pileup.analysis import alpha_pair_risk_inference as module
    if a.command=='native-features':module.native_features(a.workers)
    elif a.command=='baseline':module.baseline(a.workers)
    else:module.predict(a.domain,a.workers)
