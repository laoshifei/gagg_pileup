"""Frozen runner with metadata-only NumPy scalar serialization adapter."""
from pathlib import Path
import os,sys,json,hashlib
folder=Path(__file__).resolve().parent
root=folder.parents[2]
os.environ['GAGG_SPATIAL_ROOT']=str(root)
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
for path, expected in json.loads((folder/'runner-v3-freeze.json').read_text()).items():
    assert hashlib.sha256((root/path).read_bytes()).hexdigest()==expected,path
sys.path.insert(0,str(folder/'source'))
from pileup.analysis import alpha_pair_spatial_study as study
import numpy as np
original_save=study.save
def numpy_scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
def metadata_save(path,data):
    original_save(path,json.loads(json.dumps(data,default=numpy_scalar,allow_nan=False)))
study.save=metadata_save
study.compatibility()
original_point=study.map_point
def checkpoint_point(job):
    name,amplitude=job
    path=study.REC/'mapping-points'/f'{name}-{amplitude:.12g}.json'
    signature=hashlib.sha256((folder/'source-config-freeze.json').read_bytes()).hexdigest()
    if path.exists():
        cached=json.loads(path.read_text())
        assert cached['source_signature']==signature
        assert cached['point']['light_pe']==amplitude
        return cached['point']
    row=original_point(job)
    metadata_save(path,dict(source_signature=signature,point=row))
    return row
study.map_point=checkpoint_point
if __name__=='__main__':study.main()
