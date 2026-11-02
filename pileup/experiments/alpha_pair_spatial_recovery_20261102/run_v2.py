"""Frozen runner with metadata-only NumPy scalar serialization adapter."""
from pathlib import Path
import os,sys,json,hashlib
folder=Path(__file__).resolve().parent
root=folder.parents[2]
os.environ['GAGG_SPATIAL_ROOT']=str(root)
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
for path, expected in json.loads((folder/'runner-v2-freeze.json').read_text()).items():
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
if __name__=='__main__':study.main()
