from pathlib import Path
import os,sys,json,hashlib,importlib.util
folder=Path(__file__).resolve().parent
root=folder.parents[2]
os.environ['GAGG_SPATIAL_ROOT']=str(root)
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg-spatial-mpl')
for path,expected in json.loads((folder/'report-freeze.json').read_text()).items():
    assert hashlib.sha256((root/path).read_bytes()).hexdigest()==expected,path
sys.path.insert(0,str(folder/'source'))
spec=importlib.util.spec_from_file_location('pileup.analysis.alpha_pair_spatial_report',folder/'report/alpha_pair_spatial_report.py')
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)
if __name__=='__main__':module.main()
