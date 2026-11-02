from pathlib import Path
import os,sys,importlib.util
f=Path(__file__).resolve().parent
os.environ['GAGG_SPATIAL_ENERGY_ROOT']=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
import pileup.analysis
spec=importlib.util.spec_from_file_location('pileup.analysis.alpha_pair_spatial_energy_report',f/'report/alpha_pair_spatial_energy_report.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
if __name__=='__main__': m.main()
