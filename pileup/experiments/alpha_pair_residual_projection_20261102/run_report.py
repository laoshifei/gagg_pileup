from pathlib import Path
import os,sys,runpy
f=Path(__file__).resolve().parent
os.environ['GAGG_RESIDUAL_ROOT']=str(f.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    runpy.run_path(str(f/'report.py'),run_name='__main__')
