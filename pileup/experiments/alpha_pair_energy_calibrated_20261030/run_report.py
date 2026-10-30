from pathlib import Path
import os,sys,importlib.util
f=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_CALIBRATED_PROJECT']=str(f.parents[2])
sys.path.insert(0,str(f/'source'))
spec=importlib.util.spec_from_file_location('pileup.analysis.alpha_pair_energy_calibrated_report',f/'report-source/alpha_pair_energy_calibrated_report.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
if __name__=='__main__':m.main()
