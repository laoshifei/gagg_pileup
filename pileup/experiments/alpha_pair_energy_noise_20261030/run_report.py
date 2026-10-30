from pathlib import Path
import os,sys,importlib.util
f=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_NOISE_PROJECT']=str(f.parents[2])
sys.path.insert(0,str(f/'source-v2'))
spec=importlib.util.spec_from_file_location('pileup.analysis.alpha_pair_energy_noise_report',f/'delivery/alpha_pair_energy_noise_report.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
if __name__=='__main__': module.make_report()
