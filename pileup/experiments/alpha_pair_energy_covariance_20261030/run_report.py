from pathlib import Path
import os,sys
f=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_COVARIANCE_PROJECT']=str(f.parents[2])
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    import runpy
    runpy.run_path(str(f/'delivery/alpha_pair_energy_covariance_report.py'),run_name='pileup.analysis.alpha_pair_energy_covariance_report')
    from importlib.util import spec_from_file_location,module_from_spec
    spec=spec_from_file_location('pileup.analysis._covariance_report',f/'delivery/alpha_pair_energy_covariance_report.py')
    mod=module_from_spec(spec);spec.loader.exec_module(mod);mod.main()
