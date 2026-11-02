from pathlib import Path
import os,sys,traceback,json
f=Path(__file__).resolve().parent
root=f.parents[2]
os.environ['GAGG_CURRENT_PAIR_ENERGY_ROOT']=str(root)
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(f/'source'))
if __name__=='__main__':
    rec=root/'jupyter/run_records'/f.name
    with (rec/'pair-resume.log').open('w',buffering=1) as log:
        os.dup2(log.fileno(),1)
        os.dup2(log.fileno(),2)
        try:
            from pileup.analysis.alpha_current_pair_energy_resume import main
            main()
        except BaseException:
            traceback.print_exc()
            (rec/'resume-error.json').write_text(json.dumps({'failed':True,'traceback':traceback.format_exc()},indent=2))
            raise
