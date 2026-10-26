from pathlib import Path
import hashlib, json, os, sys
folder=Path(__file__).resolve().parent
for name in ('EMPIRICAL','WEAK','NULL','CONDITIONAL','RESIDUAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    for name,sha in json.loads((folder/'report-freeze.json').read_text())['files'].items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==sha,name
    from pileup.analysis import alpha_pair_weak_report as report
    command=sys.argv[1] if len(sys.argv)>1 else 'all'
    if command in ('all','validate'):report.validate()
    if command in ('all','figures'):
        report.figures_diagnostic();report.figures_empirical()
    if command=='notebook':
        from pileup.analysis.alpha_pair_weak_notebook import deliver
        deliver()
    if command not in ('all','validate','figures','notebook'):raise ValueError(command)
