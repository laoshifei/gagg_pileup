"""Independent energy calibration with native guard bands and interpolation tests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import sys

import numpy as np

from .sipm_energy_calibration import METHODS, fit_calibration, invert_charge, selection_score, error_summary, score_events
from .sipm_energy_guard_data import PROJECT, NAME, EXPERIMENT, OUTPUT, RECORD, OLD, load, save, digest, now, check_freeze, scalar_rows

FIT_METHODS=['model_relative_residual','joint_conserved_nll']
DIAGNOSTIC='true_charge_diagnostic'


def lines(path):return [json.loads(s) for s in Path(path).read_text().splitlines()]
def write_lines(path,rows):
    with Path(path).open('w') as f:
        for r in rows:f.write(json.dumps(r,ensure_ascii=False,allow_nan=False)+'\n')


def check_study():
    check_freeze()
    for p,sha in load(EXPERIMENT/'analysis-freeze.json')['files'].items():assert digest(EXPERIMENT/p)==sha,p


def freeze():
    check_freeze();assert load(RECORD/'pilot-validation.json')['passed']
    assert not any(r['stage']!='pilot' for r in load(EXPERIMENT/'manifest.json')['runs'])
    assert not (EXPERIMENT/'analysis-freeze.json').exists()
    src=EXPERIMENT/'analysis-source/pileup/analysis/sipm_energy_guard_study.py'
    shutil.copyfile(PROJECT/'pileup/analysis/sipm_energy_guard_study.py',src)
    (EXPERIMENT/'energy.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_GUARD_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'analysis-source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_guard_study import main
    main()
''')
    save(EXPERIMENT/'analysis-freeze.json',dict(created=now(),new_development_started=False,new_confirmation_started=False,
        files={str(p.relative_to(EXPERIMENT)):digest(p) for p in [src,EXPERIMENT/'energy.py',EXPERIMENT/'freeze.json']}))
    print('GUARD_ANALYSIS_FROZEN',flush=True)


def table():
    text=(EXPERIMENT/'source/brems/g4/src/gagg_birks_yield_curves.inc').read_text()
    def vector(name):
        s=re.search(r'\b'+name+r'\s*=\s*\{([^}]*)\}',text,re.S).group(1)
        return np.array([float(v.strip().replace('*MeV','')) for v in s.split(',') if v.strip()])
    return vector('gaggScintillationEnergy'),vector('gaggAlphaScintillationYield')


def calibrate():
    check_study();assert not (EXPERIMENT/'energy-model.json').exists()
    old=load(OUTPUT/'old-calibration-records.json');new_fit=scalar_rows('fit');new_select=scalar_rows('select')
    assert len(old)==560 and len(new_fit)==96 and len(new_select)==64
    fitting=[dict(r,stage='fit',root_sha256=r['root_original_sha256']) for r in old if r['fold']==0]+new_fit
    selecting=[dict(r,stage='select',root_sha256=r['root_original_sha256']) for r in old if r['fold']==1]+new_select
    assert len(fitting)==376 and len(selecting)==344
    assert all(r['edep_MeV']>0 and r['charge_pe']>0 for r in fitting+selecting)
    for field in ('root','seed','root_sha256'):
        assert not {r[field] for r in fitting}&{r[field] for r in selecting},field
    assert len({r['root'] for r in fitting+selecting})==720
    assert {r['source_energy_MeV'] for r in selecting}=={1.625,2.,2.5,3.5,4.5,5.486,6.5,8.}
    assert not any(r['stage']=='confirmation' for r in load(EXPERIMENT/'manifest.json')['runs'])
    save(OUTPUT/'calibration.records.json',fitting+selecting)
    def arrays(rows):return tuple(np.array([r[k] for r in rows]) for k in ['edep_MeV','charge_pe','source_energy_MeV'])
    e,q,g=arrays(fitting);se,sq,sg=arrays(selecting);te,tl=table();models={};choices=[]
    for method in METHODS:
        try:
            model=fit_calibration(e,q,g,te,tl,method)
            loss,groups=selection_score(model,se,sq,sg)
            models[method]=model;choices.append(dict(method=method,valid=True,selection_loss=loss,groups=groups))
        except (ValueError,RuntimeError) as exc:choices.append(dict(method=method,valid=False,error=str(exc)))
    valid=[r for r in choices if r['valid']];assert valid,'No valid calibration; preserve diagnostics'
    chosen=min(valid,key=lambda r:r['selection_loss'])['method']
    save(OUTPUT/'candidate-models.json',models);save(EXPERIMENT/'energy-model.json',models[chosen])
    save(EXPERIMENT/'mapping-selection.json',dict(created=now(),chosen=chosen,candidates=choices,
        fit_parents=376,select_parents=344,selection_source_groups=sorted(set(sg.tolist())),
        calibration_sha256=digest(OUTPUT/'calibration.records.json'),no_confirmation_generated=True,refit_on_selection=False))
    save(EXPERIMENT/'energy-freeze.json',dict(created=now(),confirmation_generated=False,
        files={name:digest(EXPERIMENT/name) for name in ['energy-model.json','mapping-selection.json','analysis-freeze.json','freeze.json']}))
    print('GUARD_MAPPING_SELECTED',chosen,[(r['method'],r.get('selection_loss')) for r in choices],flush=True)


def export():
    """Generation-side export: true primary Q reference, not a blind waveform estimate."""
    check_study()
    for p,sha in load(EXPERIMENT/'energy-freeze.json')['files'].items():assert digest(EXPERIMENT/p)==sha,p
    dest=OUTPUT/'single-charge-observations.jsonl';assert not dest.exists()
    records=sorted([r for r in load(OUTPUT/'native-records.json')['test_single'] if r['stage']=='confirmation'],key=lambda x:x['seed'])
    assert len(records)==320
    observations=[]
    for i,r in enumerate(records):
        with np.load(PROJECT/r['truth']) as f:
            q=float(f['charge_pe'][0])
        observations.append(dict(index=i,charge_pe=q))
    write_lines(dest,observations)
    shutil.copyfile(OLD/'fits.observations.jsonl',OUTPUT/'historical-fits.observations.jsonl')
    save(RECORD/'charge-export.json',dict(created=now(),single_count=320,
        single_input='True isolated primary-descendant Q; calibration transfer diagnostic, not blind waveform recovery',
        exported_fields=['index','charge_pe'],energy_fields_read=False,
        single_charge_sha256=digest(dest),historical_fits_sha256=digest(OUTPUT/'historical-fits.observations.jsonl'),
        source_native_index_sha256=digest(OUTPUT/'native-records.json')))
    print('GUARD_CHARGES_EXPORTED',320,flush=True)


def predict():
    check_study()
    for p,sha in load(EXPERIMENT/'energy-freeze.json')['files'].items():assert digest(EXPERIMENT/p)==sha,p
    audit=load(RECORD/'charge-export.json')
    assert audit['single_charge_sha256']==digest(OUTPUT/'single-charge-observations.jsonl')
    assert audit['historical_fits_sha256']==digest(OUTPUT/'historical-fits.observations.jsonl')
    assert not (OUTPUT/'single-energy.predictions.jsonl').exists()
    active=[True];opened=[]
    def guard(event,args):
        if active[0] and event=='open' and isinstance(args[0],(str,bytes)):
            path=str(args[0]);opened.append(path)
            if any(s in path for s in ['truth','labels','calibration.records','native-inputs']) or path.endswith('.root'):
                raise PermissionError('Energy inference cannot access labels or ROOT')
    sys.addaudithook(guard)
    for name in ['forbidden.truth.npz','forbidden.labels.json','forbidden.root']:
        try:open(name)
        except PermissionError:pass
        else:raise AssertionError('Label guard inactive')
    model=load(EXPERIMENT/'energy-model.json')
    def finite(v):return float(v) if np.isfinite(v) else None
    singles=[dict(r,energy_MeV=finite(invert_charge(model,r['charge_pe']))) for r in lines(OUTPUT/'single-charge-observations.jsonl')]
    blind=[dict(r,energy_MeV=[finite(v) for v in invert_charge(model,r['charge_pe'])]) for r in lines(OUTPUT/'historical-fits.observations.jsonl')]
    write_lines(OUTPUT/'single-energy.predictions.jsonl',singles);write_lines(OUTPUT/'historical-energy.predictions.jsonl',blind)
    active[0]=False
    save(RECORD/'energy-inference.json',dict(created=now(),passed=True,blocked_truth_probes=3,opened_paths=opened,
        model_sha256=digest(EXPERIMENT/'energy-model.json'),single_count=len(singles),historical_blind_count=len(blind),
        predictions={name:digest(OUTPUT/name) for name in ['single-energy.predictions.jsonl','historical-energy.predictions.jsonl']},
        confirmation_energies_read=False,single_mode='Known isolated primary Q reference, not blind',historical_waveform_fits_unchanged=True))
    print('GUARD_PREDICTIONS_FROZEN',len(singles),len(blind),flush=True)


def summarize(rows):
    errs=[m['relative_energy_error'] for r in rows for m in r['matches'] if m['relative_energy_error'] is not None]
    target=sum(r['target'] for r in rows)
    return dict(method=rows[0]['method'],kind=rows[0]['kind'],condition=rows[0]['condition'],gap_ns=rows[0]['gap_ns'],
        frames=len(rows),target_events=target,matched_events=sum(r['matched'] for r in rows),
        missed=sum(r['missed'] for r in rows),extra=sum(r['extra'] for r in rows),
        split=sum(r.get('split',0) for r in rows),merge=sum(r.get('merge',0) for r in rows),
        out_of_domain_predictions=sum(r['energy_out_of_domain'] for r in rows),
        recovered_5=sum(r['strict_5'] for r in rows),recovered_10=sum(r['strict_10'] for r in rows),recovered_20=sum(r['strict_20'] for r in rows),
        energy_efficiency_10=sum(r['strict_10'] for r in rows)/target if target else None,
        frame_success_10=sum(r['frame_success_10'] for r in rows),frame_efficiency_10=np.mean([r['frame_success_10'] for r in rows]).item(),
        median_relative=float(np.median(errs)) if errs else None,
        sigma68_relative=float((np.quantile(errs,.84)-np.quantile(errs,.16))/2) if errs else None)


def evaluate():
    check_study();audit=load(RECORD/'energy-inference.json')
    assert audit['passed'] and audit['model_sha256']==digest(EXPERIMENT/'energy-model.json')
    for name,sha in audit['predictions'].items():assert digest(OUTPUT/name)==sha
    assert not (OUTPUT/'report.json').exists()
    rows=sorted(scalar_rows('confirmation'),key=lambda r:r['seed']);assert len(rows)==320
    save(OUTPUT/'confirmation.labels.json',rows)
    cal=load(OUTPUT/'calibration.records.json');pil=load(OUTPUT/'pilot.records.json')
    for key in ['root','seed','root_sha256']:
        assert len({r[key] for r in rows})==320
        assert not {r[key] for r in rows}&{r[key] for r in cal+pil}
    pred=lines(OUTPUT/'single-energy.predictions.jsonl');assert len(pred)==len(rows)
    q=np.array([r['charge_pe'] for r in rows]);e=np.array([r['edep_MeV'] for r in rows]);g=np.array([r['source_energy_MeV'] for r in rows])
    assert all(p['index']==i and p['charge_pe']==q[i] for i,p in enumerate(pred))
    estimated=np.array([p['energy_MeV'] for p in pred],float)
    oldmodel=load(PROJECT/'pileup/experiments/sipm_energy_calibration_v2_20261028/energy-model.json')
    old_est= invert_charge(oldmodel,q);single=[];rng=np.random.default_rng(202610290051)
    for label,values in [('new_guard_calibration',estimated),('old_four_band_calibration',old_est)]:
        for level in sorted(set(g)):
            mask=g==level;s=error_summary(e[mask],values[mask]);finite=np.isfinite(values[mask]);rel=(values[mask]-e[mask])/e[mask]
            valid_rel=rel[finite]
            ci=np.quantile(np.median(rng.choice(valid_rel,size=(2000,len(valid_rel))),axis=1),[.025,.975]).tolist() if len(valid_rel) else [None,None]
            single.append(dict(method=label,source_energy_MeV=float(level),**s,
                median_bootstrap95=ci,edep_min_MeV=float(e[mask].min()),edep_max_MeV=float(e[mask].max()),
                charge_min_pe=float(q[mask].min()),charge_max_pe=float(q[mask].max()),
                only_in_confirmation=bool(level in [1.8,7.]),
                control_pass=bool(s['valid']/s['n']>=.95 and s['median_relative'] is not None and abs(s['median_relative'])<=.05)))
    passed=all(s['control_pass'] for s in single if s['method']=='new_guard_calibration')
    # Historical paired replay; energies and previous scores were already exposed.
    cases=load(OLD/'evaluation.cases.json')
    with np.load(OLD/'evaluation.labels.npz') as f:true_e=f['energy_MeV'];true_q=f['charge_pe'];true_t=f['times_ns'];block=f['block']
    old_scores={(r['method'],r['index']):r for r in lines(OLD/'energy.scores.jsonl')}
    scored=[]
    for p in lines(OUTPUT/'historical-energy.predictions.jsonl'):
        i=p['index'];c=cases[i];before=old_scores[p['method'],i]
        s=score_events(true_t[i],true_q[i],true_e[i],p['times_ns'],p['energy_MeV'],c['gap_ns'])
        assert all(s[k]==before[k] for k in ['matched','missed','extra','target'])
        scored.append(dict(s,index=i,method=p['method'],kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],block=int(block[i]),
            split=before['split'],merge=before['merge']))
    model=load(EXPERIMENT/'energy-model.json')
    for i,c in enumerate(cases):
        mask=true_q[i]>=500.;est=invert_charge(model,true_q[i,mask])
        s=score_events(true_t[i],true_q[i],true_e[i],true_t[i,mask],est,c['gap_ns'])
        scored.append(dict(s,index=i,method=DIAGNOSTIC,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],block=int(block[i])))
    write_lines(OUTPUT/'historical-energy.scores.jsonl',scored)
    groups={}
    for r in scored:groups.setdefault((r['method'],r['kind'],r['condition'],r['gap_ns']),[]).append(r)
    summaries=[summarize(rr) for rr in groups.values()]
    weights=np.zeros((2000,40),int);rng=np.random.default_rng(202610283501)
    for start in [0,20]:
        draws=rng.integers(start,start+20,size=(2000,20))
        for i in range(2000):weights[i]+=np.bincount(draws[i],minlength=40)
    aggregate={};boots={}
    for method in FIT_METHODS+[DIAGNOSTIC]:
        aggregate[method]={};boots[method]={}
        for weak,label in [(True,'weak'),(False,'strong')]:
            means=[];draws=[]
            for key,rr in groups.items():
                if key[0]!=method or key[1]!='pair' or ('weak' in key[2])!=weak:continue
                ids=np.array([r['block'] for r in rr]);values=np.array([r['frame_success_10'] for r in rr],float)
                num=np.bincount(ids,weights=values,minlength=40);den=np.bincount(ids,minlength=40)
                means.append(values.mean());draws.append((weights@num)/(weights@den))
            aggregate[method][label]=float(np.mean(means));boots[method][label]=np.mean(draws,axis=0)
    paired={}
    for label in ['weak','strong']:
        diff=boots[FIT_METHODS[1]][label]-boots[FIT_METHODS[0]][label]
        paired[label]=dict(difference=aggregate[FIT_METHODS[1]][label]-aggregate[FIT_METHODS[0]][label],mother_block95=np.quantile(diff,[.025,.975]).tolist())
    save(OUTPUT/'report.json',dict(created=now(),chosen=model['method'],single_confirmation=single,energy_transfer_pass=passed,
        independent_new_single_parents=320,calibration_parents=376,selection_parents=344,pilot_parents=16,
        interpolation_test_energies_MeV=[1.8,7.],mapping_selection=load(EXPERIMENT/'mapping-selection.json'),
        historical_summary=summaries,historical_aggregate=aggregate,historical_paired=paired,
        previous_aggregate=load(OLD/'report.json')['aggregate_strict_frame10'],
        new_single_mode='Truth-derived primary Q calibration transfer; not blind event detection',
        historical_mode='Previously exposed fixed fits; new energy mapping only, not a new pileup confirmation',
        new_waveform_fits=False,default_fitter_changed=False,
        width='Conditional sigma68 relative energy residual; not measured spectrum FWHM'))
    save(RECORD/'science-validation.json',dict(passed=True,created=now(),new_confirmation_after_model_freeze=True,
        new_mothers_disjoint_by=['ROOT path','original ROOT SHA256','source seed'],native_roots=496,
        independent_calibration_confirmation_parents=1056,
        historical_blind_detection_counts_unchanged=3200,truth_guard_probes=3,
        prediction_hashes=audit['predictions'],report_sha256=digest(OUTPUT/'report.json'),
        calibration_sha256=digest(OUTPUT/'calibration.records.json'),confirmation_sha256=digest(OUTPUT/'confirmation.labels.json')))
    print('GUARD_EVALUATION',dict(chosen=model['method'],pass_controls=passed,new_single_bias={str(r['source_energy_MeV']):r['median_relative'] for r in single if r['method']=='new_guard_calibration'},
                               historical=aggregate,paired=paired),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['freeze','calibrate','export','predict','evaluate'])
    globals()[parser.parse_args().command]()


if __name__=='__main__':main()
