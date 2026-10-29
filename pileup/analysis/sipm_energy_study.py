"""Frozen SiPM energy-calibration study using existing independent optical parents.

Stages run from the project root: prepare, collect, calibrate, predict, evaluate.
No waveform fitting, G4 production or changes to historical experiments.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np

from .sipm_energy_calibration import METHODS, fit_calibration, invert_charge, error_summary, selection_score, score_events

PROJECT=Path(os.environ.get('GAGG_ENERGY_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='sipm_energy_calibration_v2_20261028'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PARENT='alpha_pair_joint_likelihood_20261026'
CONFIRM='alpha_pair_joint_confirmation_20261026'
NATIVE_SOURCES=['alpha_pair_conditional_20261022','alpha_pair_residual_20261023','alpha_pair_null_gate_20261023/development']
FIT_METHODS=['model_relative_residual','joint_conserved_nll']


def load(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def lines(path):return [json.loads(s) for s in Path(path).read_text().splitlines()]
def write_lines(path,values):
    with Path(path).open('w') as f:
        for v in values:f.write(json.dumps(v,ensure_ascii=False,allow_nan=False)+'\n')


def check_freeze():
    for path,sha in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/path)==sha,path


def prepare():
    for p in (EXPERIMENT,OUTPUT,RECORD):p.mkdir(parents=True,exist_ok=False)
    source=EXPERIMENT/'source/pileup/analysis';source.mkdir(parents=True)
    (source.parent/'__init__.py').write_text('')
    (source/'__init__.py').write_text('')
    for name in ('sipm_energy_calibration.py','sipm_energy_study.py','test_sipm_energy_calibration.py'):
        shutil.copyfile(PROJECT/'pileup/analysis'/name,source/name)
    original=PROJECT/'pileup/experiments'/CONFIRM
    yield_path=Path('source/brems/g4/src/gagg_birks_yield_curves.inc')
    shutil.copyfile(original/yield_path,EXPERIMENT/'material-yield.inc')
    # All optical mother groups must share the same material response table.
    for exp in ['alpha_pair_conditional_20261022','alpha_pair_residual_20261023','alpha_pair_null_gate_20261023']:
        assert digest(PROJECT/'pileup/experiments'/exp/yield_path)==digest(EXPERIMENT/'material-yield.inc'),exp
    metadata=load(PROJECT/'jupyter/artifacts'/CONFIRM/'native-records.json')['test_single'][0]['metadata']
    save(EXPERIMENT/'sensor-config.json',load(PROJECT/metadata)['config'])
    for name in ('settings.json','protocol.json','gate-calibrations.json','legacy-gate-calibrations.json'):
        shutil.copyfile(original/name,EXPERIMENT/('parent-'+name))
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_study import main
    main()
''')
    protocol=dict(created=now(),phase='A: charge-to-energy transfer and existing-grid validation',
        scope='Pure central alpha, fixed95/351ns, 2ns sampling, full2.048us local voltage, unknownK<=2',
        adequate_dynamic_range=True,new_ROOTs=False,new_waveform_fits=False,
        main_fitting_comparison=FIT_METHODS,fit_or_gate_changed=False,
        calibration_source=PARENT,calibration_split='Existing mother replicate//2 modulo2: fold0 fit, fold1 choose mapping',
        expected_calibration_counts=[280,280],expected_evaluation_parents=320,
        evaluation_source=CONFIRM,evaluation_frames=1600,evaluation_gaps_ns=[8.,40.,200.],
        evaluation_status='Previously published Q/t confirmation; new energy endpoint, disjoint from calibration/selection mothers; not a new unseen transport campaign',
        mapping_candidates=list(METHODS),mapping_selection='Equal mean of four source-group relative-MSE; nonfinite inverse incurs loss1, all events retained; ties follow candidate order',
        fitting_weight='Equal source groups, scaled by fixed group mean Q; no per-event observed-Q weight',
        energy_target='Total deposited energy per primary, MeV; not source energy, entry energy or keVee',
        calibration_charge='Primary-descendant realized avalanche charge/q0; excludes dark ancestry',
        mapping_domain='Observed fold0 E_dep min/max; no inverse extrapolation/clipping; no extra Birks correction',
        primary_metric='All target events recovered within10% E_dep and min(10ns,delta/4) for pairs; no extra output for frame success',
        secondary_energy_tolerances_percent=[5,20],threshold_Q_pe=500.,ordinary_matching_ns=10.,
        controls='Report empty and single groups separately; inherited5% extra-output budget; no retuning',
        energy_selection_gate='Held-out single source groups: at least95% in-domain and abs median E bias<=5% in every group; report failure without retuning',
        oracle_diagnostic='True realized Q and true K/t converted only after blind predictions frozen; not blind recovery',
        uncertainties='2000 paired mother-block bootstrap draws; first20/last20 replicate//2 blocks resampled separately; same weights across methods/conditions/gaps',
        bootstrap_seed=202610283501,
        limits=['Sparse four source energies; interpolation between bands not independently validated',
                'Conditional residual sigma68 is not energy-spectrum FWHM',
                'Not mixed alpha/gamma, free-tau, random high-rate K>=3 or measured hardware'],
        next_stage='If energy-transfer controls pass, expand the frozen method condition grid; otherwise diagnose calibration on development data in a new experiment')
    save(EXPERIMENT/'protocol.json',protocol)
    plan=dict(created=now(),main_claim_candidate='Conditional performance and energy/count bias of blind GAGG-SiPM reconstruction',
        current_phase='Frozen energy transfer and three-gap native verification',
        later_phases=['New points in gap/charge/order map','Continuous0.1/1/10MHz unknown-K and energy spectra','Response shift and mixed alpha/gamma with equal species priority'],
        no_universal_time_or_rate_claim=True,no_new_architecture_search=True)
    save(EXPERIMENT/'paper-evaluation-plan.json',plan)
    # Separate frozen blind fitted parameters from labels for energy conversion.
    original_scores=PROJECT/'jupyter/artifacts'/CONFIRM/'scores.jsonl'
    fits=[]
    for r in lines(original_scores):
        if r['method'] in FIT_METHODS:
            f=r['fit'];assert all(0<=t<2048 for t in f['times_ns'])
            fits.append(dict(index=r['index'],method=r['method'],times_ns=f['times_ns'],charge_pe=f['charge_pe']))
    assert len(fits)==3200 and len({(r['index'],r['method']) for r in fits})==3200
    write_lines(OUTPUT/'fits.observations.jsonl',fits)
    save(RECORD/'fit-export.json',dict(parent_score_file=str(original_scores.relative_to(PROJECT)),
        parent_score_sha256=digest(original_scores),fits_sha256=digest(OUTPUT/'fits.observations.jsonl'),
        allowed_fields=['index','method','times_ns','charge_pe'],entries=len(fits),refitted=False,
        original_blind_validation=load(PROJECT/'jupyter/run_records'/CONFIRM/'science-validation.json')['passed']))
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in EXPERIMENT.rglob('*') if p.is_file()}
    save(EXPERIMENT/'freeze.json',dict(created=now(),files=files,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted='Previous research preserved; this snapshot includes current new energy-calibration modules',
        evaluation_energy_not_read=True))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        dependencies={n:version(n) for n in ('numpy','scipy','matplotlib','nbformat','pandas')},
        threads=1,working_directory=str(PROJECT)))
    print('ENERGY_PREPARED frozen protocol/source; 3200 cached blind fits exported',flush=True)


def native_row(record,config):
    meta_path=PROJECT/record['metadata'];truth_path=PROJECT/record['truth'];meta=load(meta_path)
    assert meta['config']==config
    assert meta['root_sha256']==record['root_original_sha256']
    with np.load(truth_path) as f:
        assert len(f['edep_MeV'])==1 and str(f['particle_name'][0])=='alpha'
        e=float(f['edep_MeV'][0]);q=float(f['charge_pe'][0]);source=float(f['source_energy_MeV'][0])
    assert e>0 and q>0 and np.isclose(source,record['energy_MeV'],rtol=2e-7,atol=1e-9)
    return dict(root=record['root'],seed=record['seed'],replicate=record['replicate'],
        source_energy_MeV=float(record['energy_MeV']),root_source_energy_MeV=source,edep_MeV=e,charge_pe=q,
        root_original_sha256=record['root_original_sha256'],truth_file=record['truth'],
        truth_sha256=digest(truth_path),metadata_file=record['metadata'],metadata_sha256=digest(meta_path),
        source_manifest=record['source_manifest'])


def collect():
    check_freeze();dest=OUTPUT/'calibration-records.json'
    if dest.exists():raise FileExistsError(dest)
    lookup={}
    for dataset in NATIVE_SOURCES:
        records=load(PROJECT/'jupyter/artifacts'/dataset/'native-records.json')['test_single']
        lookup.update({(dataset,r['seed']):r for r in records})
    config=load(EXPERIMENT/'sensor-config.json')
    cases=load(PROJECT/'jupyter/artifacts'/PARENT/'calibration.cases.json')
    with np.load(PROJECT/'jupyter/artifacts'/PARENT/'calibration-ancestry.truth.npz') as f:q=f['charge_pe']
    rows=[]
    for i,c in enumerate(cases):
        if c['kind']!='single':continue
        r=native_row(lookup[c['dataset'],c['seed']],config)
        assert np.isclose(r['charge_pe'],q[i,0],rtol=1e-10)
        assert r['replicate']//2%2==c['fold']
        rows.append(dict(r,fold=c['fold'],dataset=c['dataset'],calibration_index=i))
        if len(rows)%80==0:print('ENERGY_CALIBRATION_RECORDS',len(rows),flush=True)
    assert len(rows)==560 and len({r['root'] for r in rows})==560
    assert sum(r['fold']==0 for r in rows)==sum(r['fold']==1 for r in rows)==280
    save(dest,rows)
    save(RECORD/'calibration-provenance.json',dict(passed=True,rows=560,fold_counts=[280,280],
        charge_ancestry_checked=True,configs_identical=True,root_bytes_rehashed_this_round=False,
        root_identity='Reuse original ROOT hashes and verified export metadata; current derived truth bytes hashed',
        calibration_sha256=digest(dest)))
    print('ENERGY_CALIBRATION_COLLECTED',flush=True)


def material_table():
    text=(EXPERIMENT/'material-yield.inc').read_text()
    def vector(name):
        s=re.search(r'\b'+name+r'\s*=\s*\{([^}]*)\}',text,re.S).group(1)
        return np.array([float(v.strip().replace('*MeV','')) for v in s.split(',') if v.strip()])
    return vector('gaggScintillationEnergy'),vector('gaggAlphaScintillationYield')


def calibrate():
    check_freeze()
    if (EXPERIMENT/'energy-model.json').exists():raise FileExistsError('Energy model already frozen')
    rows=load(OUTPUT/'calibration-records.json')
    e=np.array([r['edep_MeV'] for r in rows]);q=np.array([r['charge_pe'] for r in rows]);g=np.array([r['source_energy_MeV'] for r in rows]);fold=np.array([r['fold'] for r in rows])
    te,tl=material_table();models={};results=[]
    for method in METHODS:
        try:
            model=fit_calibration(e[fold==0],q[fold==0],g[fold==0],te,tl,method)
            value,groups=selection_score(model,e[fold==1],q[fold==1],g[fold==1])
            models[method]=model;results.append(dict(method=method,selection_loss=value,groups=groups,valid=True))
        except (ValueError,RuntimeError) as exc:
            results.append(dict(method=method,valid=False,error=str(exc)))
    valid=[r for r in results if r['valid']]
    if not valid:raise RuntimeError('All prespecified energy mappings failed; preserve diagnostics')
    chosen=min(valid,key=lambda r:r['selection_loss'])['method']
    save(OUTPUT/'candidate-models.json',models)
    save(EXPERIMENT/'energy-model.json',models[chosen])
    save(EXPERIMENT/'mapping-selection.json',dict(created=now(),chosen=chosen,candidates=results,
        selection_inputs='Only fold1 independent single-event calibration, no evaluation energies',
        refit_on_all_calibration=False,calibration_sha256=digest(OUTPUT/'calibration-records.json')))
    save(EXPERIMENT/'energy-freeze.json',dict(created=now(),files={n:digest(EXPERIMENT/n) for n in
        ['energy-model.json','mapping-selection.json','freeze.json']},evaluation_energy_not_read=True))
    print('ENERGY_MAPPING_SELECTED',chosen,[(r['method'],r.get('selection_loss')) for r in results],flush=True)


def predict():
    check_freeze()
    for name,sha in load(EXPERIMENT/'energy-freeze.json')['files'].items():assert digest(EXPERIMENT/name)==sha
    path=OUTPUT/'energy.predictions.jsonl'
    if path.exists():raise FileExistsError(path)
    enabled=[True];opened=[]
    def guard(event,args):
        if enabled[0] and event=='open' and isinstance(args[0],(str,bytes)):
            name=str(args[0]);opened.append(name)
            if any(s in name for s in ['truth','labels','calibration-records']) or name.endswith('.root'):
                raise PermissionError('Energy inference cannot read truth/ROOT')
    sys.addaudithook(guard)
    for name in ['forbidden.truth.npz','forbidden.labels.json','forbidden.root']:
        try:open(name)
        except PermissionError:pass
        else:raise AssertionError('Truth guard did not fire')
    model=load(EXPERIMENT/'energy-model.json');result=[]
    for r in lines(OUTPUT/'fits.observations.jsonl'):
        energy=invert_charge(model,r['charge_pe'])
        result.append(dict(r,energy_MeV=[float(v) if np.isfinite(v) else None for v in energy]))
    write_lines(path,result);enabled[0]=False
    save(RECORD/'energy-inference.json',dict(passed=True,created=now(),entries=len(result),
        inputs=['fits.observations.jsonl','energy-model.json'],blocked_truth_probes=3,opened_paths=opened,
        predictions_sha256=digest(path),energy_model_sha256=digest(EXPERIMENT/'energy-model.json'),
        source_or_deposited_energy_inputs=False,count_time_charge_changed=False))
    print('ENERGY_PREDICTIONS_FROZEN',len(result),flush=True)


def summarize(rows):
    errors=[m['relative_energy_error'] for r in rows for m in r['matches'] if m['relative_energy_error'] is not None]
    target=sum(r['target'] for r in rows);n=len(rows)
    return dict(method=rows[0]['method'],kind=rows[0]['kind'],condition=rows[0]['condition'],gap_ns=rows[0]['gap_ns'],
        frames=n,target_events=target,matched_events=sum(r['matched'] for r in rows),
        missed=sum(r['missed'] for r in rows),extra=sum(r['extra'] for r in rows),
        split=sum(r.get('split',0) for r in rows),merge=sum(r.get('merge',0) for r in rows),
        out_of_domain_predictions=sum(r['energy_out_of_domain'] for r in rows),
        recovered_5=sum(r['strict_5'] for r in rows),recovered_10=sum(r['strict_10'] for r in rows),recovered_20=sum(r['strict_20'] for r in rows),
        energy_efficiency_10=sum(r['strict_10'] for r in rows)/target if target else None,
        frame_success_10=sum(r['frame_success_10'] for r in rows),
        frame_efficiency_10=sum(r['frame_success_10'] for r in rows)/n,
        extra_frames=sum(r['extra']>0 for r in rows),
        matched_energy_error_median=float(np.median(errors)) if errors else None,
        matched_energy_error_sigma68=float((np.quantile(errors,.84)-np.quantile(errors,.16))/2) if errors else None,
        q_strict_frame_success=sum(r.get('old_q_success',False) for r in rows))


def evaluate():
    check_freeze();assert load(RECORD/'energy-inference.json')['passed']
    assert digest(OUTPUT/'energy.predictions.jsonl')==load(RECORD/'energy-inference.json')['predictions_sha256']
    config=load(EXPERIMENT/'sensor-config.json');rows=[]
    native=load(PROJECT/'jupyter/artifacts'/CONFIRM/'native-records.json')['test_single']
    for r in native:
        rows.append(native_row(r,config))
        if len(rows)%80==0:print('ENERGY_EVALUATION_PARENTS',len(rows),flush=True)
    calibration=load(OUTPUT/'calibration-records.json')
    assert len(rows)==320 and not {r['root'] for r in calibration}&{r['root'] for r in rows}
    save(OUTPUT/'evaluation-parent-labels.json',rows)
    cases=load(PROJECT/'jupyter/artifacts'/CONFIRM/'labels.cases.json')
    with np.load(PROJECT/'jupyter/artifacts'/CONFIRM/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    e=np.zeros_like(q);block=np.full(len(cases),-1,int)
    for i,c in enumerate(cases):
        for slot,parent in enumerate(c['parents']):e[i,slot]=rows[parent]['edep_MeV']
        if c['kind']=='pair':
            blocks={rows[j]['replicate']//2 for j in c['parents']};assert len(blocks)==1
            block[i]=next(iter(blocks))
    np.savez_compressed(OUTPUT/'evaluation.labels.npz',energy_MeV=e,charge_pe=q,times_ns=t,block=block)
    save(OUTPUT/'evaluation.cases.json',cases)
    model=load(EXPERIMENT/'energy-model.json')
    old={(r['method'],r['index']):r for r in lines(PROJECT/'jupyter/artifacts'/CONFIRM/'scores.jsonl') if r['method'] in FIT_METHODS}
    predictions=lines(OUTPUT/'energy.predictions.jsonl');scored=[]
    for r in predictions:
        i=r['index'];c=cases[i];parent=old[r['method'],i]
        s=score_events(t[i],q[i],e[i],r['times_ns'],r['energy_MeV'],c['gap_ns'])
        assert s['matched']==parent['matched'] and s['missed']==parent['missed'] and s['extra']==parent['extra'],(r['method'],i,s,parent)
        assert r['times_ns']==parent['fit']['times_ns'] and r['charge_pe']==parent['fit']['charge_pe']
        scored.append(dict(s,index=i,method=r['method'],kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
                           block=int(block[i]),old_q_success=parent['success'],split=parent['split'],merge=parent['merge']))
    for i,c in enumerate(cases):
        mask=q[i]>=500.;energy=invert_charge(model,q[i,mask])
        s=score_events(t[i],q[i],e[i],t[i,mask],energy,c['gap_ns'])
        scored.append(dict(s,index=i,method='true_charge_diagnostic',kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],block=int(block[i])))
    write_lines(OUTPUT/'energy.scores.jsonl',scored)
    groups={}
    for r in scored:groups.setdefault((r['method'],r['kind'],r['condition'],r['gap_ns']),[]).append(r)
    summaries=[summarize(g) for g in groups.values()]
    # Parameter-map intervals: original independent-mother block design.
    rng=np.random.default_rng(load(EXPERIMENT/'protocol.json')['bootstrap_seed']);weights=np.zeros((2000,40),int)
    for start in (0,20):
        draws=rng.integers(start,start+20,size=(2000,20))
        for k in range(2000):weights[k]+=np.bincount(draws[k],minlength=40)
    aggregate={};bootstrap={}
    for method in FIT_METHODS+['true_charge_diagnostic']:
        values={True:[],False:[]};boots={True:[],False:[]}
        for key,rr in groups.items():
            if key[0]!=method or key[1]!='pair':continue
            ids=np.array([r['block'] for r in rr]);success=np.array([r['frame_success_10'] for r in rr],float)
            den=np.bincount(ids,minlength=40);num=np.bincount(ids,weights=success,minlength=40)
            weak='weak' in key[2];values[weak].append(float(success.mean()));boots[weak].append((weights@num)/(weights@den))
        aggregate[method]={name:float(np.mean(values[weak])) for weak,name in [(True,'weak'),(False,'strong')]}
        bootstrap[method]={name:np.mean(boots[weak],axis=0) for weak,name in [(True,'weak'),(False,'strong')]}
    paired={}
    for group in ['weak','strong']:
        diff=bootstrap[FIT_METHODS[1]][group]-bootstrap[FIT_METHODS[0]][group]
        paired[group]=dict(difference=aggregate[FIT_METHODS[1]][group]-aggregate[FIT_METHODS[0]][group],
                          mother_block95=np.quantile(diff,[.025,.975]).tolist())
    singles={}
    for group in sorted({r['source_energy_MeV'] for r in rows}):
        rr=[r for r in rows if r['source_energy_MeV']==group]
        singles[str(group)]=error_summary([r['edep_MeV'] for r in rr],invert_charge(model,[r['charge_pe'] for r in rr]))
    transfer_pass=all(r['valid']/r['n']>=.95 and r['median_relative'] is not None and abs(r['median_relative'])<=.05 for r in singles.values())
    save(OUTPUT/'report.json',dict(created=now(),chosen=load(EXPERIMENT/'mapping-selection.json')['chosen'],
        mapping_selection=load(EXPERIMENT/'mapping-selection.json')['candidates'],single_truth_charge=singles,
        single_mapping_transfer_pass=transfer_pass,summary=summaries,aggregate_strict_frame10=aggregate,paired_energy10=paired,
        domain='Native shared-cell simulated data; previously published Q/t sample, energy-calibration parents disjoint',
        independent_parent_counts=dict(calibration=280,selection=280,evaluation=320),frames=1600,
        no_refitting=True,no_gate_change=True,no_new_ROOT=True,default_changed=False,
        width_definition='Conditional relative-error sigma68=(p84-p16)/2, not measured spectral FWHM'))
    save(RECORD/'science-validation.json',dict(passed=True,created=now(),mothers_disjoint=True,parents=320,
        native_metadata_config_verified=True,truth_charge_primary_only=True,old_Q_time_detection_counts_replayed=3200,
        scored=4800,energy_predictions_frozen_before_evaluation=True,
        inputs={p.name:digest(p) for p in [OUTPUT/'calibration-records.json',OUTPUT/'evaluation-parent-labels.json',
            OUTPUT/'energy.predictions.jsonl',OUTPUT/'evaluation.labels.npz',OUTPUT/'energy.scores.jsonl']},
        report_sha256=digest(OUTPUT/'report.json')))
    print('ENERGY_EVALUATION_COMPLETE',dict(chosen=load(EXPERIMENT/'mapping-selection.json')['chosen'],
        transfer_pass=transfer_pass,aggregate=aggregate,paired=paired),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','collect','calibrate','predict','evaluate'])
    globals()[parser.parse_args().command]()


if __name__=='__main__':main()
