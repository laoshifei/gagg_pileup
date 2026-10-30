"""Snapshot delivery scripts and hash the completed diagnostic artifacts."""
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json
import shutil

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
NAME = 'alpha_pair_energy_oracle_20261029'
EXP, OUT, REC = [ROOT/p/NAME for p in ['pileup/experiments', 'jupyter/artifacts', 'jupyter/run_records']]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')


def main():
    assert json.loads((REC/'notebook-validation.json').read_text())['passed']
    assert json.loads((REC/'report-validation.json').read_text())['passed']
    assert not (EXP/'delivery-freeze.json').exists()
    (EXP/'delivery').mkdir(exist_ok=True)
    for name in ['deliver.py', 'write_docs.py', 'validate_delivery.py', 'finalize.py']:
        shutil.copyfile(REC/name, EXP/'delivery'/name)
    for name in ['input-identity-validation.json', 'provenance.json', 'notebook-validation.json', 'figure-validation.json']:
        shutil.copyfile(REC/name, EXP/name)
    created = datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
    save(EXP/'implementation-state.json', dict(created=created,
        base_commit='67796208c81000563470044f2482e9622c582f77',
        study='Truth-assisted count/time diagnostics on previously exposed development mothers',
        source_state='New uncommitted analysis modules and experiment snapshot; prior uncommitted recovery-energy work preserved',
        notebooks_local_ignored=True, new_ROOTs=0, default_changed=False,
        committed=False, pushed=False))
    outputs = {p.name: sha(p) for p in sorted(OUT.iterdir()) if p.is_file()}
    checksums = dict(created=created, files=outputs,
        per_frame_predictions='Separately frozen in run_records blind-freeze.json and oracle-freeze.json',
        inputs='Separately frozen in experiment input-hashes.json')
    save(REC/'output-checksums.json', checksums)
    save(EXP/'output-checksums.json', checksums)
    files = {str(p.relative_to(EXP)): sha(p) for p in sorted(EXP.rglob('*'))
             if p.is_file() and p.name not in ['delivery-freeze.json', 'delivery-validation.json']
             and '__pycache__' not in p.parts}
    save(EXP/'delivery-freeze.json', dict(created=created, files=files,
        scope='Completed experiment snapshot, excluding this manifest and its subsequent validation',
        inherited_inference_freezes_unchanged=True))
    print('ORACLE_DELIVERY_FROZEN', len(files), 'experiment files;', len(outputs), 'result files', flush=True)


if __name__ == '__main__': main()
