"""Sync report-only corrections; keep all computed notebook cells intact."""
from pathlib import Path
import hashlib
import json
import shutil

import nbformat as nbf

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
NAME = 'alpha_pair_recovery_energy_20261029'
EXP, OUT, REC = [ROOT/p/NAME for p in ['pileup/experiments', 'jupyter/artifacts', 'jupyter/run_records']]


def main():
    report = json.loads((OUT/'audited-report.json').read_text())
    assert report['aggregate']['static_nll']['not_converged'] is None
    old = json.loads((EXP/'result-summary.json').read_text())
    assert report['event_totals'] == old['event_totals'] and report['paired_comparisons'] == old['paired_comparisons']
    nbpath = ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb'
    nb = nbf.read(nbpath, as_version=4); assert len(nb.cells) == 55
    changed = []
    for i, cell in enumerate(nb.cells[42:], 42):
        if cell.cell_type != 'markdown': continue
        original = cell.source
        cell.source = cell.source.replace('| 原NLL＋静态能标 | 38.17% | 53.33% | 0 | 0 |',
            '| 原NLL＋静态能标 | 38.17% | 53.33% | 未复核 | 未复核 |')
        cell.source = cell.source.replace('run_report_v2.py', 'run_report_v3.py')
        cell.source = cell.source.replace('v2只修订报告审核并显式保留名义起点检查不通过，拟合输出完全未变。',
            'v2显式保留名义起点检查不通过，v3另将未保留的原基线优化器状态标为未复核，而非零失败。拟合输出完全未变。')
        if cell.source != original: changed.append(i)
    assert changed == [46, 53]
    nbf.validate(nb); nbf.write(nb, nbpath)
    shutil.copyfile(OUT/'audited-report.json', EXP/'result-summary.json')
    for path in OUT.glob('*.csv'): shutil.copyfile(path, EXP/'tables'/path.name)
    for name in ['report-validation.json', 'heldout-response-validation.json', 'plot-inputs.json']:
        shutil.copyfile(REC/name, EXP/name)
    path = REC/'notebook-validation.json'; status = json.loads(path.read_text())
    status['notebook_sha256'] = hashlib.sha256(nbpath.read_bytes()).hexdigest()
    status['report_metadata_correction_cells'] = changed
    path.write_text(json.dumps(status, indent=2)+'\n')
    print('REPORT_METADATA_SYNCED', changed, flush=True)


if __name__ == '__main__': main()
